#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "Sphere.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <bobyqa.h>
#include <direct.h>
#include <fmt/format.h>
#include <random>
#include <unordered_set>
#include <vector>

struct SphereSet {
    std::vector<Sphere> spheres;
    std::vector<float> currentSov;

    const VoxelGrid* filledGrid;
};

SphereSet preprocess(const VoxelGrid& filledGrid, int numSpheres)
{
    assert(numSpheres > 0);
    assert(numSpheres <= filledGrid.numFilledVoxels());

    SphereSet sphereSet {};
    sphereSet.filledGrid = &filledGrid;

    aabb3 bounds = filledGrid.gridBounds();
    std::uniform_real_distribution<float> uniformInGridX { bounds.min.x, bounds.max.x };
    std::uniform_real_distribution<float> uniformInGridY { bounds.min.y, bounds.max.y };
    std::uniform_real_distribution<float> uniformInGridZ { bounds.min.z, bounds.max.z };
    std::default_random_engine generator {};

    while (sphereSet.spheres.size() < numSpheres) {

        float x = uniformInGridX(generator);
        float y = uniformInGridY(generator);
        float z = uniformInGridZ(generator);

        ivec3 gridCoord = filledGrid.remapToGridSpace({ x, y, z }, std::round);
        if (filledGrid.get(gridCoord) > 0) {
            Sphere sphere { vec3(x, y, z), 0.0f };
            sphereSet.spheres.push_back(sphere);
            sphereSet.currentSov.push_back(0.0f);
        }
    }

    return sphereSet;
}

float sphereOutsideVolumeUsingVoxels(const Sphere& sphere, const VoxelGrid& grid)
{
    // TODO: Is there enough precision when using voxels?

    if (sphere.radius < 1e-6f) {
        return 0.0f;
    }

    // TODO: It's (very?) important that we don't have holes in the filled grid now..
    int numOutsideVoxels = 0;

    vec3 voxelSize = grid.voxelSize();
    vec3 first = sphere.origin - vec3(sphere.radius) - (0.05f * voxelSize);
    vec3 last = sphere.origin + vec3(sphere.radius) + (0.05f * voxelSize);

    first = max(first, grid.gridBounds().min);
    last = min(last, grid.gridBounds().max);

    for (float z = first.z; z <= last.z; z += voxelSize.z) { // NOLINT(cert-flp30-c)
        for (float y = first.y; y <= last.y; y += voxelSize.y) { // NOLINT(cert-flp30-c)
            for (float x = first.x; x <= last.x; x += voxelSize.x) { // NOLINT(cert-flp30-c)

                vec3 samplePoint = vec3(x, y, z);

                // Only count voxels that actually are inside the sphere
                if (distance2(sphere.origin, samplePoint) > sphere.radius * sphere.radius) {
                    continue;
                }

                ivec3 gridSampleCoord = grid.remapToGridSpace(samplePoint, std::round);

                // NOTE: Should be handled by the coord clamping above this loop
                //if (any(lessThan(gridSampleCoord, ivec3(0))) || any(greaterThanEqual(gridSampleCoord, grid.gridDimensions()))) {
                //    continue;
                //}

                // If this voxel (which is considered to be inside the sphere) is empty, we assume that the voxel
                // is actually outside the mesh (since the grid should be filled), and so we count this voxel.
                uint32_t value = grid.get(gridSampleCoord);
                if (value == 0) {
                    numOutsideVoxels += 1;
                }
            }
        }
    }

    float sov = float(numOutsideVoxels) * (voxelSize.x * voxelSize.y * voxelSize.z);
    return sov;
}

float totalSphereOutsideVolumeUsingVoxels(const SphereSet& set, const VoxelGrid& grid)
{
    float totalSOV = 0.0f;
    for (auto& sphere : set.spheres) {
        totalSOV += sphereOutsideVolumeUsingVoxels(sphere, grid);
    }
    return totalSOV;
}

float pointAssignment(SphereSet& set, float previousError)
{
    // TODO: From the paper "the [current] (naive) algorithm for point assignment can be accelerated"
    //  See p. 5 in the paper for more information. It's already slow and the SOV function is just a stub!

    const VoxelGrid& grid = *set.filledGrid;
    size_t numFilledVoxels = grid.numFilledVoxels();

    struct PointRef {
        PointRef(vec3 p, ivec3 c)
            : point(p)
            , gridIndex(c)
        {
        }
        vec3 point;
        ivec3 gridIndex;
    };

    std::unordered_set<uint64_t> visitedPoints {};
    visitedPoints.reserve(numFilledVoxels);

    std::vector<PointRef> stack {};

    for (auto& [sphereCenter, _] : set.spheres) {
        ivec3 gridIndex = grid.remapToGridSpace(sphereCenter, std::round);
        PointRef pointRef { sphereCenter, gridIndex };
        stack.push_back(pointRef);

        uint64_t linearIndex = grid.linearIndex(gridIndex);
        visitedPoints.insert(linearIndex);
    }

    auto oldSpheres = set.spheres;
    auto oldSphereSov = set.currentSov;

    int numAssigned = 0;

    // TODO: OpenMP would be nice here! Just make sure to mark some critical sections! However, it seems like the rounding of
    //  coordinates (std::roundf) is the main contributor of self-time according to the profiler. Worth looking into omp anyway!
    while (!stack.empty()) {
        PointRef pointRef = stack.back();
        stack.pop_back();

        int bestSphereIdx = -1;
        float smallestSovIncrease = INFINITY;

        float newPotentialRadius = INFINITY;
        float newPotentialSov = INFINITY;

        for (int si = 0; si < set.spheres.size(); ++si) {
            Sphere& sphere = set.spheres[si];

            // Create a test sphere with radius required to contain the test point
            Sphere testSphere = sphere;
            testSphere.radius = distance(sphere.origin, pointRef.point);

            // Check the increase in SOV for the test sphere
            float testSov = sphereOutsideVolumeUsingVoxels(testSphere, grid);
            float sovIncrease = testSov - set.currentSov[si];

            if (sovIncrease < smallestSovIncrease) {
                smallestSovIncrease = sovIncrease;
                bestSphereIdx = si;

                // Save these for speed (radius) and accuracy (SOV)
                newPotentialRadius = testSphere.radius;
                newPotentialSov = testSov;
            }
        }

        // Expand the best sphere to contain the new point
        assert(bestSphereIdx != -1);
        set.spheres[bestSphereIdx].radius = newPotentialRadius;
        set.currentSov[bestSphereIdx] = newPotentialSov;

        // Insert unvisited neighbors into the stack
        const auto& immediateNeighbors = grid.immediateFilledNeighbors(pointRef.gridIndex);
        for (ivec3 gridCoord : immediateNeighbors) {

            uint64_t linearIndexOfNeighbor = grid.linearIndex(gridCoord);
            if (visitedPoints.find(linearIndexOfNeighbor) != visitedPoints.end()) {
                continue;
            }

            vec3 voxelCenter = grid.voxelCenterPoint(gridCoord);
            stack.emplace_back(voxelCenter, gridCoord);

            uint64_t linearIndex = grid.linearIndex(gridCoord);
            visitedPoints.insert(linearIndex);

            assert(stack.size() <= numFilledVoxels);
        }

        numAssigned += 1;

        if (numAssigned % 100 == 0 || numAssigned == numFilledVoxels) {
            fmt::print("     assigned {}/{} points ({:.1f}%)\r", numAssigned, numFilledVoxels, 100.0f * float(numAssigned) / float(numFilledVoxels));
            fflush(stdout);
        }
    }
    fmt::print("\n");

    float newError = totalSphereOutsideVolumeUsingVoxels(set, *set.filledGrid);

    if (newError < previousError) {
        // Success, continue
        return newError;
    } else {
        // Failure, rollback the changes & trigger sphere teleportation
        set.spheres = oldSpheres;
        set.currentSov = oldSphereSov;
        return previousError;
    }
}

REAL sphereFittingObjectiveFunction(const INTEGER n, const REAL* x, void* data)
{
    // TODO: Implement the real objective function here!

    INTEGER i;
    REAL f, temp, tempa, tempb;

    f = 0.0;
    for (i = 4; i <= n; i += 2) {
        INTEGER j, j1 = i - 2;
        for (j = 2; j <= j1; j += 2) {
            tempa = x[i - 2] - x[j - 2];
            tempb = x[i - 1] - x[j - 1];
            temp = tempa * tempa + tempb * tempb;
            temp = std::max(temp, 1e-6);
            f += 1.0 / std::sqrt(temp);
        }
    }
    return f;
}

void sphereFitting(SphereSet& set)
{
    for (auto& sphere : set.spheres) {

        if (sphere.radius <= 0.0f) {
            continue;
        }

        // Origin (3) + radius (1)
        constexpr INTEGER n = 3 + 1;

        // NPT is the number of interpolation conditions.  Its value must be in the interval
        // [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.
        constexpr INTEGER npt = n + 2; // TODO: Find a good value!

        REAL x[n] = {
            sphere.origin.x,
            sphere.origin.y,
            sphere.origin.z,
            sphere.radius
        };

        // TODO: Find good values!
        REAL xLower[n] = {
            x[0] - sphere.radius,
            x[1] - sphere.radius,
            x[2] - sphere.radius,
            0.0
        };

        // TODO: Find good values!
        REAL xUpper[n] = {
            x[0] + sphere.radius,
            x[1] + sphere.radius,
            x[2] + sphere.radius,
            // (definitely shouldn't need to go higher than 10x current radius)
            sphere.radius * 10.0
        };

        // "Typically, RHOBEG should be about one tenth of the greatest expected change to a variable"
        // (and we probably should never need to move "outside" the current sphere with either the origin point or the radius.)
        // TODO: Find good value!
        REAL rhoBeg = sphere.radius / 10.0;

        // RHOEND should indicate the accuracy that is required in the final values of the variables
        // NOTE: Some models are tiny, others are huge.. So an absolute/constant accuracy is probably
        //  not very helpful here. A fraction of the sphere radius could maybe work as a scale-aware metric?
        // TODO: Find good value!
        REAL rhoEnd = sphere.radius / 100.0;

        //The array W will be used for working space.  Its length must be at least
        // (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
        // will be set to the function value at the solution. */
        REAL workingMemory[(npt + 5) * (npt + n) + 3 * n * (n + 5) / 2];

        INTEGER logLevel = 2;
        INTEGER maxObjFunCalls = 100;

        int status = bobyqa(n, npt, sphereFittingObjectiveFunction, nullptr,
                            x, xLower, xUpper, rhoBeg, rhoEnd,
                            logLevel, maxObjFunCalls, workingMemory);

        switch (status) {
        case BOBYQA_SUCCESS:
            // algorithm converged!
            return;
        case BOBYQA_BAD_NPT:
            fmt::print("bobyqa error: NPT is not in the required interval\n");
            assert(false); // npt is hardcoded to a valid value!
            break;
        case BOBYQA_TOO_CLOSE:
            fmt::print("bobyqa error: insufficient space between the bounds\n");
            break;
        case BOBYQA_ROUNDING_ERRORS:
            fmt::print("bobyqa error: too much cancellation in a denominator\n");
            break;
        case BOBYQA_TOO_MANY_EVALUATIONS:
            fmt::print("bobyqa error: maximum number of function evaluations exceeded\n");
            break;
        case BOBYQA_STEP_FAILED:
            fmt::print("bobyqa error: a trust region step has failed to reduce Q\n");
            break;
        default:
            assert(false);
        }
    }
}

float sphereOverlap(const Sphere& a, const Sphere& b)
{
    // See https://math.stackexchange.com/questions/297751/overlapping-spheres

    float d = distance(a.origin, b.origin);
    if (d >= a.radius + b.radius) {
        return 0.0f;
    }

    const float& r1 = a.radius;
    const float& r2 = b.radius;

    using namespace mathkit;
    return PI / (12.0f * d) * square(r1 + r2 - d) * (square(d) + 2.0f * d * (r1 + r2) - 3.0f * square(r1 - r2));
}

void sphereTeleportation(SphereSet& set)
{
    //
    // Find the sphere which is most overlapped by the other spheres
    //  & the sphere with highest SOV value
    //

    float totalVolume = 0.0f;
    for (auto& sphere : set.spheres) {
        totalVolume += sphere.volume();
    }

    float maxOverlapRatio = -INFINITY;
    int maxOverlapRatioIdx = -1;

    float maxSov = -INFINITY;
    int maxSovIdx = -1;

    for (size_t i = 0; i < set.spheres.size(); ++i) {

        Sphere& self = set.spheres[i];

        float sharedVolume = 0.0f;

        for (size_t k = 0; k < set.spheres.size(); ++k) {
            if (i == k) {
                continue;
            }

            Sphere& other = set.spheres[k];
            sharedVolume += sphereOverlap(self, other);
        }

        float overlapRatio = sharedVolume / totalVolume; // TODO: Correct interpretation? See p. 5 in [Wang 1].
            //  Don't think so, because the total volume here is constant!
        if (overlapRatio > maxOverlapRatio) {
            maxOverlapRatio = overlapRatio;
            maxOverlapRatioIdx = i;
        }

        float sov = sphereOutsideVolumeUsingVoxels(self, *set.filledGrid);
        if (sov > maxSov) {
            maxSov = sov;
            maxSovIdx = i;
        }
    }

    assert(maxSovIdx != -1);
    assert(maxOverlapRatioIdx != -1);

    //
    // Split the highest SOV sphere into two, with new sphere centers being
    // the two points furthest away from each other in the initial sphere,
    // and remove the most overlapped sphere
    //

    Sphere& sphereToSplit = set.spheres[maxSovIdx];
    Sphere& sphereToRemove = set.spheres[maxOverlapRatioIdx];

    const auto& splitCandidates = set.filledGrid->filledVoxelsInSphere(sphereToSplit);
    assert(!splitCandidates.empty());

    if (splitCandidates.size() == 1) {
        // Only one voxel, so pick two two corners of the voxel
        sphereToSplit.origin = splitCandidates[0] - 0.5f * set.filledGrid->voxelSize();
        sphereToRemove.origin = splitCandidates[0] + 0.5f * set.filledGrid->voxelSize();
    } else {
        float maxDistance2 = 0.0f;
        for (size_t i = 0, len = splitCandidates.size(); i < len; ++i) {
            const vec3& a = splitCandidates[i];
            for (size_t k = i + 1; k < len; ++k) {
                const vec3& b = splitCandidates[k];
                float dist2 = distance2(a, b);
                if (dist2 > maxDistance2) {
                    maxDistance2 = dist2;
                    sphereToSplit.origin = a;
                    sphereToRemove.origin = b;
                }
            }
        }

        assert(maxDistance2 > 1e-6f);
    }
}

void setApplicationWorkingDirectory(char* executableName, const std::string& workingDir)
{
    char fullPathBuf[_MAX_PATH] = {};
    assert(_fullpath(fullPathBuf, executableName, sizeof(fullPathBuf)));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(_chdir(newWorkingDir.c_str()) == 0);
}

int main(int argc, char** argv)
{
    setApplicationWorkingDirectory(argv[0], "ProxyGen");

    std::string path = "assets/Bunny/bunny.gltf";
    //std::string path = "assets/Avocado/Avocado.gltf";
    //std::string path = "assets/BoomBox/BoomBoxWithAxes.gltf";

    constexpr int numSpheres = 16;
    const size_t gridDimensions = 32;

    std::vector<SimpleMesh> simpleMeshes {};

    fmt::print("=> loading model '{}'\n", path);

    auto [basePath, model] = GltfUtil::loadModel(path);
    GltfUtil::bakeDownModelToSimpleMeshes(model, basePath, simpleMeshes);

    fmt::print("=> shell voxelization\n");

    aabb3 boundsOfAllMeshes = SimpleMesh::calculateBounds(simpleMeshes);
    VoxelGrid shellGrid { glm::ivec3(gridDimensions), boundsOfAllMeshes };

    for (size_t i = 0; i < simpleMeshes.size(); ++i) {
        fmt::print("     mesh {}/{}\n", i + 1, simpleMeshes.size());
        shellGrid.insertMesh(simpleMeshes[i], true);
    }

    fmt::print("=> color quantization\n");

    shellGrid.quantizeColors(10);

    fmt::print("=> volume filling\n");

    VoxelGrid filledGrid = shellGrid;
    filledGrid.fillVolumes(simpleMeshes);

    fmt::print("=> preprocessing sphere set\n");
    SphereSet sphereSet = preprocess(filledGrid, numSpheres);

    fmt::print("=> optimizing begin\n");

    bool didJustTeleport = false;
    float previousError = INFINITY;

    constexpr int maxIterations = 100;
    int iteration = 0;
    for (; iteration < maxIterations; ++iteration) {

        fmt::print("==> point assignment\n");
        float newError = pointAssignment(sphereSet, previousError);
        fmt::print("     error: {}\n", newError);

        if (newError < 1e-12f) {
            // in practice this probably shouldn't happen, but since our voxel grid resolution is quite low it can happen.
            fmt::print("==> error is zero\n");
            break;
        }

        if (didJustTeleport) {
            didJustTeleport = false;
            if (newError >= previousError) {
                // "if [teleportation] fails to decrease error, the algorithm terminates"
                fmt::print("==> teleportation failed to decrease error\n");
                break;
            }
        }

        if (newError < previousError) {
            fmt::print("===> sphere fitting\n");
            sphereFitting(sphereSet);
        } else {
            fmt::print("===> sphere teleportation\n");
            sphereTeleportation(sphereSet);
            didJustTeleport = true;
        }

        previousError = newError;
    }

    if (iteration == maxIterations) {
        fmt::print("==> max iterations reached\n");
    }

    fmt::print("=> optimizing done\n");

    fmt::print("=> generating sphere output\n");

    fmt::print("const vec4 spheres[] = vec4[](\n");
    for (size_t i = 0; i < numSpheres; ++i) {
        const Sphere& sphere = sphereSet.spheres[i];
        fmt::print("\tvec4({}, {}, {}, {}){}\n",
                   sphere.origin.x, sphere.origin.y, sphere.origin.z,
                   sphere.radius, (i + 1 == numSpheres) ? "" : ",");
    }
    fmt::print(");\n");
}
