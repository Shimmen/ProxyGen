#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "Sphere.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <bobyqa.h>
#include <chrono>
#include <fmt/format.h>
#include <random>
#include <unordered_set>
#include <vector>

#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

struct SphereSet {
    SphereSet(const SimpleMesh& mesh)
        : mesh(mesh)
    {
    }

    std::vector<Sphere> spheres;

    const SimpleMesh& mesh;

    struct {
        std::unique_ptr<VoxelGrid> shell;
        std::unique_ptr<VoxelGrid> filled;
        std::unique_ptr<VoxelGrid> inside;
    } grids;
};

struct Triangle {
public:
    Triangle(vec3 v0, vec3 v1, vec3 v2)
    {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;

        vec3 crossP = cross(v1 - v0, v2 - v0);
        area = length(crossP) / 2.0f;
        normal = normalize(crossP);
    }

    vec3 v[3];
    vec3 normal;
    float area;
};

vec3 getRandomSphereLocationInsideMesh(const VoxelGrid& grid, std::default_random_engine& generator)
{
    aabb3 bounds = grid.gridBounds();
    std::uniform_real_distribution<float> uniformInGridX { bounds.min.x, bounds.max.x };
    std::uniform_real_distribution<float> uniformInGridY { bounds.min.y, bounds.max.y };
    std::uniform_real_distribution<float> uniformInGridZ { bounds.min.z, bounds.max.z };

    vec3 location;
    ivec3 gridCoord;

    do {
        location = { uniformInGridX(generator),
                     uniformInGridY(generator),
                     uniformInGridZ(generator) };
        gridCoord = grid.remapToGridSpace(location, std::round);
    } while (grid.get(gridCoord) == 0);

    return location;
}

void assignSpheresRandomlyInVolume(SphereSet& set, unsigned numSpheres)
{
    auto& grid = *set.grids.inside;
    assert(numSpheres <= grid.numFilledVoxels());

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator { seed };

    for (size_t i = 0; i < numSpheres; ++i) {
        vec3 location = getRandomSphereLocationInsideMesh(grid, generator);
        Sphere sphere { location, 0.0 };
        set.spheres.push_back(sphere);
    }

    aabb3 bounds = grid.gridBounds();
    std::uniform_real_distribution<float> uniformInGridX { bounds.min.x, bounds.max.x };
    std::uniform_real_distribution<float> uniformInGridY { bounds.min.y, bounds.max.y };
    std::uniform_real_distribution<float> uniformInGridZ { bounds.min.z, bounds.max.z };

    while (set.spheres.size() < numSpheres) {

        float x = uniformInGridX(generator);
        float y = uniformInGridY(generator);
        float z = uniformInGridZ(generator);

        ivec3 gridCoord = grid.remapToGridSpace({ x, y, z }, std::round);
        if (grid.get(gridCoord) > 0) {
            Sphere sphere { vec3(x, y, z), 0.0f };
            set.spheres.push_back(sphere);
        }
    }
}

bool sphereInsideMeshVolume(const Sphere& sphere, const SimpleMesh& mesh)
{
    // NOTE: Implicit, but this assumes that the mesh is one solid blob with no islands!

    // NOTE: To avoid expensive(/not-implemented) checks to see if the sphere centers are inside the volume
    //  of the object, we will assume they are and through other means verify that it will be invariant.

    if (sphere.radius == 0.0) {
        return true;
    }

    size_t triangleCount = mesh.triangleCount();
    for (size_t i = 0; i < triangleCount; ++i) {

        // TODO: Maybe don't reconstruct these every query?!
        vec3 v0, v1, v2;
        mesh.triangle(i, v0, v1, v2);
        Triangle triangle { v0, v1, v2 };

        double h = dot(triangle.normal, sphere.origin - v0);

        // Is the sphere center in front of the triangle? Since we know that sphere
        // centers always should be inside the volume, this simply means that the
        // mesh is convex and we can ignore this triangle for this sphere.
        if (h > 0.0) {
            continue;
        }

        // Is the sphere fully behind the triangle? Since its center also it we then
        // know this triangle isn't a problem here.
        if (h + sphere.radius < 0.0) {
            continue;
        }

        // We know that we are behind the triangle and that the sphere pokes through, so
        // the sphere is *not* inside the mesh volume.
        return false;
    }

    return true;
}

bool sphereSetInsideMeshVolume(const SphereSet& set, const SimpleMesh& mesh)
{
    for (const Sphere& sphere : set.spheres) {
        if (!sphereInsideMeshVolume(sphere, mesh)) {
            return false;
        }
    }
    return true;
}

void expandSpheresMaximally(SphereSet& set)
{
    aabb3 gridBounds = set.grids.shell->gridBounds();
    double gridDiameter = distance(gridBounds.min, gridBounds.max);

    for (Sphere& sphere : set.spheres) {

        Sphere minSphere { sphere.origin, 0.0 };
        Sphere maxSphere { sphere.origin, gridDiameter };

        assert(sphereInsideMeshVolume(minSphere, set.mesh));
        assert(!sphereInsideMeshVolume(maxSphere, set.mesh));

        while (abs(maxSphere.radius - minSphere.radius) > 1e-8) {

            double midpoint = (minSphere.radius + maxSphere.radius) / 2.0;
            Sphere middleSphere { sphere.origin, midpoint };

            bool midOk = sphereInsideMeshVolume(middleSphere, set.mesh);
            if (midOk) {
                minSphere = middleSphere;
            } else {
                maxSphere = middleSphere;
            }
        }

        // Make sure to be conservative with the choice to we don't break the invariant
        sphere.radius = minSphere.radius;
        assert(sphereInsideMeshVolume(sphere, set.mesh));
    }

#if 0
    std::sort(set.spheres.begin(), set.spheres.end(), [](const Sphere& lhs, const Sphere& rhs) {
        return lhs.radius > rhs.radius;
    });
    fmt::print("   grid 'radius' {}\n", gridDiameter / 2.0);
    for (auto& sphere : set.spheres) {
        fmt::print("     sphere radius {}\n", sphere.radius);
    }
#endif
}

double sphereSetIntersectionVolumeOverestimation(const SphereSet& set)
{
    // Except some special cases there are no analytical solutions for volume intersection of many spheres. See
    // https://www.researchgate.net/publication/289097700_Volume_of_intersection_of_six_spheres_A_special_case_of_practical_interest
    // for some more information. In this case we don't really care about the exact intersection, but we do want to avoid as much
    // intersection as possible, so just overcounting it works here, I think..

    double intersection = 0.0;

    size_t sphereCount = set.spheres.size();
    for (size_t i = 0; i < sphereCount; ++i) {
        const Sphere& self = set.spheres[i];
        for (size_t k = i + 1; k < sphereCount; ++k) {
            const Sphere& other = set.spheres[k];
            intersection += self.overlapWith(other);
        }
    }

    return intersection;
}

double sphereSetUnionVolume(const SphereSet& set)
{
    double volume = 0.0;

    size_t sphereCount = set.spheres.size();
    for (const Sphere& sphere : set.spheres) {
        volume += sphere.volume();
    }

    // TODO: Maybe scale down the intersection a bit, since we know it overestimates?
    double intersectionOverestimation = sphereSetIntersectionVolumeOverestimation(set);
    //assert(intersectionOverestimation < volume); (doesn't hold!)

    volume -= intersectionOverestimation;

    return volume;
}

struct SphereFittingData {
    SphereSet* set;
    unsigned sphereIndex;
};

REAL sphereFittingObjectiveFunction(const INTEGER n, const REAL* x, void* dataPtr)
{
    assert(n == 4);

    Sphere testSphere;
    testSphere.radius = x[3];
    testSphere.origin = { x[0], x[1], x[2] };

    SphereFittingData& data = *reinterpret_cast<SphereFittingData*>(dataPtr);
    SphereSet& set = *data.set;

    // Temporarily swap in the sphere to test
    Sphere currentSphere = set.spheres[data.sphereIndex];
    set.spheres[data.sphereIndex] = testSphere;

    if (!sphereInsideMeshVolume(testSphere, set.mesh)) {
        return 99999.99; // BOBYQA doesn't seem to like infinities, but this will do
    }

    // Since BOBYQA solves a minimization problem we want to minimize
    // the *negative* volume to get a maximal volume solution.
    double metric = -sphereSetUnionVolume(set);

    // Replace the current sphere back into the sphere set
    set.spheres[data.sphereIndex] = currentSphere;

    return metric;
}

void sphereFitting(SphereSet& set)
{
    // To avoid fitting the spheres in the same order every time.
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(set.spheres.begin(), set.spheres.end(), std::default_random_engine(seed));

    for (size_t sphereIdx = 0; sphereIdx < set.spheres.size(); ++sphereIdx) {

        Sphere& sphere = set.spheres[sphereIdx];
        assert(sphere.radius > 0.0);

        // Origin (3) + radius (1)
        constexpr INTEGER n = 3 + 1;

        // NPT is the number of interpolation conditions.  Its value must be in the interval
        // [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.
        // TODO: Find a good value!
#if 0
        constexpr INTEGER npt = n + 2; // (6)
#else
        constexpr INTEGER npt = (n + 1) * (n + 2) / 2; // (15)
#endif

        REAL x[n] = {
            sphere.origin.x,
            sphere.origin.y,
            sphere.origin.z,
            sphere.radius
        };

        // To inscribe a box = {min, max} in a sphere we want to find the min & max
        // points which lie along the diagonal that goes through the box center.
        //
        //  sqrt(x^2 + y^2 + y^3) = r
        //   since diagonal -> x=y=z = d
        //             sqrt(3d^2) = r
        //                   3d^2 = r^2
        //                      d = sqrt(r^2 / 3)
        //
        double d = std::sqrt(sphere.radius * sphere.radius / 3.0);
        vec3 boxDiagonal = vec3(d, d, d) * 0.99f; // (avoid overstepping!)
        aabb3 boxInSphere { sphere.origin - boxDiagonal, sphere.origin + boxDiagonal };

        REAL xLower[n] = {
            boxInSphere.min.x,
            boxInSphere.min.y,
            boxInSphere.min.z,
            1e-8
        };

        aabb3 gridBounds = set.grids.shell->gridBounds();
        double gridMaxDistance = length(gridBounds.max - gridBounds.min);

        REAL xUpper[n] = {
            boxInSphere.max.x,
            boxInSphere.max.y,
            boxInSphere.max.z,
            gridMaxDistance
        };

        // "Typically, RHOBEG should be about one tenth of the greatest expected change to a variable"
        // TODO: Find good value!
        REAL rhoBeg = (2.0 * d) / 10.0;

        // RHOEND should indicate the accuracy that is required in the final values of the variables
        // NOTE: Some models are tiny, others are huge.. So an absolute/constant accuracy is probably
        //  not very helpful here. A fraction of the sphere radius could maybe work as a scale-aware metric?
        // TODO: Find good value!
        REAL rhoEnd = gridMaxDistance * 1e-8;

        //The array W will be used for working space.  Its length must be at least
        // (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
        // will be set to the function value at the solution. */
        REAL workingMemory[(npt + 5) * (npt + n) + 3 * n * (n + 5) / 2];

        INTEGER logLevel = 0;
        INTEGER maxObjFunCalls = 10'000;

        SphereFittingData data { &set, sphereIdx };

        int status = bobyqa(n, npt, sphereFittingObjectiveFunction, (void*)(&data),
                            x, xLower, xUpper, rhoBeg, rhoEnd,
                            logLevel, maxObjFunCalls, workingMemory);

        switch (status) {
        case BOBYQA_SUCCESS: {
            // algorithm converged!
            vec3 newOrigin { x[0], x[1], x[2] };
            double newRadius = x[3];

#if 0
            fmt::print("      - ({}, {}, {}) r={}\n", newOrigin.x, newOrigin.y, newOrigin.z, newRadius);

            vec3 diffOrigin = newOrigin - sphere.origin;
            double diffRadius = newRadius - sphere.radius;
            //fmt::print("      - ({}, {}, {}) r={}\n", diffOrigin.x, diffOrigin.y, diffOrigin.z, diffRadius);
            //assert(length(diffOrigin) > 1e-6 && abs(diffRadius > 1e-6));
#endif

            sphere.origin = newOrigin;
            sphere.radius = newRadius;
            return;
        }
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

vec3 randomPointInSphere(const Sphere& sphere)
{
    std::uniform_real_distribution<double> randomOffset { -sphere.radius, +sphere.radius };

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator { seed };

    vec3 offset;
    do {
        offset = { randomOffset(generator), randomOffset(generator), randomOffset(generator) };
    } while (length(offset) > sphere.radius);

    return sphere.origin + offset;
}

void sphereTeleportation(SphereSet& set)
{
    if (set.spheres.size() == 1) {
        // TODO: Maybe pick a new starting point?
        // Or let the teleportation fail like this and restart from the main loop?
        return;
    }

    std::sort(set.spheres.begin(), set.spheres.end(), [](const Sphere& lhs, const Sphere& rhs) {
        return lhs.radius < rhs.radius;
    });

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator { seed };

    size_t numToReplace = std::ceil(set.spheres.size() / 2.0);
    for (size_t i = 0; i < numToReplace; ++i) {

        Sphere& sphere = set.spheres[i];
        sphere.origin = getRandomSphereLocationInsideMesh(*set.grids.inside, generator);

        // NOTE: A sphere expansion will happen after this is done, so we don't need to take care of that now
        sphere.radius = 0.0;
    }
}

void generateOutput(const SphereSet& sphereSet)
{
    size_t sphereCount = sphereSet.spheres.size();

    fmt::print("const vec4 spheres[] = vec4[](\n");
    for (size_t i = 0; i < sphereCount; ++i) {
        const Sphere& sphere = sphereSet.spheres[i];
        fmt::print("\tvec4({}, {}, {}, {}){}\n",
                   sphere.origin.x, sphere.origin.y, sphere.origin.z,
                   sphere.radius, (i + 1 == sphereCount) ? "" : ",");
    }
    fmt::print(");\n");
}

void setApplicationWorkingDirectory(char* executableName, const std::string& workingDir)
{
#ifdef WIN32
    char fullPathBuf[_MAX_PATH] = {};
    assert(_fullpath(fullPathBuf, executableName, sizeof(fullPathBuf)));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(_chdir(newWorkingDir.c_str()) == 0);
#else
    char fullPathBuf[PATH_MAX] = {};
    assert(realpath(executableName, fullPathBuf));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(chdir(newWorkingDir.c_str()) == 0);
#endif
}

int main(int argc, char** argv)
{
    setApplicationWorkingDirectory(argv[0], "ProxyGen");

    std::string path = "assets/Bunny/bunny_lowres.gltf";
    //std::string path = "assets/Avocado/Avocado.gltf";
    //std::string path = "assets/Sphere/sphere.gltf";
    //std::string path = "assets/BoomBox/BoomBoxWithAxes.gltf";

    constexpr unsigned numSpheres = 10;
    const size_t gridDimensions = 64;

    fmt::print("=> loading model '{}'\n", path);
    auto [basePath, model] = GltfUtil::loadModel(path);
    std::vector<SimpleMesh> simpleMeshes {};
    GltfUtil::bakeDownModelToSimpleMeshes(model, basePath, simpleMeshes);
    if (simpleMeshes.size() > 1) {
        fmt::print("==> model has more than one meshes; only the first will be used!\n");
    }

    const SimpleMesh& mesh = simpleMeshes[0];
    SphereSet set { mesh };

    fmt::print("=> creating voxel grids\n");

    aabb3 meshBounds { vec3(+INFINITY), vec3(-INFINITY) };
    mesh.extendAABB(meshBounds.min, meshBounds.max);

    set.grids.shell = std::make_unique<VoxelGrid>(glm::ivec3(gridDimensions), meshBounds);
    set.grids.shell->insertMesh(mesh, true);
    set.grids.shell->quantizeColors(10);

    set.grids.filled = std::make_unique<VoxelGrid>(*set.grids.shell);
    set.grids.filled->fillVolumes({ mesh });

    set.grids.inside = std::make_unique<VoxelGrid>(*set.grids.filled);
    set.grids.inside->subtractGrid(*set.grids.shell);

    fmt::print("=> assigning initial spheres\n");
    assignSpheresRandomlyInVolume(set, numSpheres);

    fmt::print("=> optimizing begin\n");

    bool didJustTeleport = false;
    double previousVolume = -std::numeric_limits<double>::infinity();

    constexpr int maxIterations = 1000;
    int iteration = 0;
    for (; iteration < maxIterations; ++iteration) {

        fmt::print("==> expanding spheres\n");
        expandSpheresMaximally(set);
        double newVolume = sphereSetUnionVolume(set);
        fmt::print("   volume: {}\n", newVolume);

        if (didJustTeleport) {
            didJustTeleport = false;
            if (newVolume <= previousVolume) {
                // TODO: Maybe consider calling assignSpheresRandomlyInVolume K times (i.e. restart) a few times before giving up!
                // If teleportation fails to increase volume, the algorithm terminates
                fmt::print("==> teleportation failed to increase volume\n");
                break;
            }
        }

        if (newVolume > previousVolume) {
            fmt::print("===> sphere fitting\n");
            sphereFitting(set);
        } else {
            fmt::print("===> sphere teleportation\n");
            sphereTeleportation(set);
            didJustTeleport = true;
        }

        previousVolume = newVolume;
    }

    if (iteration == maxIterations) {
        fmt::print("==> max iterations reached\n");
    }

    fmt::print("=> optimizing done\n");

    fmt::print("=> generating sphere output\n");
    generateOutput(set);
}
