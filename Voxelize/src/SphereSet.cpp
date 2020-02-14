#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <fmt/format.h>
#include <random>
#include <unordered_set>
#include <vector>

struct Sphere {
    vec3 origin;
    float radius;
};

struct SphereSet {
    std::vector<Sphere> spheres;
    std::vector<float> currentSov;
    float currentTotalSov;

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

enum class PointAssignmentResult {
    PerformSphereFitting,
    PerformSphereTeleportation,
};

PointAssignmentResult pointAssignment(SphereSet& set)
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
    set.currentTotalSov = totalSphereOutsideVolumeUsingVoxels(set, grid);

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
            fmt::print("     assigned {}/{} points ({:.1f}%)\n", numAssigned, numFilledVoxels, 100.0f * float(numAssigned) / float(numFilledVoxels));
        }
    }

    float newSov = totalSphereOutsideVolumeUsingVoxels(set, *set.filledGrid);

    if (newSov < set.currentTotalSov) {
        // Success, update the total SOV and continue
        set.currentTotalSov = newSov;
        return PointAssignmentResult::PerformSphereFitting;
    } else {
        // Failure, rollback the changes & trigger sphere teleportation
        set.spheres = oldSpheres;
        set.currentSov = oldSphereSov;
        return PointAssignmentResult::PerformSphereTeleportation;
    }
}

void sphereFitting(SphereSet& set)
{
    fmt::print("TODO: Implement sphere fitting!\n");
    // TODO!
}

void sphereTeleportation(SphereSet& set)
{
    fmt::print("TODO: Implement sphere teleportation!\n");
    // TODO!
}

int main()
{
    std::string path = "../assets/Avocado/Avocado.gltf";
    //std::string path = "../assets/BoomBox/BoomBoxWithAxes.gltf";

    const size_t gridDimensions = 16;
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
    SphereSet sphereSet = preprocess(filledGrid, 32);

    fmt::print("=> TEST: initial point assignment\n");
    auto result = pointAssignment(sphereSet);

    switch (result) {
    case PointAssignmentResult::PerformSphereFitting:
        fmt::print("==> TEST: sphere fitting\n");
        sphereFitting(sphereSet);
        break;
    case PointAssignmentResult::PerformSphereTeleportation:
        fmt::print("==> TEST: sphere teleportation\n");
        sphereTeleportation(sphereSet);
        break;
    }
}
