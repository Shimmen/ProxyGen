#pragma once

#include "SimpleMesh.h"
#include "Sphere.h"
#include "mathkit.h"
#include <functional>
#include <unordered_map>
#include <vector>

class VoxelGrid {
public:
    VoxelGrid(glm::ivec3 size, aabb3 bounds);

    ivec3 remapToGridSpace(vec3, float (*roundingFunc)(float)) const;
    vec3 voxelCenterPoint(ivec3) const;

    ivec3 gridDimensions() const;
    aabb3 gridBounds() const;
    vec3 voxelSize() const;

    struct TriangleRef {
        const SimpleMesh& mesh;
        uint64_t triangleIndex;
    };

    void forEachFilledVoxel(std::function<void(aabb3 aabb, const std::vector<TriangleRef>&)>) const;

    [[nodiscard]] size_t numFilledVoxels() const;
    [[nodiscard]] std::vector<vec3> filledVoxelsInSphere(const Sphere& sphere) const;

    [[nodiscard]] uint64_t linearIndex(ivec3 gridIndex) const;
    [[nodiscard]] uint64_t linearIndex(int x, int y, int z) const;

    [[nodiscard]] uint32_t get(int x, int y, int z) const;
    [[nodiscard]] uint32_t get(ivec3) const;
    void set(int x, int y, int z, uint32_t);
    void set(vec3 point, uint32_t);

    std::vector<ivec3> immediateFilledNeighbors(ivec3) const;

    void insertMesh(const SimpleMesh&, bool assignVoxelColorsToSurface);
    void fillVolumes();
    void subtractGrid(const VoxelGrid&);

    void quantizeColors(uint32_t numBins);
    void writeToVox(const std::string& path) const;

private:
    aabb3 m_bounds;
    size_t m_sx, m_sy, m_sz;

    std::vector<uint32_t> m_grid;
    std::vector<vec3> m_colors;

    std::unordered_map<uint64_t, std::vector<TriangleRef>> m_trianglesForVoxelIndex {};
};
