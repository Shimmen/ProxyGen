#pragma once

#include "SimpleMesh.h"
#include "mathkit.h"
#include <unordered_map>
#include <vector>

class VoxelGrid {
public:
    static constexpr uint32_t Empty = 0u;

    VoxelGrid(glm::ivec3 size, aabb3 bounds);

    ivec3 remapToGridSpace(vec3, float (*roundingFunc)(float)) const;

    [[nodiscard]] size_t numFilledVoxels() const;

    [[nodiscard]] uint64_t linearIndex(int x, int y, int z) const;

    [[nodiscard]] uint32_t get(int x, int y, int z) const;
    void set(int x, int y, int z, uint32_t);
    void set(vec3 point, uint32_t);

    void insertMesh(const SimpleMesh&, bool assignVoxelColorsToSurface);
    void fillVolumes(const std::vector<SimpleMesh>&);

    void quantizeColors(uint32_t numBins);
    void writeToVox(const std::string& path) const;

private:
    aabb3 m_bounds;
    size_t m_sx, m_sy, m_sz;

    std::vector<uint32_t> m_grid;
    std::vector<vec3> m_colors;

    struct TriangleRef {
        const SimpleMesh* mesh;
        uint64_t triangleIndex;
    };
    std::unordered_map<uint64_t, std::vector<TriangleRef>> m_trianglesForVoxelIndex {};
};
