#pragma once

#include "SimpleMesh.h"
#include "mathkit.h"
#include <vector>

class VoxelGrid {
public:
    VoxelGrid(glm::ivec3 size, aabb3 bounds);

    [[nodiscard]] size_t numFilledVoxels() const;

    [[nodiscard]] uint8_t get(int x, int y, int z) const;
    void set(int x, int y, int z, uint8_t);
    void set(vec3 point, uint8_t);

    void insertMesh(const SimpleMesh&);

    void writeToVox(const std::string& path) const;

private:
    aabb3 m_bounds;
    size_t m_sx, m_sy, m_sz;
    std::vector<uint8_t> m_grid;
};
