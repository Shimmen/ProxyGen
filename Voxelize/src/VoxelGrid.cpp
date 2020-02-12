#include "VoxelGrid.h"

#include "tribox3.h"
#include <fmt/format.h>
#include <fstream>
#include <sstream>

VoxelGrid::VoxelGrid(glm::ivec3 size, aabb3 bounds)
    : m_bounds(bounds)
    , m_sx(size.x)
    , m_sy(size.y)
    , m_sz(size.z)
{
    size_t totalSize = size.x * size.y * size.z;
    assert(totalSize < 4ul * 1024ul * 1024ul * 1024ul);

    // note: this will allocate the whole dense grid
    m_grid.resize(totalSize);
}

ivec3 VoxelGrid::remapToGridSpace(vec3 point, float (*roundingFunc)(float)) const
{
    vec3 normalized = (point - m_bounds.min) / (m_bounds.max - m_bounds.min);
    vec3 floatGridSpace = normalized * vec3(m_sx - 1, m_sy - 1, m_sz - 1);

    int x = (int)roundingFunc(floatGridSpace.x);
    int y = (int)roundingFunc(floatGridSpace.y);
    int z = (int)roundingFunc(floatGridSpace.z);
    return { x, y, z };
}

size_t VoxelGrid::numFilledVoxels() const
{
    size_t count = 0;

    size_t totalSize = m_sx * m_sy * m_sz;
    for (size_t i = 0; i < totalSize; ++i) {
        if (m_grid[i] > 0) {
            count += 1;
        }
    }

    return count;
}

uint64_t VoxelGrid::linearIndex(int x, int y, int z) const
{
    assert(x >= 0 && x < m_sx && y >= 0 && y < m_sy && z >= 0 && z < m_sz);
    return x + y * (m_sx) + z * (m_sx * m_sy);
}

uint32_t VoxelGrid::get(int x, int y, int z) const
{
    return m_grid[linearIndex(x, y, z)];
}

void VoxelGrid::set(int x, int y, int z, uint32_t value)
{
    m_grid[linearIndex(x, y, z)] = value;
}

void VoxelGrid::set(vec3 point, uint32_t value)
{
    ivec3 grid = remapToGridSpace(point, std::round);
    set(grid.x, grid.y, grid.z, value);
}

void VoxelGrid::insertMesh(const SimpleMesh& mesh, bool assignVoxelColorsToSurface)
{
    vec3 voxelSize = (m_bounds.max - m_bounds.min) / vec3(m_sx, m_sy, m_sz);
    vec3 voxelHalfSize = 0.5f * voxelSize;

    vec3 centerOfFirst = m_bounds.min + voxelHalfSize;

    size_t numTris = mesh.triangleCount();
    for (size_t ti = 0; ti < numTris; ++ti) {

        vec3 v[3];
        mesh.triangle(ti, v[0], v[1], v[2]);

        vec3 triangleMin = min(v[0], min(v[1], v[2]));
        vec3 triangleMax = max(v[0], max(v[1], v[2]));

        ivec3 gridFirst = remapToGridSpace(triangleMin, std::floor);
        ivec3 gridLast = remapToGridSpace(triangleMax, std::ceil);

        for (size_t z = gridFirst.z; z <= gridLast.z; ++z) {
            for (size_t y = gridFirst.y; y <= gridLast.y; ++y) {
                for (size_t x = gridFirst.x; x <= gridLast.x; ++x) {

                    vec3 voxelCenter = centerOfFirst + (vec3(x, y, z) * voxelSize);

                    if (triangleBoxIntersection(voxelCenter, voxelHalfSize, v)) {
                        uint32_t value = 1;
                        if (assignVoxelColorsToSurface) {

                            // Compute barycentric coordinates of point projected onto the triangle
                            // https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle

                            vec3 e1 = v[1] - v[0];
                            vec3 e2 = v[2] - v[0];

                            vec3 normal = normalize(cross(e1, e2));
                            vec3 w = voxelCenter - v[0];

                            float gamma = dot(normal, cross(e1, w));
                            float beta = dot(normal, cross(w, e2));
                            float alpha = 1.0f - gamma - beta;

                            if (alpha >= 0.0f && alpha <= 1.0f && beta >= 0.0f && beta <= 1.0f && gamma >= 0.0f && gamma <= 1.0f) {
                                vec2 uv0, uv1, uv2;
                                mesh.triangleTexcoords(ti, uv0, uv1, uv2);

                                vec2 uv = alpha * uv0 + beta * uv1 + gamma * uv2;
                                vec3 color = mesh.texture().sample(uv);

                                value = m_colors.size();
                                m_colors.push_back(color);

                                fmt::print("got color {}, {}, {}\n", color.r, color.g, color.b);
                            } else {
                                fmt::print("no color luck\n");
                            }
                        }

                        set(x, y, z, value);

                        uint64_t index = linearIndex(x, y, z);
                        m_trianglesForVoxelIndex[index].push_back({ &mesh, ti });
                    }
                }
            }
        }
    }

    fmt::print("Num colors: {}\n", m_colors.size());
}

void VoxelGrid::fillVolumes(const std::vector<SimpleMesh>& meshes)
{
    const vec3 stepDirection = { 1, 0, 0 };

    for (size_t z = 0; z < m_sz; ++z) {
        for (size_t y = 0; y < m_sy; ++y) {

            bool filling = false;

            for (size_t x = 0; x < m_sx; ++x) {

                uint8_t value = get(x, y, z);
                if (value == 0) {
                    if (filling) {
                        set(x, y, z, 1);
                    }
                    continue;
                }

                vec3 normal { 0.0f };

                auto& triangleRefs = m_trianglesForVoxelIndex[linearIndex(x, y, z)];
                for (auto& [mesh, triangleIdx] : triangleRefs) {

                    vec3 v0, v1, v2;
                    mesh->triangle(triangleIdx, v0, v1, v2);

                    vec3 triNormal = normalize(cross(v1 - v0, v2 - v0));
                    normal += triNormal;
                }

                normal = normalize(normal);
                float d = dot(normal, stepDirection);

                // Be liberal when it comes to stopping the filling ...
                if (filling && d > 0.0f) {
                    filling = false;
                    continue;
                }

                // ... but be conservative when it comes to starting it.
                // This avoids errors where we don't stop and fill entire lines along the x-axis
                // which is quite egregious. On the other hand, under filling is probably okay.
                if (!filling && d < -0.4f) {
                    filling = true;
                    continue;
                }
            }
        }
    }
}

void VoxelGrid::writeToVox(const std::string& path) const
{
    if (m_sx > 126 || m_sy > 126 || m_sz > 126) {
        fmt::print("VoxelGrid::writeToVox(): the .vox format only supports resultions up to 126x126x126, so parts of the voxel grid might be cut off\n");
    }

    auto writeLabel = [](std::ostream& stream, const char* label) {
        stream.write(label, strlen(label));
    };
    auto writeUInt32 = [](std::ostream& stream, uint32_t value) {
        stream.write((const char*)&value, 4);
    };
    auto writeUInt8 = [](std::ostream& stream, uint8_t value) {
        stream.write((const char*)&value, 1);
    };
    auto writeChunkHeader = [&](std::ostream& stream, const char* chunkId, uint32_t sizeOfSelf, uint32_t sizeOfChildren) {
        writeLabel(stream, chunkId);
        writeUInt32(stream, sizeOfSelf);
        writeUInt32(stream, sizeOfChildren);
    };
    auto writeBuffer = [&](std::ostream& stream, std::ostringstream& source) {
        std::string content = source.str();
        stream << content;
    };

    uint32_t sx = std::min(m_sx, 126ul);
    uint32_t sy = std::min(m_sy, 126ul);
    uint32_t sz = std::min(m_sz, 126ul);

    std::ostringstream sizeBuffer {};
    writeUInt32(sizeBuffer, sx);
    writeUInt32(sizeBuffer, sy);
    writeUInt32(sizeBuffer, sz);

    std::ostringstream xyziBuffer {};
    writeUInt32(xyziBuffer, numFilledVoxels());

    for (size_t z = 0; z < m_sz; ++z) {
        for (size_t y = 0; y < m_sy; ++y) {
            for (size_t x = 0; x < m_sx; ++x) {
                uint8_t voxel = get(x, y, z);
                if (voxel > 0) {
                    if (x > 126 || y > 126 || z > 126) {
                        continue;
                    }
                    writeUInt8(xyziBuffer, sx - x - 1);
                    writeUInt8(xyziBuffer, z);
                    writeUInt8(xyziBuffer, y);
                    writeUInt8(xyziBuffer, voxel);
                }
            }
        }
    }

    std::ostringstream mainBuffer {};
    writeChunkHeader(mainBuffer, "SIZE", sizeBuffer.tellp(), 0);
    writeBuffer(mainBuffer, sizeBuffer);
    writeChunkHeader(mainBuffer, "XYZI", xyziBuffer.tellp(), 0);
    writeBuffer(mainBuffer, xyziBuffer);

    std::ofstream file;
    file.open(path, std::ios::binary | std::ios::out);

    writeLabel(file, "VOX ");
    writeUInt32(file, 150u);

    writeChunkHeader(file, "MAIN", 0, mainBuffer.tellp());
    writeBuffer(file, mainBuffer);

    file.flush();
    file.close();
}
