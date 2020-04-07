#include "VoxelGrid.h"

#include "tribox3.h"
#include <algorithm>
#include <fmt/format.h>
#include <fstream>
#include <sstream>

#define MIN(a, b) ((a < b) ? a : b)

VoxelGrid::VoxelGrid(glm::ivec3 size, aabb3 bounds)
    : m_bounds(bounds)
    , m_sx(size.x)
    , m_sy(size.y)
    , m_sz(size.z)
{
    size_t totalSize = size.x * size.y * size.z;
    assert(totalSize < 4ull * 1024ull * 1024ull * 1024ull);

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

vec3 VoxelGrid::voxelCenterPoint(ivec3 gridCoord) const
{
    const vec3 voxSize = voxelSize();
    return m_bounds.min + vec3(gridCoord) * voxSize + (0.5f * voxSize);
}

ivec3 VoxelGrid::gridDimensions() const
{
    return ivec3(m_sx, m_sy, m_sz);
}

aabb3 VoxelGrid::gridBounds() const
{
    return m_bounds;
}

vec3 VoxelGrid::voxelSize() const
{
    return (m_bounds.max - m_bounds.min) / vec3(m_sx, m_sy, m_sz);
}

void VoxelGrid::forEachFilledVoxel(std::function<void(aabb3 aabb, const std::vector<TriangleRef>&)> callback) const
{
    vec3 voxelHalfSize = 0.5f * voxelSize();

    for (size_t z = 0; z < m_sz; ++z) {
        for (size_t y = 0; y < m_sy; ++y) {
            for (size_t x = 0; x < m_sx; ++x) {

                if (get(x, y, z) > 0) {

                    uint64_t linearIdx = linearIndex(x, y, z);

                    auto entry = m_trianglesForVoxelIndex.find(linearIdx);
                    if (entry != m_trianglesForVoxelIndex.end()) {

                        const std::vector<TriangleRef>& triRefs = entry->second;

                        vec3 gridCenter = voxelCenterPoint({ x, y, z });
                        vec3 min = gridCenter - voxelHalfSize;
                        vec3 max = gridCenter + voxelHalfSize;

                        callback({ min, max }, triRefs);
                    }
                }
            }
        }
    }
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

std::vector<vec3> VoxelGrid::filledVoxelsInSphere(const Sphere& sphere) const
{
    std::vector<vec3> voxels {};

    vec3 voxSize = voxelSize();
    vec3 first = sphere.origin - vec3(sphere.radius) - (0.05f * voxSize);
    vec3 last = sphere.origin + vec3(sphere.radius) + (0.05f * voxSize);

    first = max(first, m_bounds.min);
    last = min(last, m_bounds.max);

    for (float z = first.z; z <= last.z; z += voxSize.z) { // NOLINT(cert-flp30-c)
        for (float y = first.y; y <= last.y; y += voxSize.y) { // NOLINT(cert-flp30-c)
            for (float x = first.x; x <= last.x; x += voxSize.x) { // NOLINT(cert-flp30-c)

                vec3 samplePoint = vec3(x, y, z);

                // Only count voxels that actually are inside the sphere
                if (distance2(sphere.origin, samplePoint) > sphere.radius * sphere.radius) {
                    continue;
                }

                ivec3 gridSampleCoord = remapToGridSpace(samplePoint, std::round);

                if (get(gridSampleCoord) > 0) {
                    voxels.emplace_back(samplePoint);
                }
            }
        }
    }

    return voxels;
}

uint64_t VoxelGrid::linearIndex(ivec3 gridIndex) const
{
    return linearIndex(gridIndex.x, gridIndex.y, gridIndex.z);
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

uint32_t VoxelGrid::get(ivec3 coord) const
{
    return get(coord.x, coord.y, coord.z);
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

std::vector<ivec3> VoxelGrid::immediateFilledNeighbors(ivec3 gridCoord) const
{
    std::vector<ivec3> neighbors {};
    neighbors.reserve(3 * 3 * 3 - 1);

    for (int zz = -1; zz <= 1; ++zz) {
        for (int yy = -1; yy <= 1; ++yy) {
            for (int xx = -1; xx <= 1; ++xx) {

                if (xx == 0 && yy == 0 && zz == 0) {
                    continue;
                }

                int x = gridCoord.x + xx;
                int y = gridCoord.y + yy;
                int z = gridCoord.z + zz;

                if (x < 0 || x >= m_sx || y < 0 || y >= m_sy || z < 0 || z >= m_sz) {
                    continue;
                }

                if (get(x, y, z) > 0) {
                    neighbors.emplace_back(x, y, z);
                }
            }
        }
    }

    return neighbors;
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

                            vec3 normal = cross(e1, e2);
                            vec3 w = voxelCenter - v[0];

                            float gamma = dot(normal, cross(e1, w)) / length2(normal);
                            float beta = dot(normal, cross(w, e2)) / length2(normal);
                            float alpha = 1.0f - gamma - beta;

                            // NOTE: In theory we might want to check that the barycentric constraints apply, but since we just take a single
                            //  sample from the voxelCenter it's not going to be accurate anyway, and the UVs are probably going to be fine..
                            //constexpr float eps = 0.5f;
                            //constexpr float lo = 0.0f - eps;
                            //constexpr float hi = 1.0f + eps;
                            //if (alpha >= lo && alpha <= hi && beta >= lo && beta <= hi && gamma >= lo && gamma <= hi)

                            vec2 uv0, uv1, uv2;
                            mesh.triangleTexcoords(ti, uv0, uv1, uv2);

                            vec2 uv = alpha * uv0 + beta * uv1 + gamma * uv2;

                            vec3 color = vec3(1, 0, 1);
                            if (mesh.hasTexture()) {
                                color = mesh.texture().sample(uv);
                            }

                            value = m_colors.size() + 1;
                            m_colors.push_back(color);
                        }

                        set(x, y, z, value);

                        uint64_t index = linearIndex(x, y, z);
                        m_trianglesForVoxelIndex[index].push_back({ mesh, ti });
                    }
                }
            }
        }
    }
}

void VoxelGrid::fillVolumes()
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
                    mesh.triangle(triangleIdx, v0, v1, v2);

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

void VoxelGrid::subtractGrid(const VoxelGrid& other)
{
    assert(gridDimensions() == other.gridDimensions());
    assert(gridBounds().min == other.gridBounds().min);
    assert(gridBounds().max == other.gridBounds().max);

    for (size_t z = 0; z < m_sz; ++z) {
        for (size_t y = 0; y < m_sy; ++y) {
            for (size_t x = 0; x < m_sx; ++x) {

                uint8_t sub = other.get(x, y, z);
                if (sub > 0) {
                    uint8_t current = get(x, y, z);
                    set(x, y, z, (sub >= current) ? current - sub : 0);
                }
            }
        }
    }
}

void VoxelGrid::quantizeColors(uint32_t numBins)
{
    assert(numBins > 0);
    constexpr size_t numIterations = 4; // TODO!

    //
    // Perform k-means clustering
    //

    struct QColor {
        vec3 color;
        int clusterIndex;
    };

    std::vector<vec3> clusters {};
    for (uint32_t k = 0; k < numBins; ++k) {
        size_t pointIdx = rand() % m_colors.size(); // NOLINT(cert-msc30-c,cert-msc50-cpp)
        clusters.push_back(m_colors[pointIdx]);
    }

    std::vector<QColor> quantizedColors {};
    for (const vec3& color : m_colors) {
        quantizedColors.push_back({ color, -1 });
    }

    for (size_t it = 0; it < numIterations; ++it) {

        // Find the closest cluster for each point

        for (QColor& qColor : quantizedColors) {

            int closestIndex = -1;
            float closestVal = 999.99f;

            for (uint32_t k = 0; k < numBins; ++k) {
                const vec3& clusterCenter = clusters[k];
                float dist2 = distance2(qColor.color, clusterCenter);
                if (dist2 < closestVal) {
                    closestVal = dist2;
                    closestIndex = k;
                }
            }

            qColor.clusterIndex = closestIndex;
        }

        // Move cluster centers to average position of all members

        struct AvgPair {
            vec3 sum { 0.0f };
            int count { 0 };
        };
        std::vector<AvgPair> averagePairs { numBins };

        for (QColor& qColor : quantizedColors) {
            AvgPair& avgPair = averagePairs[qColor.clusterIndex];
            avgPair.sum += qColor.color;
            avgPair.count += 1;
        }

        for (uint32_t k = 0; k < numBins; ++k) {
            AvgPair& avgPair = averagePairs[k];
            clusters[k] = avgPair.sum / vec3(avgPair.count);
        }
    }

    //
    // Reassign color indices in the grid
    //

    size_t totalSize = m_sx * m_sy * m_sz;
    for (size_t i = 0; i < totalSize; ++i) {
        uint32_t value = m_grid[i];
        if (value > 0) {
            uint32_t colorIndex = value - 1;
            QColor& quantized = quantizedColors[colorIndex];
            m_grid[i] = quantized.clusterIndex + 1;
        }
    }

    m_colors = clusters;
}

void VoxelGrid::writeToVox(const std::string& path) const
{
    if (m_sx > 126 || m_sy > 126 || m_sz > 126) {
        fmt::print("VoxelGrid::writeToVox(): the .vox format only supports resultions up to 126x126x126, so parts of the voxel grid might be cut off\n");
    }

    auto writeLabel = [](std::ostream& stream, const char* label) {
        stream.write(label, strlen(label));
    };
    auto writeInt32 = [](std::ostream& stream, int32_t value) {
        stream.write((const char*)&value, 4);
    };
    auto writeUInt8 = [](std::ostream& stream, uint8_t value) {
        stream.write((const char*)&value, 1);
    };
    auto writeChunkHeader = [&](std::ostream& stream, const char* chunkId, int32_t sizeOfSelf, int32_t sizeOfChildren) {
        writeLabel(stream, chunkId);
        writeInt32(stream, sizeOfSelf);
        writeInt32(stream, sizeOfChildren);
    };
    auto writeBuffer = [&](std::ostream& stream, std::ostringstream& source) {
        std::string content = source.str();
        stream << content;
    };

    int32_t sx = MIN(m_sx, 126ull);
    int32_t sy = MIN(m_sy, 126ull);
    int32_t sz = MIN(m_sz, 126ull);

    std::ostringstream sizeBuffer {};
    writeInt32(sizeBuffer, sx);
    writeInt32(sizeBuffer, sy);
    writeInt32(sizeBuffer, sz);

    std::ostringstream xyziBuffer {};
    writeInt32(xyziBuffer, numFilledVoxels());

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

    if (!m_colors.empty()) {
        if (m_colors.size() > 256) {
            fmt::print("VoxelGrid::writeToVox(): more than 256 colors defined, will be truncated!\n");
        }
        std::ostringstream rgbaBuffer {};

        for (int i = 0; i < 256; ++i) {
            if (i < m_colors.size()) {
                ivec3 color = ivec3(m_colors[i] * 255.99f);
                int32_t colorInt = (0xFF << 24) | (color.b << 16) | (color.g << 8) | (color.r); // NOLINT(hicpp-signed-bitwise)
                writeInt32(rgbaBuffer, colorInt);
            } else {
                writeInt32(rgbaBuffer, 0xFFFF00FF);
            }
        }

        writeChunkHeader(mainBuffer, "RGBA", rgbaBuffer.tellp(), 0);
        writeBuffer(mainBuffer, rgbaBuffer);
    }

    std::ofstream file;
    file.open(path, std::ios::binary | std::ios::out);

    writeLabel(file, "VOX ");
    writeInt32(file, 150);

    writeChunkHeader(file, "MAIN", 0, mainBuffer.tellp());
    writeBuffer(file, mainBuffer);

    file.flush();
    file.close();
}
