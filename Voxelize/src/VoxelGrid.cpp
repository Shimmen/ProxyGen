#include "VoxelGrid.h"

#include <fmt/printf.h>
#include <fstream>
#include <sstream>

VoxelGrid::VoxelGrid(glm::ivec3 size, aabb3 bounds)
    : m_bounds(bounds)
    , m_sx(size.x)
    , m_sy(size.y)
    , m_sz(size.z)
{
    size_t totalSize = size.x * size.y * size.z;

    constexpr size_t maxReasonableSize { 4ul * 1024ul * 1024ul * 1024ul };
    assert(totalSize < maxReasonableSize);

    // note: this will allocate the whole dense grid
    m_grid.resize(totalSize);
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

uint8_t VoxelGrid::get(int x, int y, int z) const
{
    assert(x >= 0 && x < m_sx && y >= 0 && y < m_sy && z >= 0 && z < m_sz);
    return m_grid[x + y * (m_sx) + z * (m_sx * m_sy)];
}

void VoxelGrid::set(int x, int y, int z, uint8_t value)
{
    assert(x >= 0 && x < m_sx && y >= 0 && y < m_sy && z >= 0 && z < m_sz);
    m_grid[x + y * (m_sx) + z * (m_sx * m_sy)] = value;
}

void VoxelGrid::set(vec3 point, uint8_t value)
{
    vec3 normalized = (point - m_bounds.min) / (m_bounds.max - m_bounds.min);
    vec3 gridSpace = normalized * vec3(m_sx - 1, m_sy - 1, m_sz - 1);

    int x = (int)round(gridSpace.x);
    int y = (int)round(gridSpace.y);
    int z = (int)round(gridSpace.z);
    set(x, y, z, value);
}

void VoxelGrid::insertMesh(const SimpleMesh& mesh)
{
    // TODO: Perform proper triangle-aabb intersection!

    // TODO: Use colors instead!
    uint8_t value = 1;

    for (size_t i = 0; i < mesh.vertexCount(); ++i) {
        vec3 position = mesh.positions()[i];
        set(position, value);
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

    std::ostringstream sizeBuffer {};
    writeUInt32(sizeBuffer, m_sx);
    writeUInt32(sizeBuffer, m_sy);
    writeUInt32(sizeBuffer, m_sz);

    std::ostringstream xyziBuffer {};
    writeUInt32(xyziBuffer, numFilledVoxels());

    for (size_t z = 0; z < m_sz; ++z) {
        int numWritten = 0;
        for (size_t y = 0; y < m_sy; ++y) {
            for (size_t x = 0; x < m_sx; ++x) {
                uint8_t voxel = get(x, y, z);
                if (voxel > 0) {
                    if (x > 126 || y > 126 || z > 126) {
                        continue;
                    }
                    writeUInt8(xyziBuffer, m_sx - x - 1);
                    writeUInt8(xyziBuffer, z);
                    writeUInt8(xyziBuffer, y);
                    writeUInt8(xyziBuffer, voxel);
                    numWritten += 1;
                }
            }
        }
        fmt::print("Wrote {} voxels for z={}\n", numWritten, z);
        numWritten = 0;
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
