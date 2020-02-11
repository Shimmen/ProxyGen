#include "SimpleMesh.h"

SimpleMesh::SimpleMesh(std::vector<vec3>&& positions, std::vector<vec2>&& texcoords, std::vector<size_t>&& indices, const std::string& texturePath)
    : m_positions(positions)
    , m_texcoords(texcoords)
    , m_indices(indices)
    , m_texture(texturePath)
{
}

const Texture& SimpleMesh::texture() const
{
    return m_texture;
}

size_t SimpleMesh::vertexCount() const
{
    assert(m_positions.size() == m_texcoords.size());
    return m_positions.size();
}

const std::vector<vec3>& SimpleMesh::positions() const
{
    return m_positions;
}

const std::vector<vec2>& SimpleMesh::texcoords() const
{
    return m_texcoords;
}

size_t SimpleMesh::triangleCount() const
{
    assert(m_indices.size() % 3 == 0);
    return m_indices.size() / 3;
}

void SimpleMesh::triangle(size_t triangleIndex, float triangleVerts[3][3]) const
{
    uint32_t i0 = m_indices[3 * triangleIndex + 0];
    triangleVerts[0][0] = m_positions[i0].x;
    triangleVerts[0][1] = m_positions[i0].y;
    triangleVerts[0][2] = m_positions[i0].z;

    uint32_t i1 = m_indices[3 * triangleIndex + 1];
    triangleVerts[1][0] = m_positions[i1].x;
    triangleVerts[1][1] = m_positions[i1].y;
    triangleVerts[1][2] = m_positions[i1].z;

    uint32_t i2 = m_indices[3 * triangleIndex + 2];
    triangleVerts[2][0] = m_positions[i2].x;
    triangleVerts[2][1] = m_positions[i2].y;
    triangleVerts[2][2] = m_positions[i2].z;
}

void SimpleMesh::triangle(size_t triangleIndex, vec3& v0, vec3& v1, vec3& v2) const
{
    uint32_t i0 = m_indices[3 * triangleIndex + 0];
    uint32_t i1 = m_indices[3 * triangleIndex + 0];
    uint32_t i2 = m_indices[3 * triangleIndex + 0];
    v0 = m_positions[i0];
    v1 = m_positions[i1];
    v2 = m_positions[i2];
}

void SimpleMesh::extendAABB(vec3& min, vec3& max) const
{
    for (const vec3& point : m_positions) {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);

        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);
    }
}
