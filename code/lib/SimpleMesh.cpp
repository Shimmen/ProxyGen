#include "SimpleMesh.h"

#include <algorithm>

SimpleMesh::SimpleMesh(std::vector<vec3>&& positions, std::vector<vec2>&& texcoords, std::vector<size_t>&& indices, Texture texture)
    : m_positions(positions)
    , m_texcoords(texcoords)
    , m_indices(indices)
    , m_texture(texture)
{
}

bool SimpleMesh::hasTexture() const
{
    return true;
}

const Texture& SimpleMesh::texture() const
{
    assert(hasTexture());
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

void SimpleMesh::triangle(size_t triangleIndex, vec3& v0, vec3& v1, vec3& v2) const
{
    uint32_t i0 = m_indices[3 * triangleIndex + 0];
    uint32_t i1 = m_indices[3 * triangleIndex + 1];
    uint32_t i2 = m_indices[3 * triangleIndex + 2];
    v0 = m_positions[i0];
    v1 = m_positions[i1];
    v2 = m_positions[i2];
}

void SimpleMesh::triangleTexcoords(size_t triangleIndex, vec2& uv0, vec2& uv1, vec2& uv2) const
{
    uint32_t i0 = m_indices[3 * triangleIndex + 0];
    uint32_t i1 = m_indices[3 * triangleIndex + 1];
    uint32_t i2 = m_indices[3 * triangleIndex + 2];
    uv0 = m_texcoords[i0];
    uv1 = m_texcoords[i1];
    uv2 = m_texcoords[i2];
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

aabb3 SimpleMesh::calculateBounds(std::vector<SimpleMesh>& meshes)
{
    if (meshes.empty()) {
        return { vec3(-INFINITY), vec3(+INFINITY) };
    }

    vec3 aabbMin { +INFINITY };
    vec3 aabbMax { -INFINITY };

    for (auto& mesh : meshes) {
        mesh.extendAABB(aabbMin, aabbMax);
    }

    return { aabbMin, aabbMax };
}
