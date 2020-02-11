#include "SimpleMesh.h"

SimpleMesh::SimpleMesh(std::vector<vec3>&& positions, std::vector<vec2>&& texcoords, const std::string& texturePath)
    : m_positions(positions)
    , m_texcoords(texcoords)
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
