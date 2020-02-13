#pragma once

#include "Texture.h"
#include "mathkit.h"
#include <vector>

class SimpleMesh {
public:
    SimpleMesh(std::vector<vec3>&& positions, std::vector<vec2>&& texcoords, std::vector<size_t>&& indices, const std::string& texturePath);

    [[nodiscard]] const Texture& texture() const;
    [[nodiscard]] size_t vertexCount() const;

    [[nodiscard]] const std::vector<vec3>& positions() const;
    [[nodiscard]] const std::vector<vec2>& texcoords() const;

    [[nodiscard]] size_t triangleCount() const;
    void triangle(size_t triangleIndex, vec3& v0, vec3& v1, vec3& v2) const;
    void triangleTexcoords(size_t triangleIndex, vec2& uv0, vec2& uv1, vec2& uv2) const;

    void extendAABB(vec3& min, vec3& max) const;

    static aabb3 calculateBounds(std::vector<SimpleMesh>& meshes);

private:
    std::vector<vec3> m_positions;
    std::vector<vec2> m_texcoords;
    std::vector<size_t> m_indices;
    Texture m_texture;
};
