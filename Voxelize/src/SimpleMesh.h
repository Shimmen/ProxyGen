#pragma once

#include "Texture.h"
#include "mathkit.h"
#include <vector>

class SimpleMesh {
public:
    SimpleMesh(std::vector<vec3>&& positions, std::vector<vec2>&& texcoords, const std::string& texturePath);

    [[nodiscard]] const Texture& texture() const;
    [[nodiscard]] size_t vertexCount() const;

    [[nodiscard]] const std::vector<vec3>& positions() const;
    [[nodiscard]] const std::vector<vec2>& texcoords() const;

    void extendAABB(vec3& min, vec3& max) const;

private:
    std::vector<vec3> m_positions;
    std::vector<vec2> m_texcoords;
    Texture m_texture;
};
