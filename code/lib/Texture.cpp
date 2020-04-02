#include "Texture.h"

#include <fmt/format.h>

Texture::Texture(std::string path)
    : m_path(std::move(path))
{
    m_pixels = stbi_load(m_path.c_str(), &m_width, &m_height, &m_numComponents, STBI_default);
    if (!m_pixels) {
        fmt::print("Texture: could not load image '{}'\n", m_path);
        m_pixels = nullptr;
    }
}

Texture::~Texture()
{
    // todo: maybe don't leak the stuff?
    //if (m_pixels) {
    //    stbi_image_free(m_pixels);
    //}
}

vec3 Texture::sample(vec2 uv) const
{
    uv = vec2(std::max(0.0f, std::min(uv.x, 1.0f)),
              std::max(0.0f, std::min(uv.y, 1.0f)));

    vec2 pixelCoords = uv * vec2(m_width - 1, m_height - 1);
    ivec2 nearest = round(pixelCoords);

    int idx = m_numComponents * (nearest.x + nearest.y * m_width);

    float r = float(m_pixels[idx + 0]) / UINT8_MAX;
    if (m_numComponents == 1) {
        return { r, r, r };
    }

    float g = float(m_pixels[idx + 1]) / UINT8_MAX;
    if (m_numComponents == 2) {
        return { r, g, 0.0f };
    }

    if (m_numComponents >= 3) {
        float b = float(m_pixels[idx + 2]) / UINT8_MAX;
        return { r, g, b };
    }

    assert(false);
    return {};
}
