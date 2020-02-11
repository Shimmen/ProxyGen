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
