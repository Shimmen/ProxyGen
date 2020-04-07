#pragma once

#include "mathkit.h"
#include <stb_image.h>
#include <string>

class Texture {
public:
    explicit Texture(vec4 color);
    explicit Texture(std::string path);
    ~Texture();

    vec3 sample(vec2 uv) const;

private:
    std::string m_path;

    int m_width { 0 };
    int m_height { 0 };
    int m_numComponents { 0 };
    stbi_uc* m_pixels { nullptr };
};
