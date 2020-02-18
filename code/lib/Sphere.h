#pragma once

#include "mathkit.h"

struct Sphere {
public:
    vec3 origin;
    float radius;

    [[nodiscard]] float volume() const
    {
        return (4.0f / 3.0f) * mathkit::PI * (radius * radius * radius);
    }
};
