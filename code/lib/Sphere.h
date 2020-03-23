#pragma once

#include "mathkit.h"

struct Sphere {
public:
    vec3 origin;
    double radius;

    [[nodiscard]] double volume() const
    {
        return (4.0 / 3.0) * mathkit::PI * (radius * radius * radius);
    }

    double overlapWith(const Sphere& other) const
    {
        // See https: //www.researchgate.net/publication/289097700_Volume_of_intersection_of_six_spheres_A_special_case_of_practical_interest
        // Maybe also https://math.stackexchange.com/questions/297751/overlapping-spheres.

        double d = distance(origin, other.origin);

        const double& r1 = radius;
        const double& r2 = other.radius;

        if (d >= r1 + r2) {
            return 0.0;
        }

        if (d <= r2 - r1) {
            return this->volume();
        }

        if (d <= r1 - r2) {
            return other.volume();
        }

        using namespace mathkit;
        return (PI * square(r1 + r2 - d) * (square(r1 + r2 + d) - 4.0 * (r1 * r1 - r1 * r2 + r2 * r2))) / (12.0 * d);
    }
};
