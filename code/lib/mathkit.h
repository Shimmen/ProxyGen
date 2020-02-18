#pragma once

#include <cmath>
#include <vector>

#define GLM_FORCE_CXX17

#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

// (exposed for GLSL shader shared headers & convenience)
using ivec2 = glm::ivec2;
using ivec3 = glm::ivec3;
using vec2 = glm::vec2;
using vec3 = glm::vec3;
using vec4 = glm::vec4;
using mat3 = glm::mat3;
using mat4 = glm::mat4;
using quat = glm::quat;

struct aabb3 {
    vec3 min;
    vec3 max;
};

namespace mathkit {

#include <glm/gtc/matrix_transform.hpp>

constexpr float PI = 3.14159265359f;
constexpr float TWO_PI = 2.0f * PI;

constexpr vec3 globalX = vec3 { 1.0f, 0.0f, 0.0f };
constexpr vec3 globalY = vec3 { 0.0f, 1.0f, 0.0f };
constexpr vec3 globalZ = vec3 { 0.0f, 0.0f, 1.0f };

constexpr vec3 globalRight = globalX;
constexpr vec3 globalUp = globalY;
constexpr vec3 globalForward = -globalZ;

constexpr float radians(float degrees)
{
    return degrees / 180.0f * PI;
}

constexpr float square(float x)
{
    return x * x;
}

inline float clamp(float val, float minVal, float maxVal)
{
    return glm::clamp(val, minVal, maxVal);
}

inline float mix(float a, float b, float blend)
{
    return glm::mix(a, b, blend);
}

inline float length2(const vec3& v)
{
    return glm::dot(v, v);
}

inline quat axisAngle(const vec3& axis, float angle)
{
    return glm::angleAxis(angle, axis);
}

inline mat4 axisAngleMatrix(const vec3& axis, float angle)
{
    quat quaternion = axisAngle(axis, angle);
    mat4 matrix = glm::mat4_cast(quaternion);
    return matrix;
}

inline mat4 lookAt(const vec3& eye, const vec3& target, const vec3& up = globalY)
{
    return glm::lookAt(eye, target, up);
}

inline quat quatLookAt(const vec3& direction, const vec3& up = globalY)
{
    return glm::quatLookAt(direction, up);
}

inline mat4 translate(float x, float y, float z)
{
    return glm::translate(mat4 { 1.0f }, { x, y, z });
}

inline mat4 scale(float x, float y, float z)
{
    return glm::scale(mat4 { 1.0f }, { x, y, z });
}

inline mat4 scale(float factor)
{
    return glm::scale(mat4 { 1.0f }, { factor, factor, factor });
}

inline mat4 rotate(quat quaternion)
{
    return glm::toMat4(quaternion);
}

inline mat4 infinitePerspective(float fieldOfViewY, float aspectRatio, float zNear)
{
    mat4 matrix = glm::infinitePerspective(fieldOfViewY, aspectRatio, zNear);
    matrix[1][1] *= -1.0f; // (flip for OpenGL -> Vulkan conventions)
    return matrix;
}

inline mat4 orthographicProjection(float size, float aspectRatio, float zNear, float zFar)
{
    float half = size / 2.0f;
    mat4 matrix = glm::orthoRH_ZO(-half * aspectRatio, half * aspectRatio, half, -half, zNear, zFar);
    return matrix;
}

inline vec3 rotateWithQuaternion(vec3 vector, quat rotation)
{
    return glm::rotate(rotation, vector);
}

inline mat4 linearToMat4(const std::vector<double>& vals)
{
    if (vals.empty()) {
        return mat4(1.0f);
    }

    assert(vals.size() == 16);
    return glm::make_mat4(vals.data());
}

}