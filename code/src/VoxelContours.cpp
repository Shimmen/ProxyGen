#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <fmt/format.h>
#include <fstream>
#include <iomanip>
#include <json.hpp>
#include <unordered_map>
#include <vector>

#ifdef WIN32
#define NOMINMAX
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

void setApplicationWorkingDirectory(char* executableName, const std::string& workingDir)
{
#ifdef WIN32
    char fullPathBuf[_MAX_PATH] = {};
    assert(_fullpath(fullPathBuf, executableName, sizeof(fullPathBuf)));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(_chdir(newWorkingDir.c_str()) == 0);
#else
    char fullPathBuf[PATH_MAX] = {};
    assert(realpath(executableName, fullPathBuf));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(chdir(newWorkingDir.c_str()) == 0);
#endif
}

struct Triangle {
    Triangle(vec3 v0, vec3 v1, vec3 v2)
    {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;

        vec3 crossP = cross(v1 - v0, v2 - v0);
        area = length(crossP) / 2.0f;
        normal = normalize(crossP);
    }

    vec3 v[3];
    vec3 normal;
    float area;
};

struct Plane {
    vec3 normal;
    float distance;
};

struct VoxelContour {
    aabb3 aabb;
    Plane plane;
    uint32_t colorIdx;
};

bool triangleNormalsInSameHemisphere(const std::vector<Triangle>& triangles)
{
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& ti = triangles[i];
        for (size_t k = i + 1; k < triangles.size(); ++k) {
            const Triangle& tk = triangles[k];
            float a = dot(ti.normal, tk.normal);
            if (a < 0.0f) {
                return false;
            }
        }
    }

    return true;
}

std::optional<std::pair<float, float>> aabbRayIntersection(aabb3 aabb, vec3 rayOrigin, vec3 rayDirection)
{
    float tmin = (aabb.min.x - rayOrigin.x) / rayDirection.x;
    float tmax = (aabb.max.x - rayOrigin.x) / rayDirection.x;

    if (tmin > tmax)
        std::swap(tmin, tmax);

    float tymin = (aabb.min.y - rayOrigin.y) / rayDirection.y;
    float tymax = (aabb.max.y - rayOrigin.y) / rayDirection.y;

    if (tymin > tymax)
        std::swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
        return {};

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (aabb.min.z - rayOrigin.z) / rayDirection.z;
    float tzmax = (aabb.max.z - rayOrigin.z) / rayDirection.z;

    if (tzmin > tzmax)
        std::swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
        return {};

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

    return { { tmin, tmax } };
}

std::vector<vec3> aabbLineSegmentIntersection(aabb3 aabb, vec3 a, vec3 b)
{
    vec3 rayOrigin = a;
    vec3 rayDirection = normalize(b - a);
    float maxAllowedT = distance(a, b);

    auto intersection = aabbRayIntersection(aabb, rayOrigin, rayDirection);

    if (!intersection.has_value()) {
        return {};
    }

    std::vector<vec3> points {};
    auto& [tmin, tmax] = intersection.value();

    if (tmin >= 0.0f && tmin <= maxAllowedT) {
        points.push_back(rayOrigin + tmin * rayDirection);
    }

    if (tmax >= 0.0f && tmax <= maxAllowedT) {
        points.push_back(rayOrigin + tmax * rayDirection);
    }

    return points;
}

std::vector<vec3> clipTriangleToInternalPoint(aabb3 aabb, const Triangle& triangle)
{
    glm::bvec3 inside = { aabb.contains(triangle.v[0]), aabb.contains(triangle.v[1]), aabb.contains(triangle.v[2]) };

    std::vector<vec3> points {};

    // Add all triangle vertices already inside the aabb
    for (int i = 0; i < 3; ++i) {
        if (inside[i]) {
            points.push_back(triangle.v[i]);
        }
    }

    if (all(inside)) {
        return points;
    }

    for (int i = 0; i < 3; ++i) {

        int iNext = (i + 1) % 3;
        std::vector<vec3> ps = aabbLineSegmentIntersection(aabb, triangle.v[i], triangle.v[iNext]);
        points.insert(points.end(), ps.begin(), ps.end());
    }

    return points;
}

VoxelContour createContourFromConsensusTriangles(aabb3 aabb, const std::vector<Triangle>& triangles)
{
    vec3 sharedNormal = vec3(0.0f);
    for (auto& triangle : triangles) {
        sharedNormal += triangle.normal;
    }
    sharedNormal = normalize(sharedNormal);

    // Find the point furthest back along the shared normal
    float minDistance = std::numeric_limits<float>::max();
    for (auto& triangle : triangles) {
        for (vec3 point : clipTriangleToInternalPoint(aabb, triangle)) {
            float distance = dot(point, sharedNormal);
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
    }

    Plane plane;
    plane.normal = sharedNormal;
    plane.distance = minDistance;

    VoxelContour contour { aabb, plane };
    return contour;
}

VoxelContour createContourForSingleTriangle(aabb3 aabb, const Triangle& triangle)
{
    float distance = dot(triangle.v[0], triangle.normal);
    Plane plane { triangle.normal, distance };

    VoxelContour contour { aabb, plane };
    return contour;
}

std::optional<VoxelContour> createContour(aabb3 aabb, const std::vector<VoxelGrid::TriangleRef>& refs)
{
    std::vector<Triangle> triangles;
    for (auto& ref : refs) {
        vec3 v0, v1, v2;
        ref.mesh.triangle(ref.triangleIndex, v0, v1, v2);
        triangles.emplace_back(v0, v1, v2);
    }

    if (triangles.size() == 1) {
        return createContourForSingleTriangle(aabb, triangles.front());
    }

    if (triangleNormalsInSameHemisphere(triangles)) {
        return createContourFromConsensusTriangles(aabb, triangles);
    }

    /*
    static int numWithTwo = 0;
    static int numWithConsensus = 0;
    numWithTwo += 1;
    numWithConsensus += triangleNormalsInSameHemisphere(triangles) ? 1 : 0;
    fmt::print("       -> {}, consensus {}\n", numWithTwo, numWithConsensus);
    */

    return {};
}

void writeResult(const std::string& outPath, const std::vector<VoxelContour>& contours, const std::vector<vec3>& colors)
{
    using namespace nlohmann;

    json j;
    j["proxy"] = "voxel-contours";

    j["contours"] = {};
    for (const VoxelContour& contour : contours) {
        json jsonContour = {
            { "aabbMin", { contour.aabb.min.x, contour.aabb.min.y, contour.aabb.min.z } },
            { "aabbMax", { contour.aabb.max.x, contour.aabb.max.y, contour.aabb.max.z } },
            { "normal", { contour.plane.normal.x, contour.plane.normal.y, contour.plane.normal.z } },
            { "distance", contour.plane.distance },
            { "colorIndex", contour.colorIdx }
        };
        j["contours"].push_back(jsonContour);
    }

    j["colors"] = {};
    for (const vec3& color : colors) {
        json jsonColor = { color.r, color.g, color.b };
        j["colors"].push_back(jsonColor);
    }

    std::ofstream outstream(outPath);
    outstream << std::setw(4) << j << std::endl;
}

int main(int argc, char** argv)
{
    setApplicationWorkingDirectory(argv[0], "ProxyGen");

#if 1
    std::string path = "assets/Barrel/barrel.gltf";
    std::string outPath = "assets/Barrel/barrel_contours.json";
    constexpr uint32_t gridDimensions = 22;
#else
    std::string path = "assets/Sponza/glTF/Sponza.gltf";
    std::string outPath = "assets/Sponza/sponza_contours.json";
    constexpr uint32_t gridDimensions = 64;
#endif

    constexpr uint32_t numColors = 40;

    fmt::print("=> loading model '{}'\n", path);
    auto [basePath, model] = GltfUtil::loadModel(path);
    std::vector<SimpleMesh> meshes {};
    GltfUtil::bakeDownModelToSimpleMeshes(model, basePath, meshes);

    fmt::print("=> creating voxel grid of resultion {}x{}x{}\n", gridDimensions, gridDimensions, gridDimensions);

    aabb3 meshBounds = SimpleMesh::calculateBounds(meshes);
    VoxelGrid shell { gridDimensions, meshBounds };

    size_t triangleCount = 0;
    for (const SimpleMesh& mesh : meshes) {
        shell.insertMesh(mesh, true);
        triangleCount += mesh.triangleCount();
    }

    fmt::print("==> num triangles:    {}\n", triangleCount);
    fmt::print("==> num shell voxels: {}\n", shell.numFilledVoxels());

    fmt::print("=> quantizing colors (to {} unique)\n", numColors);
    shell.quantizeColors(numColors);

    fmt::print("=> defining contours\n");

    std::vector<VoxelContour> contours {};
    size_t numIgnored = 0;

    shell.forEachFilledVoxel([&](aabb3 aabb, uint32_t value, const std::vector<VoxelGrid::TriangleRef>& triangles) {
        auto maybeContour = createContour(aabb, triangles);
        if (maybeContour.has_value()) {
            VoxelContour contour = maybeContour.value();
            contour.colorIdx = value - 1;
            contours.push_back(contour);
        } else {
            numIgnored += 1;
        }
    });

    fmt::print("=> {} contours generator\n", contours.size());
    fmt::print("=> {} contours ignored\n", numIgnored);

    writeResult(outPath, contours, shell.colors());
    //shell.writeToVox(outPath);
}
