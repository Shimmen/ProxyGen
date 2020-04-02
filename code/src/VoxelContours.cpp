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

struct VoxelContour {
    aabb3 aabb;
    vec3 normal;
    float centerOffset;
    //uint32_t colorIdx;
};

VoxelContour createContourForSingleTriangle(aabb3 aabb, const Triangle& triangle)
{
    // Project voxel center on the triangle (or at least the distance)
    vec3 center = 0.5f * (aabb.min + aabb.max);
    vec3 relPoint = center - triangle.v[0];
    float triToCenterOffset = dot(relPoint, triangle.normal);

    VoxelContour contour {};
    contour.aabb = aabb;
    contour.normal = triangle.normal;
    contour.centerOffset = -triToCenterOffset; // -1e-6f; (maybe give some margin too?)

    return contour;
}

void writeResult(const std::string& outPath, const std::vector<VoxelContour>& contours)
{
    using namespace nlohmann;

    json j;
    j["proxy"] = "voxel-contours";
    j["contours"] = {};

    for (const VoxelContour& contour : contours) {
        json jsonContour = {
            { "aabbMin", { contour.aabb.min.x, contour.aabb.min.y, contour.aabb.min.z } },
            { "aabbMax", { contour.aabb.max.x, contour.aabb.max.y, contour.aabb.max.z } },
            { "normal", { contour.normal.x, contour.normal.y, contour.normal.z } },
            { "centerOffset", contour.centerOffset }
        };
        j["contours"].push_back(jsonContour);
    }

    std::ofstream outstream(outPath);
    outstream << std::setw(4) << j << std::endl;
}

int main(int argc, char** argv)
{
    setApplicationWorkingDirectory(argv[0], "ProxyGen");

#if 1
    std::string path = "assets/Bunny/bunny.gltf";
    std::string outPath = "assets/Bunny/bunny_contours.json";
    constexpr uint32_t gridDimensions = 64;
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
    VoxelGrid shell { glm::ivec3(gridDimensions), meshBounds };

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

    shell.forEachFilledVoxel([&](aabb3 aabb, const std::vector<VoxelGrid::TriangleRef>& triangles) {
        if (triangles.size() == 1) {

            auto ref = triangles.front();

            vec3 v0, v1, v2;
            ref.mesh.triangle(ref.triangleIndex, v0, v1, v2);

            Triangle triangle { v0, v1, v2 };
            VoxelContour contour = createContourForSingleTriangle(aabb, triangle);
            contours.push_back(contour);

        } else {
            // TODO: Implement this case!
            numIgnored += 1;

            VoxelContour contour {};
            contour.aabb = aabb;
            contour.normal = vec3(0, 1, 0);
            contour.centerOffset = 0.0f;
            contours.push_back(contour);
        }
    });

    fmt::print("=> {} contours generator\n", contours.size());
    fmt::print("=> {} contours ignored\n", numIgnored);

    writeResult(outPath, contours);
    //shell.writeToVox(outPath);
}
