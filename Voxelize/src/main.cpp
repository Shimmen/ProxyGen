#include "SimpleMesh.h"
#include "Texture.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <fmt/format.h>
#include <tiny_gltf.h>
#include <vector>

using Node = tinygltf::Node;
using Model = tinygltf::Model;
using Mesh = tinygltf::Mesh;
using Primitive = tinygltf::Primitive;

mat4 createMatrix(const Node& node)
{
    if (!node.matrix.empty()) {
        return mathkit::linearToMat4(node.matrix);
    } else {
        mat4 translation = node.translation.empty()
            ? mat4(1.0f)
            : mathkit::translate(node.translation[0], node.translation[1], node.translation[2]);
        mat4 rotation = node.rotation.empty()
            ? mat4(1.0f)
            : mathkit::rotate(quat(node.rotation[0], node.rotation[1], node.rotation[2], node.rotation[3]));
        mat4 scale = node.scale.empty()
            ? mat4(1.0f)
            : mathkit::scale(node.scale[0], node.scale[1], node.scale[2]);
        return translation * rotation * scale;
    }
}

std::pair<std::string, Model> loadModel(const std::string& path)
{
    int lastSlash = path.rfind('/');
    if (lastSlash == -1) {
        return {};
    }
    auto basePath = path.substr(0, lastSlash + 1);

    tinygltf::TinyGLTF loader {};

    std::string error;
    std::string warning;

    Model model;
    bool result = loader.LoadASCIIFromFile(&model, &error, &warning, path);

    if (!warning.empty()) {
        fmt::print("glTF loader warning:\n\t{}\n", warning);
    }

    if (!error.empty()) {
        fmt::print("glTF loader error:\n\t{}\n", error);
    }

    if (!result) {
        fmt::print("glTF loader: could not load file: '{}'\n", path);
        exit(1);
    }

    if (model.defaultScene == -1 && model.scenes.size() > 1) {
        fmt::print("glTF loader: scene ambiguity in model '{}'\n", path);
    }

    return { basePath, model };
}

std::vector<vec3> bakedPositionData(const Model& model, const tinygltf::Primitive& primitive, mat4 matrix)
{
    auto entry = primitive.attributes.find("POSITION");
    assert(entry != primitive.attributes.end());
    const tinygltf::Accessor accessor = model.accessors[entry->second];

    assert(accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT);
    assert(accessor.type == TINYGLTF_TYPE_VEC3);

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    assert(view.byteStride == 0); // (i.e. tightly packed)

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];

    const unsigned char* start = buffer.data.data() + view.byteOffset;
    auto* first = reinterpret_cast<const vec3*>(start);

    std::vector<vec3> bakedPositionData;
    bakedPositionData.reserve(accessor.count);

    for (size_t i = 0; i < accessor.count; ++i) {
        vec4 point = vec4(*(first + i), 1.0f);
        bakedPositionData.emplace_back(matrix * point);
    }

    return bakedPositionData;
}

std::vector<vec2> texcoordData(const Model& model, const tinygltf::Primitive& primitive)
{
    auto entry = primitive.attributes.find("TEXCOORD_0");
    assert(entry != primitive.attributes.end());
    const tinygltf::Accessor accessor = model.accessors[entry->second];

    assert(accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT);
    assert(accessor.type == TINYGLTF_TYPE_VEC2);

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    assert(view.byteStride == 0); // (i.e. tightly packed)

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];

    const unsigned char* start = buffer.data.data() + view.byteOffset;
    auto* first = reinterpret_cast<const vec2*>(start);

    std::vector<vec2> vec { first, first + accessor.count };
    return vec;
}

std::vector<size_t> indexData(const Model& model, const tinygltf::Primitive& primitive)
{
    assert(primitive.indices != -1);
    const tinygltf::Accessor accessor = model.accessors[primitive.indices];
    assert(accessor.type == TINYGLTF_TYPE_SCALAR);

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    assert(view.byteStride == 0); // (i.e. tightly packed)

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];
    const unsigned char* start = buffer.data.data() + view.byteOffset;

    std::vector<size_t> indexData;
    indexData.reserve(accessor.count);

    if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
        auto* first = reinterpret_cast<const uint16_t*>(start);
        for (size_t i = 0; i < accessor.count; ++i) {
            size_t index = *(first + i);
            indexData.emplace_back(index);
        }
    } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
        auto* first = reinterpret_cast<const uint32_t*>(start);
        for (size_t i = 0; i < accessor.count; ++i) {
            size_t index = *(first + i);
            indexData.emplace_back(index);
        }
    } else {
        assert(false);
    }

    return indexData;
}

std::string baseColorTextureURI(const Model& model, const Primitive& primitive)
{
    assert(primitive.material != -1);
    auto& material = model.materials[primitive.material];
    int texIndex = material.pbrMetallicRoughness.baseColorTexture.index;
    assert(texIndex != -1);
    auto& texture = model.textures[texIndex];
    assert(texture.source != -1);
    auto& image = model.images[texture.source];
    assert(!image.uri.empty());
    return image.uri;
}

void bakeDownMesh(const Model& model, const std::string& basePath, std::vector<SimpleMesh>& simpleMeshes)
{
    std::function<void(const Node&, mat4)> findMeshesRecursively = [&](const Node& node, mat4 matrix) {
        matrix = createMatrix(node) * matrix;

        if (node.mesh != -1) {
            auto& mesh = model.meshes[node.mesh];
            for (auto& primitive : mesh.primitives) {
                std::string texturePath = basePath + baseColorTextureURI(model, primitive);
                simpleMeshes.emplace_back(bakedPositionData(model, primitive, matrix), texcoordData(model, primitive), indexData(model, primitive), texturePath);
            }
        }

        for (int childIdx : node.children) {
            auto& child = model.nodes[childIdx];
            findMeshesRecursively(child, matrix);
        }
    };

    const tinygltf::Scene& scene = (model.defaultScene != -1)
        ? model.scenes[model.defaultScene]
        : model.scenes.front();

    for (int nodeIdx : scene.nodes) {
        auto& node = model.nodes[nodeIdx];
        findMeshesRecursively(node, createMatrix(node));
    }
}

aabb3 calculateMeshBounds(std::vector<SimpleMesh>& meshes)
{
    vec3 aabbMin;
    vec3 aabbMax;

    for (auto& mesh : meshes) {
        mesh.extendAABB(aabbMin, aabbMax);
    }

    return { aabbMin, aabbMax };
}

int main()
{
    // TODO: Take these as command line parameters!
    std::string path = "../assets/BoomBox/BoomBoxWithAxes.gltf";
    std::string outfile = "../assets/BoomBox.vox";
    size_t gridDimensions = 32;

    auto [basePath, model] = loadModel(path);

    std::vector<SimpleMesh> simpleMeshes {};
    bakeDownMesh(model, basePath, simpleMeshes);

    aabb3 boundsOfAllMeshes = calculateMeshBounds(simpleMeshes);
    VoxelGrid grid { glm::ivec3(gridDimensions), boundsOfAllMeshes };

    fmt::print("= voxelization begin = \n");
    {
        size_t numMeshes = simpleMeshes.size();
        for (size_t i = 0; i < numMeshes; ++i) {
            fmt::print(" mesh {}/{}\n", i + 1, numMeshes);
            grid.insertMesh(simpleMeshes[i]);
        }
    }
    fmt::print("\n= voxelization done  =\n");

    grid.writeToVox(outfile);
}
