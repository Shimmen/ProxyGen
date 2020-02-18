#pragma once

#include "SimpleMesh.h"
#include "mathkit.h"
#include <fmt/format.h>
#include <tiny_gltf.h>

namespace GltfUtil {

std::pair<std::string, tinygltf::Model> loadModel(const std::string& path);
void bakeDownModelToSimpleMeshes(const tinygltf::Model& model, const std::string& basePath, std::vector<SimpleMesh>& simpleMeshes);

mat4 createMatrixForNode(const tinygltf::Node& node);

std::vector<vec3> bakedPositionData(const tinygltf::Model& model, const tinygltf::Primitive& primitive, mat4 matrix);
std::vector<vec2> texcoordData(const tinygltf::Model& model, const tinygltf::Primitive& primitive);
std::vector<size_t> indexData(const tinygltf::Model& model, const tinygltf::Primitive& primitive);

std::string baseColorTextureURI(const tinygltf::Model& model, const tinygltf::Primitive& primitive);

}