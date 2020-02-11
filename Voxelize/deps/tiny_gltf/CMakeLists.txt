cmake_minimum_required(VERSION 3.0)

project(tiny_gltf LANGUAGES CXX)
add_library(tiny_gltf tiny_gltf.cpp)
target_include_directories(tiny_gltf PRIVATE ../nlohmann_json)
target_include_directories(tiny_gltf PRIVATE ../stb_image)

target_include_directories(tiny_gltf PUBLIC .)

