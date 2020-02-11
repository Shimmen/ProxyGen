cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
project(Voxelize)

set(CMAKE_CXX_STANDARD 17)

add_executable(Voxelize
        src/main.cpp
        src/Texture.cpp
        src/SimpleMesh.cpp
        src/VoxelGrid.cpp)

target_compile_options(Voxelize PRIVATE -Wall)
target_compile_options(Voxelize PRIVATE -Wimplicit-fallthrough)

target_include_directories(Voxelize PRIVATE src/)
target_include_directories(Voxelize PRIVATE deps/glm-0.9.9.6)
target_include_directories(Voxelize PRIVATE deps/nlohmann_json)

add_subdirectory(deps/fmt)
target_link_libraries(Voxelize fmt)

add_subdirectory(deps/stb_image)
target_link_libraries(Voxelize stb_image)

add_subdirectory(deps/tiny_gltf)
target_link_libraries(Voxelize tiny_gltf)
