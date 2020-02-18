
add_library(ProxyGenLib
        Texture.cpp
        SimpleMesh.cpp
        VoxelGrid.cpp
        GltfUtil.cpp)

target_compile_options(ProxyGenLib PRIVATE ${PROXYGEN_WARNINGS})

target_include_directories(ProxyGenLib PUBLIC .)
target_include_directories(ProxyGenLib PUBLIC ../deps/glm-0.9.9.6)

target_link_libraries(ProxyGenLib PUBLIC fmt)
target_link_libraries(ProxyGenLib PUBLIC stb_image)
target_link_libraries(ProxyGenLib PUBLIC tiny_gltf)
