cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
project(ProxyGen)

set(CMAKE_CXX_STANDARD 17)

set(PROXYGEN_WARNINGS
        -Wall
        -Wpedantic
        -Wimplicit-fallthrough)

add_subdirectory(deps/)
add_subdirectory(lib/)

add_executable(Voxelize src/Voxelize.cpp)
target_link_libraries(Voxelize PRIVATE ProxyGenLib)
target_compile_options(Voxelize PRIVATE ${PROXYGEN_WARNINGS})

add_executable(SphereSet src/SphereSet.cpp)
target_link_libraries(SphereSet PRIVATE ProxyGenLib)
target_compile_options(SphereSet PRIVATE ${PROXYGEN_WARNINGS})
