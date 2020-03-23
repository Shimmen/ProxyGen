cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
project(ProxyGen)

set(CMAKE_CXX_STANDARD 17)

if(NOT MSVC)
    set(PROXYGEN_WARNINGS
        -Wall
        -Wpedantic
        -Wimplicit-fallthrough)
endif()

add_subdirectory(deps/)
add_subdirectory(lib/)

add_executable(Voxelize src/Voxelize.cpp)
target_link_libraries(Voxelize PRIVATE ProxyGenLib)
target_compile_options(Voxelize PRIVATE ${PROXYGEN_WARNINGS})

add_executable(SphereSet src/SphereSet.cpp)
target_link_libraries(SphereSet PRIVATE ProxyGenLib)
target_link_libraries(SphereSet PRIVATE bobyqa)
target_link_libraries(SphereSet PRIVATE nanoflann::nanoflann)
target_compile_options(SphereSet PRIVATE ${PROXYGEN_WARNINGS})

add_executable(BoundedSphereSet src/BoundedSphereSet.cpp)
target_link_libraries(BoundedSphereSet PRIVATE ProxyGenLib)
target_link_libraries(BoundedSphereSet PRIVATE bobyqa)
target_compile_options(BoundedSphereSet PRIVATE ${PROXYGEN_WARNINGS})
