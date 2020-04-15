cmake_minimum_required(VERSION 3.0)
project(spherical_harmonics LANGUAGES CXX)
add_library(spherical_harmonics spherical_harmonics.cc spherical_harmonics.h image.h default_image.cc default_image.h)
target_include_directories(spherical_harmonics PUBLIC .)

target_link_libraries(spherical_harmonics PUBLIC eigen)
target_include_directories(spherical_harmonics PUBLIC ../eigen-3.3.7/)
