cmake_minimum_required(VERSION 3.0)
project(bobyqa LANGUAGES C)
add_library(bobyqa bobyqa.c)
target_include_directories(bobyqa PUBLIC .)
