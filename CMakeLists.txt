cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 11)

project(MPM2DProject)

add_library(MPM2DLib "")

add_subdirectory(ace_mate)
add_subdirectory(include)
add_subdirectory(src)

target_include_directories(
  MPM2DLib
  PUBLIC
    "include"
    "ace_mate"
)


add_executable(MPM2D  main.cpp)

target_link_libraries(MPM2D MPM2DLib)
