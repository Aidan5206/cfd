cmake_minimum_required(VERSION 3.30)
project(cfd)

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_TOOLCHAIN_FILE "/Users/aidanhsu/.vcpkg-clion/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")

add_executable(cfd main.cpp
        main.h
        Grid.cpp
        Grid.h
        functions.cpp
        functions.h)

find_package(Matplot++ CONFIG REQUIRED)
target_link_libraries(cfd PRIVATE Matplot++::cimg  Matplot++::matplot)
