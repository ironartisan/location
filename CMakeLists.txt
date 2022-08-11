cmake_minimum_required(VERSION 3.21)
project(location)

set(CMAKE_CXX_STANDARD 14)

add_library(main1lib src/algorithims/Point.h src/algorithims/Point.cpp src/algorithims/Kmeans.h src/algorithims/Kmeans.cpp)
add_executable(main1 src/algorithims/test_kmeans.cpp)
target_link_libraries(main1 main1lib)

add_executable(main src/algorithims/main.cpp)