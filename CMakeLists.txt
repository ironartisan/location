cmake_minimum_required(VERSION 3.21)
project(location)

set(CMAKE_CXX_STANDARD 14)

add_library(main1lib Point.h Point.cpp Kmeans.h Kmeans.cpp)
add_executable(main1 test_kmeans.cpp)
target_link_libraries(main1 main1lib)

add_executable(main main.cpp)