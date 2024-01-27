cmake_minimum_required(VERSION 3.27)
project(BlackHoleRayTracer)
include_directories(imgui-master)
set(CMAKE_CXX_STANDARD 17)
set(IMGUI_DIR "C:/Users/alvel/Desktop/Dev/C++/BlackHoleRayTracer/imgui-master")

find_package(OpenCL REQUIRED)
set(SFML_DIR "c:/Users/alvel/SFML-2.6.1/lib/cmake/SFML")

find_package(SFML COMPONENTS system graphics window REQUIRED)
add_subdirectory(imgui-sfml-2.6.x)

add_executable(BlackHoleRayTracer main.cpp)

target_link_libraries(BlackHoleRayTracer sfml-system sfml-graphics sfml-window OpenCL::OpenCL ImGui-SFML::ImGui-SFML )