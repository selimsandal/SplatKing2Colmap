cmake_minimum_required(VERSION 3.20)

project(sk2colmap LANGUAGES CXX)

add_executable(sk2colmap
    src/main.cpp
    src/sk2cm.cpp
)

target_include_directories(sk2colmap PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)

target_compile_features(sk2colmap PRIVATE cxx_std_20)

if(MSVC)
    target_compile_options(sk2colmap PRIVATE /W4 /permissive-)
else()
    target_compile_options(sk2colmap PRIVATE -Wall -Wextra -Wpedantic)
endif()
