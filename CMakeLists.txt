cmake_minimum_required(VERSION 3.12)

project(ising++ VERSION 0.1 LANGUAGES CXX)

# C++ language settings
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Linting settings
set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)
# Code editor settings
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build structure settings
# add_subdirectory(src build)
add_executable(ising++ ising.cpp main.cpp)

# First, we set global error flags that should work with everything
target_compile_options(ising++ PRIVATE -Wall -Werror)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(ising++ PRIVATE -g0 -O3)
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_options(ising++ PRIVATE -g  -O3)
elseif (CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    target_compile_options(ising++ PRIVATE -g0 -Os)
    target_link_options(ising++ PRIVATE -s)
else() # Assume debug
    target_compile_options(ising++ PRIVATE -g  -O0 -ftime-trace)
endif()

# External library settings
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(boost_program_options 1.70 REQUIRED COMPONENTS program_options)

target_link_libraries(ising++ Boost::program_options png)

