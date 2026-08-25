# tests/gates/p2buildtest — rework-cycle local driver (GRAFT G1, F8).
# Name must not match root .gitignore `build*/`.
# Compiles src/refit_build.cpp directly. Does not add the TU to stl2step_core.

cmake_minimum_required(VERSION 3.16)

if(NOT TARGET stl2step::core)
    project(stl2step_p2buildtest LANGUAGES CXX)
    get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
    if(EXISTS "${_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        set(STL2STEP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(STL2STEP_INSTALL OFF CACHE BOOL "" FORCE)
        add_subdirectory("${_root}" "${CMAKE_CURRENT_BINARY_DIR}/_engine")
    else()
        find_package(stl2step REQUIRED)
    endif()
    unset(_root)
endif()

if(NOT TARGET stl2step::harness)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../harness"
                     "${CMAKE_CURRENT_BINARY_DIR}/_harness")
endif()

add_executable(stl2step_p2buildtest
    build_local.cpp
    "${CMAKE_CURRENT_SOURCE_DIR}/../../../src/refit_build.cpp")
target_compile_features(stl2step_p2buildtest PRIVATE cxx_std_17)
target_include_directories(stl2step_p2buildtest PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../../../src")
target_link_libraries(stl2step_p2buildtest PRIVATE
    stl2step::core stl2step::harness)
set_target_properties(stl2step_p2buildtest PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

if(MSVC)
    target_compile_definitions(stl2step_p2buildtest PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(stl2step_p2buildtest PRIVATE /EHsc)
endif()

# Live re-gate (lane p2-real). CTest entry lives under tests/gates/p2real/.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../p2real/CMakeLists.txt")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../p2real"
                     "${CMAKE_CURRENT_BINARY_DIR}/_p2real")
endif()
