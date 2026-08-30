# tests/tools/prismscan — corpus routing census (AC3-P1).
cmake_minimum_required(VERSION 3.16)

get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)

if(NOT TARGET stl2step::core)
    project(stl2step_prismscan LANGUAGES CXX)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL
            "prismscan standalone: do not recurse into tests/" FORCE)
        set(STL2STEP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(STL2STEP_INSTALL OFF CACHE BOOL "" FORCE)
        add_subdirectory("${_stl2step_root}" "${CMAKE_CURRENT_BINARY_DIR}/_engine")
    else()
        find_package(stl2step REQUIRED)
    endif()
endif()

if(NOT TARGET stl2step::harness)
    add_subdirectory("${_stl2step_root}/tests/harness"
                     "${CMAKE_CURRENT_BINARY_DIR}/_harness")
endif()

add_executable(prismscan prismscan.cpp)
target_compile_features(prismscan PRIVATE cxx_std_17)
set_target_properties(prismscan PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(prismscan PRIVATE stl2step::harness)
target_include_directories(prismscan PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/tests/harness")
if(MSVC)
    target_compile_definitions(prismscan PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(prismscan)
endif()
unset(_stl2step_root)
