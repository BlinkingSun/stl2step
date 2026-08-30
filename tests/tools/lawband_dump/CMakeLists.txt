# tests/tools/lawband_dump — throwaway dump of lawChainAccept on labeled bands.
cmake_minimum_required(VERSION 3.16)

get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)

if(NOT TARGET stl2step::core)
    project(stl2step_lawband_dump LANGUAGES CXX)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL
            "lawband_dump standalone: do not recurse into tests/" FORCE)
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

add_executable(stl2step_lawband_dump lawband_dump.cpp)
target_compile_features(stl2step_lawband_dump PRIVATE cxx_std_17)
set_target_properties(stl2step_lawband_dump PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(stl2step_lawband_dump PRIVATE stl2step::harness)
target_include_directories(stl2step_lawband_dump PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/tests/harness")
if(MSVC)
    target_compile_definitions(stl2step_lawband_dump PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(stl2step_lawband_dump)
endif()
unset(_stl2step_root)
