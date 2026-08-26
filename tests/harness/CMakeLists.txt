# tests/harness — weld/split MeshView harness (lane p-meshview-harness).
#
# Nested (after p-cmake-scaffold lands):
#   tests/CMakeLists.txt does add_subdirectory(harness) once stl2step::core exists.
#
# Standalone (this worktree, no tests/CMakeLists.txt yet):
#   cmake -S tests/harness -B build-harness -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-harness -j
# The guard below pulls in the parent engine project when stl2step::core is
# missing, with STL2STEP_BUILD_TESTS=OFF so a later root add_subdirectory(tests)
# cannot recurse back into this file.

cmake_minimum_required(VERSION 3.16)

if(NOT TARGET stl2step::core)
    project(stl2step_meshharness LANGUAGES CXX)
    get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL
            "Harness standalone: do not recurse into tests/" FORCE)
        set(STL2STEP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(STL2STEP_INSTALL OFF CACHE BOOL "" FORCE)
        # CLI stays ON so the same build tree can run engine-vs-harness comparison.
        add_subdirectory("${_stl2step_root}" "${CMAKE_CURRENT_BINARY_DIR}/_engine")
    else()
        find_package(stl2step REQUIRED)
    endif()
    unset(_stl2step_root)
endif()

add_library(stl2step_harness STATIC mesh_harness.cpp)
add_library(stl2step::harness ALIAS stl2step_harness)
target_compile_features(stl2step_harness PUBLIC cxx_std_17)
set_target_properties(stl2step_harness PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF)
target_include_directories(stl2step_harness PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(stl2step_harness PUBLIC stl2step::core)

add_executable(stl2step_meshharness mesh_dump.cpp)
target_compile_features(stl2step_meshharness PRIVATE cxx_std_17)
target_link_libraries(stl2step_meshharness PRIVATE stl2step_harness)
set_target_properties(stl2step_meshharness PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_executable(stl2step_harness_test mesh_harness_test.cpp)
target_compile_features(stl2step_harness_test PRIVATE cxx_std_17)
target_link_libraries(stl2step_harness_test PRIVATE stl2step_harness)

enable_testing()

set(_harness_deg  "${CMAKE_CURRENT_SOURCE_DIR}/testdata/degenerate.stl")
set(_harness_two  "${CMAKE_CURRENT_SOURCE_DIR}/testdata/two_cubes.stl")
set(_harness_cube "${CMAKE_CURRENT_SOURCE_DIR}/../cube.stl")

add_test(NAME harness_fixtures
    COMMAND stl2step_harness_test
        "${_harness_deg}" "${_harness_two}" "${_harness_cube}")

add_test(NAME harness_dump_cube
    COMMAND stl2step_meshharness "${_harness_cube}")
set_tests_properties(harness_dump_cube PROPERTIES
    PASS_REGULAR_EXPRESSION
    "\"triangles\":12.*\"vertices\":8.*\"tris\":12.*\"vtx\":8.*\"edges\":18.*\"open\":0.*\"conflict\":0.*\"nonManifold\":0.*\"clean\":true")

add_test(NAME harness_determinism
    COMMAND "${CMAKE_COMMAND}"
        -DDUMP=$<TARGET_FILE:stl2step_meshharness>
        -DFILE1=${_harness_cube}
        -DFILE2=${_harness_deg}
        -DFILE3=${_harness_two}
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_determinism.cmake")

unset(_harness_deg)
unset(_harness_two)
unset(_harness_cube)

# ad-hoc signing (see root helper)
if(TARGET stl2step_meshharness)
  stl2step_adhoc_sign(stl2step_meshharness)
endif()
if(TARGET stl2step_harness_test)
  stl2step_adhoc_sign(stl2step_harness_test)
endif()
