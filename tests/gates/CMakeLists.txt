# tests/gates — shared gate tree (lane p0-gates owns this file only).
#
# Nested (after tests/CMakeLists.txt lands):
#   add_subdirectory(gates) picks up census/, regiondump/, buildtest/ as they land.
#
# Standalone:
#   cmake -S tests/gates -B build-gates
#   cmake --build build-gates -j

cmake_minimum_required(VERSION 3.16)

if(NOT TARGET stl2step::core)
    project(stl2step_gates LANGUAGES CXX)
    get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL
            "Gates standalone: do not recurse into tests/" FORCE)
        set(STL2STEP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(STL2STEP_INSTALL OFF CACHE BOOL "" FORCE)
        add_subdirectory("${_stl2step_root}" "${CMAKE_CURRENT_BINARY_DIR}/_engine")
    else()
        find_package(stl2step REQUIRED)
    endif()
    unset(_stl2step_root)
endif()

foreach(_stl2step_gate_dir census regiondump buildtest)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_stl2step_gate_dir}/CMakeLists.txt")
        add_subdirectory(${_stl2step_gate_dir})
    endif()
endforeach()
unset(_stl2step_gate_dir)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(_stl2step_gates_runner "${CMAKE_CURRENT_SOURCE_DIR}/run_gates.py")
set(_stl2step_gates_cube "${CMAKE_CURRENT_SOURCE_DIR}/../cube.stl")
set(_stl2step_gates_corpus "${CMAKE_CURRENT_SOURCE_DIR}/../corpus")

if(EXISTS "${_stl2step_gates_corpus}/S01.stl" AND TARGET stl2step)
    add_test(NAME gates_full
        COMMAND "${Python3_EXECUTABLE}" "${_stl2step_gates_runner}"
            --binary "$<TARGET_FILE:stl2step>"
            --corpus "${_stl2step_gates_corpus}"
            --no-verify
            --jobs 4)
    set_tests_properties(gates_full PROPERTIES
        LABELS "gates"
        TIMEOUT 600)
endif()

if(TARGET stl2step AND EXISTS "${_stl2step_gates_runner}" AND EXISTS "${_stl2step_gates_cube}")
    add_test(NAME gates_smoke
        COMMAND "${Python3_EXECUTABLE}" "${_stl2step_gates_runner}"
            --smoke
            --binary "$<TARGET_FILE:stl2step>"
            --no-verify)
    set_tests_properties(gates_smoke PROPERTIES
        LABELS "gates;smoke"
        TIMEOUT 120)
endif()

unset(_stl2step_gates_runner)
unset(_stl2step_gates_cube)
unset(_stl2step_gates_corpus)
