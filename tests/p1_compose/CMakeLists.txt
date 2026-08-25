cmake_minimum_required(VERSION 3.16)
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)

if(NOT TARGET stl2step::harness)
    add_subdirectory("${_root}/tests/harness" "${CMAKE_CURRENT_BINARY_DIR}/_harness")
endif()

set(_p1_tus
    refit_math.cpp refit_grow.cpp refit_fillet.cpp refit_chains.cpp refit_segment.cpp)
set(_p1_sources "")
foreach(_tu IN LISTS _p1_tus)
    list(APPEND _p1_sources "${_root}/src/${_tu}")
endforeach()

add_executable(p1_compose_diag diag_stages.cpp ${_p1_sources})
target_compile_features(p1_compose_diag PRIVATE cxx_std_17)
set_target_properties(p1_compose_diag PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(p1_compose_diag PRIVATE stl2step::harness)
target_include_directories(p1_compose_diag PRIVATE "${_root}/src")
if(MSVC)
    target_compile_definitions(p1_compose_diag PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()

find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(_compose_script "${CMAKE_CURRENT_SOURCE_DIR}/run_compose_gate.py")
add_test(NAME p1_compose_gate
    COMMAND ${Python3_EXECUTABLE} ${_compose_script}
            --repo ${_root}
            --build-dir ${CMAKE_BINARY_DIR}
            --dump ${CMAKE_BINARY_DIR}/stl2step_regiondump
            --ichecker ${_root}/tests/gates/check_regionset.py
            --corpus ${_root}/tests/corpus
            --cube ${_root}/tests/cube.stl)
set_tests_properties(p1_compose_gate PROPERTIES
    LABELS "p1_compose"
    TIMEOUT 600)

# AC2-A1: ctest -E gates_full is the green suite; this guard asserts the
# inherited gates_full G1 list has not grown.
add_test(NAME p1_ac2_gates_baseline
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_gates_full_baseline.py
            --build-dir ${CMAKE_BINARY_DIR})
set_tests_properties(p1_ac2_gates_baseline PROPERTIES
    LABELS "p1_compose"
    TIMEOUT 600)

unset(_root)
unset(_p1_tus)
unset(_p1_sources)
unset(_tu)
unset(_compose_script)
