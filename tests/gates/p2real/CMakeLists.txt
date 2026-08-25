# tests/gates/p2real — live RegionSet re-gate of refit::buildFaces.
# Picked up from tests/gates/p2buildtest/CMakeLists.txt (this lane owns both).
# Do not edit root CMakeLists.txt.

cmake_minimum_required(VERSION 3.16)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)

set(_p2real_script "${CMAKE_CURRENT_SOURCE_DIR}/run_p2real.py")

if(TARGET stl2step_p2buildtest AND TARGET stl2step_regiondump)
    add_test(NAME p2real_live
        COMMAND ${Python3_EXECUTABLE} ${_p2real_script}
                --repo ${_stl2step_root}
                --build-dir ${CMAKE_BINARY_DIR}
                --dump ${CMAKE_BINARY_DIR}/stl2step_regiondump
                --p2 ${CMAKE_BINARY_DIR}/stl2step_p2buildtest
                --corpus ${_stl2step_root}/tests/corpus
                --cube ${_stl2step_root}/tests/cube.stl
                --dumps-dir ${_stl2step_root}/_team/p1-dumps
                --sidecars ${_stl2step_root}/../p0-gen-rw1/tests/corpus
                --compose-dump ${_stl2step_root}/../p1-compose-fix-rw1/build/stl2step_regiondump)
    set_tests_properties(p2real_live PROPERTIES
        LABELS "gates;p2real"
        TIMEOUT 600)
endif()

unset(_stl2step_root)
unset(_p2real_script)
