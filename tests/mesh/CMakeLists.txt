# tests/mesh — --mesh mode ctest fixtures (U3a-2).
#
# Committed: seamed_cylinder.step (R=10 mm, H=30 mm analytic cylinder).
# Generator source: gen_seamed_cylinder.cpp (not built by default).

if(NOT TARGET stl2step)
    return()
endif()

find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(_mesh_py "${CMAKE_CURRENT_SOURCE_DIR}/mesh_tests.py")
set(_mesh_cyl "${CMAKE_CURRENT_SOURCE_DIR}/seamed_cylinder.step")
set(_mesh_step "${CMAKE_CURRENT_SOURCE_DIR}/../corpus/handle-lock.step")
set(_mesh_stl "${CMAKE_CURRENT_SOURCE_DIR}/../corpus/handle-lock.stl")

if(NOT EXISTS "${_mesh_py}")
    return()
endif()

add_test(NAME mesh_step_smoke
    COMMAND ${Python3_EXECUTABLE} "${_mesh_py}"
        --binary "$<TARGET_FILE:stl2step>"
        --test smoke
        --step "${_mesh_step}"
        --make-step-from "${_mesh_stl}")
set_tests_properties(mesh_step_smoke PROPERTIES
    LABELS "mesh"
    TIMEOUT 180
    REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_mesh_py};${_mesh_stl}")

add_test(NAME mesh_edges_seam
    COMMAND ${Python3_EXECUTABLE} "${_mesh_py}"
        --binary "$<TARGET_FILE:stl2step>"
        --test edges_seam
        --cylinder-step "${_mesh_cyl}")
set_tests_properties(mesh_edges_seam PROPERTIES
    LABELS "mesh"
    TIMEOUT 60
    REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_mesh_py};${_mesh_cyl}")

add_test(NAME mesh_fail_closed
    COMMAND ${Python3_EXECUTABLE} "${_mesh_py}"
        --binary "$<TARGET_FILE:stl2step>"
        --test fail_closed
        --step "${_mesh_step}"
        --make-step-from "${_mesh_stl}")
set_tests_properties(mesh_fail_closed PROPERTIES
    LABELS "mesh"
    TIMEOUT 60
    REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_mesh_py};${_mesh_stl}")

unset(_mesh_py)
unset(_mesh_cyl)
unset(_mesh_step)
unset(_mesh_stl)
