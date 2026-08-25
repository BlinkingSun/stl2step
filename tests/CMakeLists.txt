# stl2step test tree.
#
# Deliberately near-empty: tools arrive one subdirectory at a time (corpus
# generator, STEP census, RegionSet dump, build-fixture harness, mesh harness).
# Each of those lanes adds its OWN CMakeLists.txt under its OWN directory and
# is picked up by the guarded blocks below - nothing else in the build system
# needs to change when a tool lands.
#
# Conventions for tools added here (binding on later lanes):
#   * one directory per tool, e.g. tests/harness/, tests/corpus/, tests/gates/
#   * link stl2step::core (never the raw stl2step_core name)
#   * C++17, no new third-party dependencies
#   * EXCLUDE_FROM_ALL is NOT used: a tool that does not build is a red build

foreach(_stl2step_test_dir harness corpus gates)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_stl2step_test_dir}/CMakeLists.txt")
        add_subdirectory(${_stl2step_test_dir})
    endif()
endforeach()
unset(_stl2step_test_dir)
