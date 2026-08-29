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

# SPEC AC-J6HEAL API probe — calls stl2step::convert() directly. Not a ctest
# (lane AC pins 27/27). Built whenever the core library is.
if(TARGET stl2step::core)
    add_executable(j6heal_api_probe j6heal_api_probe.cpp)
    target_compile_features(j6heal_api_probe PRIVATE cxx_std_17)
    target_link_libraries(j6heal_api_probe PRIVATE stl2step::core)
    if(COMMAND stl2step_adhoc_sign)
        stl2step_adhoc_sign(j6heal_api_probe)
    endif()
endif()
