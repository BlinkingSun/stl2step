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

foreach(_stl2step_test_dir harness corpus gates unit)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_stl2step_test_dir}/CMakeLists.txt")
        add_subdirectory(${_stl2step_test_dir})
    endif()
endforeach()
unset(_stl2step_test_dir)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/lawband_dump/CMakeLists.txt")
    add_subdirectory(tools/lawband_dump)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/dxf_check/CMakeLists.txt")
    add_subdirectory(tools/dxf_check)
endif()

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

# SPEC AC-CASCADE API probe — calls stl2step::convert() directly. Not a ctest
# (lane pins 27/27 + hl_census_ratchet + cascade_no_waiver).
if(TARGET stl2step::core)
    add_executable(cascade_api_probe cascade_api_probe.cpp)
    target_compile_features(cascade_api_probe PRIVATE cxx_std_17)
    target_link_libraries(cascade_api_probe PRIVATE stl2step::core)
    if(COMMAND stl2step_adhoc_sign)
        stl2step_adhoc_sign(cascade_api_probe)
    endif()
endif()

# SPEC AC2-L3 RULE 2.2 scoped-adoption probe — convert() + C4/C5/C6. Not a ctest
# (C8 pins the clean count).
if(TARGET stl2step::core)
    add_executable(rule22_adopt_api rule22_adopt_api.cpp)
    target_compile_features(rule22_adopt_api PRIVATE cxx_std_17)
    target_link_libraries(rule22_adopt_api PRIVATE stl2step::core)
    if(COMMAND stl2step_adhoc_sign)
        stl2step_adhoc_sign(rule22_adopt_api)
    endif()
endif()
