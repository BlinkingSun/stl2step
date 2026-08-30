# tests/unit — law-band API tests (AC2-L1).
cmake_minimum_required(VERSION 3.16)

get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)

if(NOT TARGET stl2step::core)
    project(stl2step_unit LANGUAGES CXX)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
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

add_executable(stl2step_lawband_test lawband_test.cpp)
target_compile_features(stl2step_lawband_test PRIVATE cxx_std_17)
set_target_properties(stl2step_lawband_test PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(stl2step_lawband_test PRIVATE stl2step::harness)
target_include_directories(stl2step_lawband_test PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/tests/harness")
if(MSVC)
    target_compile_definitions(stl2step_lawband_test PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(stl2step_lawband_test)
endif()

set(_anat "${_stl2step_root}/_team/reports/ac2/band-anatomy.json")
if(NOT EXISTS "${_anat}")
    set(_anat "/Users/jroberts/Desktop/Internal Development/Tools/stl2step/_team/reports/ac2/band-anatomy.json")
endif()
set(_stl "${_stl2step_root}/tests/corpus/handle-lock.stl")
set(_mis "${_stl2step_root}/_team/reports/ac2/misfit-map.json")
if(NOT EXISTS "${_mis}")
    set(_mis "/Users/jroberts/Desktop/Internal Development/Tools/stl2step/_team/reports/ac2/misfit-map.json")
endif()
set(_lab "${_stl2step_root}/_team/reports/ac2/tri-labels.json")
if(NOT EXISTS "${_lab}")
    set(_lab "/Users/jroberts/Desktop/Internal Development/Tools/stl2step/_team/reports/ac2/tri-labels.json")
endif()

add_test(NAME lawband_unit
    COMMAND stl2step_lawband_test
        "${_stl}" "${_anat}" "${_lab}" "${_mis}")
set_tests_properties(lawband_unit PROPERTIES
    LABELS "unit;lawband"
    TIMEOUT 180)

# A16: RULE 4.2a — banned recognition literals must not appear as thresholds.
add_test(NAME lawband_rule42a_grep
    COMMAND "${CMAKE_COMMAND}" -E env
        "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_rule42a.cmake")
set_tests_properties(lawband_rule42a_grep PROPERTIES
    LABELS "unit;lawband"
    TIMEOUT 30)

add_executable(stl2step_prism_test prism_test.cpp)
target_compile_features(stl2step_prism_test PRIVATE cxx_std_17)
set_target_properties(stl2step_prism_test PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(stl2step_prism_test PRIVATE stl2step::harness)
target_include_directories(stl2step_prism_test PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/tests/harness")
if(MSVC)
    target_compile_definitions(stl2step_prism_test PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(stl2step_prism_test)
endif()

set(_body11 "${_stl2step_root}/tests/corpus/Body11.stl")
set(_neg "${_stl2step_root}/tests/corpus/nonprismatic-control.stl")
add_test(NAME prism_unit
    COMMAND stl2step_prism_test "${_stl}" "${_body11}" "${_neg}")
set_tests_properties(prism_unit PROPERTIES
    LABELS "unit;prism"
    TIMEOUT 180)

add_test(NAME prism_literal_ban
    COMMAND "${CMAKE_COMMAND}" -E env
        "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_prism_literals.cmake")
set_tests_properties(prism_literal_ban PROPERTIES
    LABELS "unit;prism"
    TIMEOUT 30)

add_executable(stl2step_profile_test
    profile_test.cpp
    "${_stl2step_root}/tests/tools/profile_dump/prism_levels_stub.cpp")
target_compile_features(stl2step_profile_test PRIVATE cxx_std_17)
set_target_properties(stl2step_profile_test PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(stl2step_profile_test PRIVATE stl2step::harness)
target_include_directories(stl2step_profile_test PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/tests/harness")
if(MSVC)
    target_compile_definitions(stl2step_profile_test PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(stl2step_profile_test)
endif()

add_test(NAME profile_unit
    COMMAND stl2step_profile_test "${_stl}")
set_tests_properties(profile_unit PROPERTIES
    LABELS "unit;profile"
    TIMEOUT 180)

add_executable(stl2step_prism_build_test prism_build_test.cpp)
target_compile_features(stl2step_prism_build_test PRIVATE cxx_std_17)
set_target_properties(stl2step_prism_build_test PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(stl2step_prism_build_test PRIVATE stl2step::harness stl2step::core)
target_include_directories(stl2step_prism_build_test PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/tests/harness")
if(MSVC)
    target_compile_definitions(stl2step_prism_build_test PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(stl2step_prism_build_test)
endif()

add_test(NAME prism_build_unit
    COMMAND stl2step_prism_build_test "${_stl}")
set_tests_properties(prism_build_unit PROPERTIES
    LABELS "unit;prism-build"
    TIMEOUT 180)

unset(_stl2step_root)
unset(_anat)
unset(_stl)
unset(_mis)
unset(_lab)
unset(_body11)
unset(_neg)
