# tests/gates/regiondump — RegionSet JSON emitter (lane p1-dump).
#
# Nested (after p0-gates lands tests/gates/CMakeLists.txt):
#   foreach(census regiondump buildtest) add_subdirectory(...)
#
# Standalone (this worktree, no tests/gates/CMakeLists.txt yet):
#   cmake -S tests/gates/regiondump -B build-dump -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-dump -j
#
# Depth from this directory to the repo root is ../../..
# The EXISTS guard on src/refit_segment.cpp compiles stub_segment.cpp only while
# that TU is absent; when it lands, this file does not need to change.

cmake_minimum_required(VERSION 3.16)

get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)

if(NOT TARGET stl2step::core)
    project(stl2step_regiondump LANGUAGES CXX)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL
            "regiondump standalone: do not recurse into tests/" FORCE)
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

# Reconfigure when any P1 TU appears or disappears so the EXISTS switch below
# is not stale across a no-edit drop-in.
set(_p1_tus
    refit_math.cpp
    refit_grow.cpp
    refit_fillet.cpp
    refit_chains.cpp
    refit_segment.cpp)
foreach(_tu IN LISTS _p1_tus)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${_stl2step_root}/src/${_tu}")
endforeach()

set(_regiondump_sources dump_regionset.cpp)
# One predicate over the five-name list: compile every P1 TU iff all five
# exist. A 1–4 half-link (FINDINGS-0 D1) is the defect this lane repairs.
set(_p1_all_present TRUE)
set(_p1_sources "")
foreach(_tu IN LISTS _p1_tus)
    set(_path "${_stl2step_root}/src/${_tu}")
    if(EXISTS "${_path}")
        list(APPEND _p1_sources "${_path}")
    else()
        set(_p1_all_present FALSE)
    endif()
endforeach()

if(_p1_all_present)
    list(APPEND _regiondump_sources ${_p1_sources})
    message(STATUS
        "stl2step_regiondump: using 5/5 P1 TUs (stub not compiled)")
else()
    list(APPEND _regiondump_sources stub_segment.cpp)
    message(WARNING
        "stl2step_regiondump: not all five P1 TUs present; compiling stub_segment.cpp")
endif()

add_executable(stl2step_regiondump ${_regiondump_sources})
target_compile_features(stl2step_regiondump PRIVATE cxx_std_17)
set_target_properties(stl2step_regiondump PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(stl2step_regiondump PRIVATE stl2step::harness)
target_include_directories(stl2step_regiondump PRIVATE
    "${_stl2step_root}/src")

if(NOT _p1_all_present)
    target_compile_definitions(stl2step_regiondump PRIVATE STL2STEP_REGIONDUMP_STUB=1)
endif()

if(MSVC)
    target_compile_definitions(stl2step_regiondump PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()

# CTest hook lives here (this file is owned). Do not edit tests/CMakeLists.txt.
if(EXISTS "${_stl2step_root}/tests/p1_compose/CMakeLists.txt")
    add_subdirectory("${_stl2step_root}/tests/p1_compose"
                     "${CMAKE_CURRENT_BINARY_DIR}/_p1_compose")
endif()

unset(_regiondump_sources)
unset(_p1_tus)
unset(_p1_sources)
unset(_p1_all_present)
unset(_tu)
unset(_path)
unset(_stl2step_root)

# ad-hoc signing (see root helper)
if(TARGET stl2step_regiondump)
  stl2step_adhoc_sign(stl2step_regiondump)
endif()
