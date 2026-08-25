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

# Reconfigure when the real P1 TU appears or disappears so the EXISTS
# switch below is not stale across a no-edit drop-in.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${_stl2step_root}/src/refit_segment.cpp")
file(GLOB _refit_src_watch "${_stl2step_root}/src/refit_*.cpp")
if(_refit_src_watch)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_refit_src_watch})
endif()
unset(_refit_src_watch)

set(_regiondump_sources dump_regionset.cpp)
set(_refit_segment "${_stl2step_root}/src/refit_segment.cpp")
if(EXISTS "${_refit_segment}")
    list(APPEND _regiondump_sources "${_refit_segment}")
    message(STATUS "stl2step_regiondump: using ${_refit_segment} (stub not compiled)")
else()
    list(APPEND _regiondump_sources stub_segment.cpp)
    message(WARNING
        "stl2step_regiondump: src/refit_segment.cpp absent; compiling stub_segment.cpp (vanishes with no file edit when the real TU lands)")
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

if(NOT EXISTS "${_refit_segment}")
    target_compile_definitions(stl2step_regiondump PRIVATE STL2STEP_REGIONDUMP_STUB=1)
endif()

unset(_regiondump_sources)
unset(_refit_segment)
unset(_stl2step_root)
