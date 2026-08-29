cmake_minimum_required(VERSION 3.16)

# Ad-hoc codesign built executables on macOS: a freshly built unsigned binary
# doing heavy test I/O trips SentinelOne behavioral remediation (2026-08-25).
function(stl2step_adhoc_sign tgt)
  if(APPLE)
    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND codesign --force --sign - $<TARGET_FILE:${tgt}>
      COMMENT "ad-hoc codesign ${tgt}" VERBATIM)
  endif()
endfunction()


project(stl2step
    VERSION 1.1.0
    DESCRIPTION "Universal STL (mesh) -> STEP (B-Rep solid) conversion engine"
    LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(STL2STEP_BUILD_CLI      "Build the stl2step command-line tool" ON)
option(STL2STEP_BUILD_EXAMPLES "Build the integration examples"       OFF)
option(STL2STEP_BUILD_TESTS    "Register the smoke test with CTest"
       ${PROJECT_IS_TOP_LEVEL})
option(STL2STEP_INSTALL        "Generate install + package-config rules"
       ${PROJECT_IS_TOP_LEVEL})

# Default to an optimised build when the caller did not pick one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Threads REQUIRED)
find_package(OpenCASCADE REQUIRED)

# ---------------------------------------------------------------------------
# Engine library
# ---------------------------------------------------------------------------
add_library(stl2step_core
    src/stl2step.cpp
    src/refit_math.cpp
    src/refit_grow.cpp
    src/refit_fillet.cpp
    src/refit_chains.cpp
    src/refit_segment.cpp
    src/refit_build.cpp
    src/refit_lawband.cpp)
add_library(stl2step::core ALIAS stl2step_core)

target_compile_features(stl2step_core PUBLIC cxx_std_17)
set_target_properties(stl2step_core PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON)

# On Windows an executable named `stl2step` also emits an import `stl2step.lib`,
# which collides with the static library's `stl2step.lib` (LNK1149). Keep the
# tidy `libstl2step.{a,so,dylib}` name on Unix; fall back to the target name
# (`stl2step_core.lib`) on Windows. CMake consumers link the `stl2step::core`
# target, so the on-disk name is transparent to them either way.
if(NOT WIN32)
    set_target_properties(stl2step_core PROPERTIES OUTPUT_NAME stl2step)
endif()

target_include_directories(stl2step_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(stl2step_core
    PUBLIC  OpenCASCADE::OpenCASCADE
    PRIVATE Threads::Threads)

if(MSVC)
    target_compile_definitions(stl2step_core PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(stl2step_core PRIVATE /EHsc)
endif()

# ---------------------------------------------------------------------------
# Command-line tool
# ---------------------------------------------------------------------------
if(STL2STEP_BUILD_CLI)
    add_executable(stl2step src/main.cpp)
    target_link_libraries(stl2step PRIVATE stl2step::core)
    set_target_properties(stl2step PROPERTIES OUTPUT_NAME stl2step)
endif()

# ---------------------------------------------------------------------------
# Smoke test (converts the committed cube and expects a clean exit 0)
# ---------------------------------------------------------------------------
if(STL2STEP_BUILD_TESTS AND STL2STEP_BUILD_CLI)
    enable_testing()
    add_test(NAME convert_cube
        COMMAND stl2step "${CMAKE_CURRENT_SOURCE_DIR}/tests/cube.stl"
                         "${CMAKE_CURRENT_BINARY_DIR}/cube.step")
    set_tests_properties(convert_cube PROPERTIES
        PASS_REGULAR_EXPRESSION "\"ok\":true.*\"solids\":1.*\"volumeDeltaPct\":0")
endif()

# ---------------------------------------------------------------------------
# Test / gate tree (corpus generator, censuses, region dumps, gate harnesses).
# Each tool owns its own CMakeLists so the lanes that add them never collide
# in this file.
# ---------------------------------------------------------------------------
if(STL2STEP_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------
if(STL2STEP_BUILD_EXAMPLES)
    add_executable(convert_min examples/convert_min.cpp)
    target_link_libraries(convert_min PRIVATE stl2step::core)
endif()

# ---------------------------------------------------------------------------
# Install + package config  (lets downstream projects find_package(stl2step))
# ---------------------------------------------------------------------------
if(STL2STEP_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    set(_targets stl2step_core)
    if(STL2STEP_BUILD_CLI)
        list(APPEND _targets stl2step)
    endif()

    install(TARGETS ${_targets}
        EXPORT stl2stepTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(EXPORT stl2stepTargets
        NAMESPACE stl2step::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/stl2step)

    # Ship the OCCT finder so downstream find_package(stl2step) can locate OCCT.
    install(FILES cmake/FindOpenCASCADE.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/stl2step)

    configure_package_config_file(
        cmake/stl2stepConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/stl2stepConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/stl2step)

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/stl2stepConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/stl2stepConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/stl2stepConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/stl2step)
endif()

# ad-hoc signing (see root helper)
if(TARGET stl2step)
  stl2step_adhoc_sign(stl2step)
endif()
if(TARGET convert_min)
  stl2step_adhoc_sign(convert_min)
endif()
