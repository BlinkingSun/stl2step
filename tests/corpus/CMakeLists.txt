# tests/corpus — synthetic STL corpus + ground-truth sidecars (lane p0-gen).
#
# Standalone:
#   cmake -S tests/corpus -B build-corpus -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-corpus -j
# Nested: picked up by tests/CMakeLists.txt once stl2step::core exists.

cmake_minimum_required(VERSION 3.16)

if(NOT TARGET stl2step::core)
    project(stl2step_corpus LANGUAGES CXX)
    get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        set(STL2STEP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(STL2STEP_INSTALL OFF CACHE BOOL "" FORCE)
        add_subdirectory("${_stl2step_root}" "${CMAKE_CURRENT_BINARY_DIR}/_engine")
    else()
        find_package(stl2step REQUIRED)
    endif()
    unset(_stl2step_root)
endif()

set(_corpus_out_dir "${CMAKE_CURRENT_SOURCE_DIR}")
set(_corpus_stamp "${CMAKE_CURRENT_BINARY_DIR}/corpus.stamp")

add_executable(stl2step_corpus_gen gen_corpus.cpp)
target_compile_features(stl2step_corpus_gen PRIVATE cxx_std_17)
target_link_libraries(stl2step_corpus_gen PRIVATE stl2step::core)
if(MSVC)
    target_compile_definitions(stl2step_corpus_gen PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
# Generator uses MakeFillet/MakeChamfer (TKFillet), ThruSections (TKOffset),
# and booleans (TKBool). Those are NOT in the engine's OpenCASCADE::OpenCASCADE
# toolkit list. Resolve them in THIS directory: nested builds see the root
# find_package imported targets; standalone (`cmake -S tests/corpus`) does not,
# because that find_package ran inside add_subdirectory(_engine).
# Graft: re-run find_package HERE so OpenCASCADE_LIBRARY_DIR is populated in
# this directory's scope (standalone) before the find_library fallback.
if(NOT TARGET TKFillet OR NOT TARGET TKOffset OR NOT TARGET TKBool OR NOT TARGET TKMesh)
    find_package(OpenCASCADE REQUIRED)
endif()
# TKMesh: BRepMesh_IncrementalMesh in tessellate(). Homebrew's OCCT config
# links it transitively; Ubuntu's libocct 7.6 does not — name it explicitly.
foreach(_tk IN ITEMS TKFillet TKOffset TKBool TKMesh)
    if(TARGET ${_tk})
        target_link_libraries(stl2step_corpus_gen PRIVATE ${_tk})
    else()
        find_library(_corpus_${_tk}
            NAMES ${_tk}
            HINTS
                "${OpenCASCADE_LIBRARY_DIR}"
                "${OpenCASCADE_BINARY_DIR}"
                "${OpenCASCADE_ROOT}"
                "$ENV{OpenCASCADE_ROOT}"
                "$ENV{CASROOT}"
                ${CASROOT}
            PATH_SUFFIXES lib lib64 win64/vc14/lib
            PATHS
                /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
                /opt/homebrew/lib /opt/local/lib)
        if(NOT _corpus_${_tk})
            message(FATAL_ERROR
                "stl2step_corpus_gen requires OCCT toolkit ${_tk} (imported target or library)")
        endif()
        target_link_libraries(stl2step_corpus_gen PRIVATE "${_corpus_${_tk}}")
    endif()
endforeach()
set_target_properties(stl2step_corpus_gen PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_executable(stl2step_corpus_check check_corpus.cpp)
target_compile_features(stl2step_corpus_check PRIVATE cxx_std_17)
target_link_libraries(stl2step_corpus_check PRIVATE stl2step::core)
set_target_properties(stl2step_corpus_check PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_custom_command(
    OUTPUT "${_corpus_stamp}"
    COMMAND stl2step_corpus_gen "${_corpus_out_dir}"
    COMMAND "${CMAKE_COMMAND}" -E touch "${_corpus_stamp}"
    DEPENDS stl2step_corpus_gen gen_corpus.cpp corpus_util.hpp
    COMMENT "Generating synthetic STL corpus fixtures"
    VERBATIM)

add_custom_target(stl2step_corpus_data ALL DEPENDS "${_corpus_stamp}")

enable_testing()

add_test(NAME corpus_generate
    COMMAND stl2step_corpus_gen "${_corpus_out_dir}")
add_test(NAME corpus_sidecar_check
    COMMAND stl2step_corpus_check "${_corpus_out_dir}")
set_tests_properties(corpus_sidecar_check PROPERTIES
    DEPENDS corpus_generate)

add_test(NAME corpus_determinism
    COMMAND "${CMAKE_COMMAND}"
        -DGEN=$<TARGET_FILE:stl2step_corpus_gen>
        -DDIR_A=${CMAKE_CURRENT_BINARY_DIR}/_det_a
        -DDIR_B=${CMAKE_CURRENT_BINARY_DIR}/_det_b
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_corpus_determinism.cmake")

    if(TARGET stl2step)
    # A bare .sh has no shebang execution on Windows (ctest "Not Run");
    # invoke through bash when one exists (Git Bash on GitHub runners).
    find_program(STL2STEP_BASH bash)
    if(STL2STEP_BASH)
        set(_engine_check_cmd "${STL2STEP_BASH}"
            "${CMAKE_CURRENT_SOURCE_DIR}/run_engine_check.sh")
    else()
        set(_engine_check_cmd "${CMAKE_CURRENT_SOURCE_DIR}/run_engine_check.sh")
    endif()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(NAME corpus_engine_convert
        COMMAND ${_engine_check_cmd}
                $<TARGET_FILE:stl2step>)
    # PYTHON: Git Bash `python3` is the Windows Store stub (Permission denied).
    # PYTHONUTF8=1: CMake ENVIRONMENT (not `set X=1 &&`, which poisons the value).
    set_tests_properties(corpus_engine_convert PROPERTIES
        ENVIRONMENT "CORPUS=${_corpus_out_dir};CENSUS=$<TARGET_FILE:stl2step_census>;PYTHON=${Python3_EXECUTABLE};PYTHONUTF8=1"
        DEPENDS corpus_generate)
endif()
