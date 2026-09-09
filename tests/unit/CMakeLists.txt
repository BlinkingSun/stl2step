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

set(_anat "${_stl2step_root}/tests/gates/labels/handle-lock.band-anatomy.json")
set(_stl "${_stl2step_root}/tests/corpus/handle-lock.stl")
set(_mis "${_stl2step_root}/tests/gates/labels/handle-lock.misfit-map.json")
set(_lab "${_stl2step_root}/tests/gates/labels/handle-lock.tri-labels.json")

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

# 130-CONE-MATH (SPEC-130-cone-math) — the closed-form cone math in
# src/refit_cone_math.cpp, certified against the OCCT parametric equation,
# GeomAPI_ProjectPointOnSurf and IntAna_QuadQuadGeo.
#
# This target deliberately COMPILES src/refit_cone_math.cpp instead of linking
# stl2step::core: the math under test is pure gp_ geometry (refit_cone_math.hpp
# pulls in no engine header), and `--target cone_math_unit` must not drag in
# refit_build.cpp, which another lane owns. That one file is the whole TU — the
# cone math sits outside src/refit_math.cpp precisely because the D5.3 include
# allowlist covers the P1 sources and admits neither gp.hxx nor a third project
# header. Nothing else in the tree changes how it builds.
# `OpenCASCADE::OpenCASCADE` is an IMPORTED target and therefore DIRECTORY
# scoped: the sibling targets only see it transitively through stl2step::core.
# When this file is configured standalone (the `NOT TARGET stl2step::core`
# branch at the top) the engine's find_package ran in a child scope and the name
# does not exist here, so naming it directly fails at generate time. Re-find it
# in this scope — the module path and the cache are already primed, so it costs
# nothing — rather than linking stl2step::core, which carries its own
# refit_cone_math.cpp and would duplicate every symbol in this executable.
if(NOT TARGET OpenCASCADE::OpenCASCADE)
    list(APPEND CMAKE_MODULE_PATH "${_stl2step_root}/cmake")
    find_package(OpenCASCADE REQUIRED)
endif()

add_executable(cone_math_unit
    cone_math_test.cpp
    "${_stl2step_root}/src/refit_cone_math.cpp")
target_compile_features(cone_math_unit PRIVATE cxx_std_17)
set_target_properties(cone_math_unit PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_include_directories(cone_math_unit PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/include")
target_link_libraries(cone_math_unit PRIVATE OpenCASCADE::OpenCASCADE)
if(MSVC)
    target_compile_definitions(cone_math_unit PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(cone_math_unit PRIVATE /EHsc)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(cone_math_unit)
endif()

add_test(NAME cone_math_unit COMMAND cone_math_unit)
set_tests_properties(cone_math_unit PROPERTIES
    LABELS "unit;cone-math"
    TIMEOUT 60)

# 130-CONE-BUILD — the certified `circle-on-cone` / `line-on-cone` BIND classes
# in src/refit_cone_bind.cpp: the supremum of |C3d(t) - S(pcurve(t))|, which is
# the surface supremum cone_math_unit certifies COMPOSED with the
# parametrisation. Same construction as cone_math_unit for the same reason: it
# compiles the two cone TUs rather than linking stl2step::core, so
# `--target cone_bind_unit` never drags in refit_build.cpp.
add_executable(cone_bind_unit
    cone_bind_test.cpp
    "${_stl2step_root}/src/refit_cone_bind.cpp"
    "${_stl2step_root}/src/refit_cone_math.cpp")
target_compile_features(cone_bind_unit PRIVATE cxx_std_17)
set_target_properties(cone_bind_unit PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_include_directories(cone_bind_unit PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/include")
target_link_libraries(cone_bind_unit PRIVATE OpenCASCADE::OpenCASCADE)
if(MSVC)
    target_compile_definitions(cone_bind_unit PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(cone_bind_unit PRIVATE /EHsc)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(cone_bind_unit)
endif()

add_test(NAME cone_bind_unit COMMAND cone_bind_unit)
set_tests_properties(cone_bind_unit PROPERTIES
    LABELS "unit;cone-math"
    TIMEOUT 60)

# torus-math (SPEC-torus PART B, D-140-6 §3(4)) — the closed-form
# `circle-on-torus` supremum in src/refit_torus_math.cpp, certified against a
# numeric distance oracle built from the OCCT parametric equation (meridian
# circles), GeomAPI_ProjectPointOnSurf on a Geom_ToroidalSurface, and the S04 /
# S19 tangency circles. Same construction as cone_math_unit for the same reason:
# it COMPILES the one pure-gp_ TU instead of linking stl2step::core, so
# `--target torus_math_unit` never drags in refit_build.cpp.
add_executable(torus_math_unit
    torus_math_test.cpp
    "${_stl2step_root}/src/refit_torus_math.cpp")
target_compile_features(torus_math_unit PRIVATE cxx_std_17)
set_target_properties(torus_math_unit PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_include_directories(torus_math_unit PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/include")
target_link_libraries(torus_math_unit PRIVATE OpenCASCADE::OpenCASCADE)
if(MSVC)
    target_compile_definitions(torus_math_unit PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(torus_math_unit PRIVATE /EHsc)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(torus_math_unit)
endif()

add_test(NAME torus_math_unit COMMAND torus_math_unit)
set_tests_properties(torus_math_unit PROPERTIES
    LABELS "unit;torus-math"
    TIMEOUT 120)

# ellipse-math (SPEC-ellipse PART A, D-140-9 §2) — the closed forms behind the
# ELLIPSE bind class in src/refit_ellipse_math.cpp: the EXACT `ellipse-on-plane`
# supremum (quartic critical set) and the certified `ellipse-on-cylinder-bound`
# (section identity + Bernstein hull + Chebyshev remainder), certified against a
# >= 10^4-sample refined oracle on OCCT's own pcurves (GeomAPI::To2d,
# GeomProjLib::Curve2d) and S11-b's measured degree-8 poles. Same construction as
# cone_math_unit for the same reason: it COMPILES the one pure-gp_ TU instead of
# linking stl2step::core, so `--target ellipse_math_unit` never drags in
# refit_build.cpp.
add_executable(ellipse_math_unit
    ellipse_math_test.cpp
    "${_stl2step_root}/src/refit_ellipse_math.cpp")
target_compile_features(ellipse_math_unit PRIVATE cxx_std_17)
set_target_properties(ellipse_math_unit PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_include_directories(ellipse_math_unit PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/include")
target_link_libraries(ellipse_math_unit PRIVATE OpenCASCADE::OpenCASCADE)
if(MSVC)
    target_compile_definitions(ellipse_math_unit PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(ellipse_math_unit PRIVATE /EHsc)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(ellipse_math_unit)
endif()

add_test(NAME ellipse_math_unit COMMAND ellipse_math_unit)
set_tests_properties(ellipse_math_unit PROPERTIES
    LABELS "unit;ellipse-math"
    TIMEOUT 300)

# 130-ABSORB-3 (D-130-12) — the mesh file's own coordinate quantization.
# Compiles src/stl_quant.cpp directly rather than linking stl2step::core: the TU
# is pure <cstdio>/<cmath> with no OCCT and no engine header, so the test needs
# no OpenCASCADE at all and `--target quant_floor_unit` drags in nothing else.
add_executable(quant_floor_unit
    quant_floor_test.cpp
    "${_stl2step_root}/src/stl_quant.cpp")
target_compile_features(quant_floor_unit PRIVATE cxx_std_17)
set_target_properties(quant_floor_unit PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_include_directories(quant_floor_unit PRIVATE "${_stl2step_root}/src")
if(MSVC)
    target_compile_definitions(quant_floor_unit PRIVATE _USE_MATH_DEFINES NOMINMAX)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(quant_floor_unit)
endif()

add_test(NAME quant_floor_unit COMMAND quant_floor_unit "${CMAKE_BINARY_DIR}")
set_tests_properties(quant_floor_unit PROPERTIES
    LABELS "unit;quant"
    TIMEOUT 60)

# SPEC-volid §4.3 / D-140-10 — S01/S02/S03 volume-attribution identity.
# Spawns the CLI with STL2STEP_VOLID=1; asserts the §2.5 closesA/closesB
# bracket on the three synthetic fixtures only.
if(TARGET stl2step)
    add_executable(volid_identity_test volid_identity.cpp)
    target_compile_features(volid_identity_test PRIVATE cxx_std_17)
    set_target_properties(volid_identity_test PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
    target_include_directories(volid_identity_test PRIVATE
        "${_stl2step_root}/tests/support")
    if(MSVC)
        target_compile_definitions(volid_identity_test PRIVATE _USE_MATH_DEFINES NOMINMAX)
    endif()
    if(COMMAND stl2step_adhoc_sign)
        stl2step_adhoc_sign(volid_identity_test)
    endif()
    add_test(NAME volid_identity
        COMMAND volid_identity_test
            $<TARGET_FILE:stl2step>
            "${_stl2step_root}/tests/corpus/S01.stl"
            "${_stl2step_root}/tests/corpus/S02.stl"
            "${_stl2step_root}/tests/corpus/S03.stl")
    set_tests_properties(volid_identity PROPERTIES
        LABELS "unit;volid"
        TIMEOUT 180
        DEPENDS corpus_generate)
endif()

# 140-UNION-CENSUS (SPEC-union PART B, D-140-8 U-R5..U-R8, U-R13) -- the
# same-surface union census in src/refit_union_census.cpp: `claimChordSagitta`
# against a closed form and a numerical supremum along every welded edge, and
# the U-R6/U-R7/U-R8 counts on hand-built rings. Same construction as
# cone_math_unit for the same reason: it compiles the one TU rather than
# linking stl2step::core, so `--target union_census_unit` never drags in
# refit_grow.cpp or refit_build.cpp, which other lanes own.
add_executable(union_census_unit
    union_census_test.cpp
    "${_stl2step_root}/src/refit_union_census.cpp")
target_compile_features(union_census_unit PRIVATE cxx_std_17)
set_target_properties(union_census_unit PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_include_directories(union_census_unit PRIVATE
    "${_stl2step_root}/src"
    "${_stl2step_root}/include")
target_link_libraries(union_census_unit PRIVATE OpenCASCADE::OpenCASCADE)
if(MSVC)
    target_compile_definitions(union_census_unit PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_compile_options(union_census_unit PRIVATE /EHsc)
endif()
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(union_census_unit)
endif()

add_test(NAME union_census_unit COMMAND union_census_unit)
set_tests_properties(union_census_unit PROPERTIES
    LABELS "unit;union-census"
    TIMEOUT 60)

# D-142-1: --threads 1 constructs zero std::thread; threads=4 spawns and
# produces identical canonical STEP DATA.
add_executable(threads_spawn_test threads_spawn_test.cpp)
target_compile_features(threads_spawn_test PRIVATE cxx_std_17)
set_target_properties(threads_spawn_test PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(threads_spawn_test PRIVATE stl2step::core)
target_include_directories(threads_spawn_test PRIVATE
    "${_stl2step_root}/include")
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(threads_spawn_test)
endif()

set(_threads_scratch "${CMAKE_BINARY_DIR}/threads1-test")
set(_mesh_cyl "${_stl2step_root}/tests/mesh/seamed_cylinder.step")
add_test(NAME threads1_no_spawn
    COMMAND threads_spawn_test no_spawn
        "${_stl2step_root}/tests/corpus"
        "${_mesh_cyl}"
        "${_threads_scratch}")
set_tests_properties(threads1_no_spawn PROPERTIES
    LABELS "unit;threads"
    TIMEOUT 300)
add_test(NAME threads_default_spawns
    COMMAND threads_spawn_test default_spawns
        "${_stl2step_root}/tests/corpus"
        "${_mesh_cyl}"
        "${_threads_scratch}")
set_tests_properties(threads_default_spawns PROPERTIES
    LABELS "unit;threads"
    TIMEOUT 300)
unset(_threads_scratch)
unset(_mesh_cyl)

# D-I-4: Result::toJson must not wrap a negative snprintf count into a 0-size
# string (MSVC SSO 16). Portable: the wrap itself is the unit condition.
add_executable(result_tojson_test result_tojson_test.cpp)
target_compile_features(result_tojson_test PRIVATE cxx_std_17)
set_target_properties(result_tojson_test PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_link_libraries(result_tojson_test PRIVATE stl2step::core)
target_include_directories(result_tojson_test PRIVATE
    "${_stl2step_root}/include")
if(COMMAND stl2step_adhoc_sign)
    stl2step_adhoc_sign(result_tojson_test)
endif()
add_test(NAME result_tojson_unit COMMAND result_tojson_test)
set_tests_properties(result_tojson_unit PROPERTIES
    LABELS "unit;tojson"
    TIMEOUT 30)

unset(_stl2step_root)
unset(_anat)
unset(_stl)
unset(_mis)
unset(_lab)
unset(_body11)
unset(_neg)
