# tests/gates — shared gate tree (lane p0-gates owns this file only).
#
# Nested (after tests/CMakeLists.txt lands):
#   add_subdirectory(gates) picks up census/, regiondump/, p2buildtest/ as they land.
#
# Standalone:
#   cmake -S tests/gates -B build-gates
#   cmake --build build-gates -j

cmake_minimum_required(VERSION 3.16)

if(NOT TARGET stl2step::core)
    project(stl2step_gates LANGUAGES CXX)
    get_filename_component(_stl2step_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    if(EXISTS "${_stl2step_root}/CMakeLists.txt")
        set(STL2STEP_BUILD_TESTS OFF CACHE BOOL
            "Gates standalone: do not recurse into tests/" FORCE)
        set(STL2STEP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(STL2STEP_INSTALL OFF CACHE BOOL "" FORCE)
        add_subdirectory("${_stl2step_root}" "${CMAKE_CURRENT_BINARY_DIR}/_engine")
    else()
        find_package(stl2step REQUIRED)
    endif()
    unset(_stl2step_root)
endif()

foreach(_stl2step_gate_dir census regiondump buildtest p2buildtest)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_stl2step_gate_dir}/CMakeLists.txt")
        add_subdirectory(${_stl2step_gate_dir})
    endif()
endforeach()
unset(_stl2step_gate_dir)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(_stl2step_gates_runner "${CMAKE_CURRENT_SOURCE_DIR}/run_gates.py")
set(_stl2step_gates_cube "${CMAKE_CURRENT_SOURCE_DIR}/../cube.stl")
set(_stl2step_gates_corpus "${CMAKE_CURRENT_SOURCE_DIR}/../corpus")

# CTest default path: pass --census when the witness is in this build so G1
# uses BRepCheck, not the watertight proxy. DEPENDS corpus_generate so the
# STLs exist even if they were absent at configure time.
if(TARGET stl2step AND EXISTS "${_stl2step_gates_runner}")
    set(_gates_full_cmd
        "${Python3_EXECUTABLE}" "${_stl2step_gates_runner}"
        --binary "$<TARGET_FILE:stl2step>"
        --corpus "${_stl2step_gates_corpus}"
        --no-verify
        --jobs 4)
    set(_gates_full_depends "")
    set(_gates_full_required "$<TARGET_FILE:stl2step>")
    if(TARGET stl2step_census)
        list(APPEND _gates_full_cmd --census "$<TARGET_FILE:stl2step_census>")
        set(_gates_full_required "${_gates_full_required};$<TARGET_FILE:stl2step_census>")
    endif()
    if(TARGET stl2step_regiondump)
        list(APPEND _gates_full_cmd --dump "$<TARGET_FILE:stl2step_regiondump>")
        set(_gates_full_required "${_gates_full_required};$<TARGET_FILE:stl2step_regiondump>")
    endif()
    add_test(NAME gates_full COMMAND ${_gates_full_cmd})
    # corpus_generate lives in the sibling tests/corpus/ dir; if(TEST) is not
    # visible here, but CTest resolves DEPENDS by name at run time.
    if(TARGET stl2step_corpus_gen)
        list(APPEND _gates_full_depends corpus_generate)
    endif()
    set_tests_properties(gates_full PROPERTIES
        LABELS "gates"
        # G0.1 twin-builds the pinned baseline inside this test; an MSVC
        # cold build alone can exceed 600s, so the budget covers build+gates.
        TIMEOUT 3600
        REQUIRED_FILES "${_gates_full_required}")
    if(_gates_full_depends)
        set_tests_properties(gates_full PROPERTIES DEPENDS "${_gates_full_depends}")
    endif()
    unset(_gates_full_cmd)
    unset(_gates_full_depends)
    unset(_gates_full_required)
endif()

if(TARGET stl2step AND EXISTS "${_stl2step_gates_runner}" AND EXISTS "${_stl2step_gates_cube}")
    set(_gates_smoke_cmd
        "${Python3_EXECUTABLE}" "${_stl2step_gates_runner}"
        --smoke
        --binary "$<TARGET_FILE:stl2step>"
        --no-verify)
    set(_gates_smoke_required "$<TARGET_FILE:stl2step>")
    if(TARGET stl2step_census)
        list(APPEND _gates_smoke_cmd --census "$<TARGET_FILE:stl2step_census>")
        set(_gates_smoke_required "${_gates_smoke_required};$<TARGET_FILE:stl2step_census>")
    endif()
    add_test(NAME gates_smoke COMMAND ${_gates_smoke_cmd})
    set_tests_properties(gates_smoke PROPERTIES
        LABELS "gates;smoke"
        TIMEOUT 600
        REQUIRED_FILES "${_gates_smoke_required}")
    unset(_gates_smoke_cmd)
    unset(_gates_smoke_required)
endif()

# P2 standing gate: stl2step_p2buildtest over the canonical hand-fixture tree.
# Directory is p2buildtest (not buildtest) to dodge the historical unanchored
# build*/ ignore rule. Driver CLI: <regionset.json> <stl> [<expected.json>]
if(TARGET stl2step_p2buildtest)
    set(_p2_fixtures_root "${CMAKE_CURRENT_SOURCE_DIR}/build_fixtures")
    set(_p2_canonical_fixtures
        counterbore
        fillet_strip_2tri
        full_360_hole
        plate_half_hole
        quarter_round_1tri
        r1_round2_cascade
        r1_success
        r2_forced
        r8_360_explode_caps
        s09_mixed
        slotted_stadium)
    foreach(_p2_fx ${_p2_canonical_fixtures})
        set(_p2_rs  "${_p2_fixtures_root}/${_p2_fx}/${_p2_fx}.regionset.json")
        set(_p2_stl "${_p2_fixtures_root}/${_p2_fx}/${_p2_fx}.stl")
        set(_p2_exp "${_p2_fixtures_root}/${_p2_fx}/${_p2_fx}.expected.json")
        if(EXISTS "${_p2_rs}" AND EXISTS "${_p2_stl}")
            if(EXISTS "${_p2_exp}")
                add_test(NAME p2buildtest_${_p2_fx}
                    COMMAND stl2step_p2buildtest "${_p2_rs}" "${_p2_stl}" "${_p2_exp}")
            else()
                add_test(NAME p2buildtest_${_p2_fx}
                    COMMAND stl2step_p2buildtest "${_p2_rs}" "${_p2_stl}")
            endif()
            set_tests_properties(p2buildtest_${_p2_fx} PROPERTIES
                LABELS "gates;p2buildtest"
                TIMEOUT 120)
        endif()
        unset(_p2_rs)
        unset(_p2_stl)
        unset(_p2_exp)
    endforeach()
    unset(_p2_fx)
    unset(_p2_canonical_fixtures)
    unset(_p2_fixtures_root)
endif()

unset(_stl2step_gates_runner)
unset(_stl2step_gates_cube)
unset(_stl2step_gates_corpus)

# --- gates-ratchet ---
# Census-growth ratchet + waiver_audit (lane gates-ratchet). Additive only.
if(TARGET stl2step AND TARGET stl2step_census)
    set(_gr_hl_script "${CMAKE_CURRENT_SOURCE_DIR}/hl_census_ratchet.py")
    set(_gr_hl_stl "${CMAKE_CURRENT_SOURCE_DIR}/../corpus/handle-lock.stl")
    set(_gr_hl_ratchet "${CMAKE_CURRENT_SOURCE_DIR}/baseline/hl-ratchet.json")
    if(EXISTS "${_gr_hl_script}" AND EXISTS "${_gr_hl_stl}" AND EXISTS "${_gr_hl_ratchet}")
        add_test(NAME hl_census_ratchet
            COMMAND ${Python3_EXECUTABLE} "${_gr_hl_script}"
                --binary "$<TARGET_FILE:stl2step>"
                --census "$<TARGET_FILE:stl2step_census>"
                --stl "${_gr_hl_stl}"
                --ratchet "${_gr_hl_ratchet}")
        set_tests_properties(hl_census_ratchet PROPERTIES
            LABELS "gates"
            TIMEOUT 180
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;$<TARGET_FILE:stl2step_census>;${_gr_hl_stl};${_gr_hl_ratchet}")
        add_test(NAME hl_census_ratchet_selftest
            COMMAND ${Python3_EXECUTABLE} "${_gr_hl_script}" --self-test)
        set_tests_properties(hl_census_ratchet_selftest PROPERTIES
            LABELS "gates"
            TIMEOUT 30)
    endif()
    unset(_gr_hl_script)
    unset(_gr_hl_stl)
    unset(_gr_hl_ratchet)
endif()
if(TARGET stl2step AND TARGET stl2step_census AND TARGET stl2step_regiondump)
    set(_gr_waiver "${CMAKE_CURRENT_SOURCE_DIR}/../diag/handle-lock/waiver_audit.py")
    set(_gr_hl_script "${CMAKE_CURRENT_SOURCE_DIR}/hl_census_ratchet.py")
    if(EXISTS "${_gr_waiver}" AND EXISTS "${_gr_hl_script}")
        add_test(NAME waiver_audit
            COMMAND ${Python3_EXECUTABLE} "${_gr_hl_script}"
                --waiver-audit
                --waiver "${_gr_waiver}"
                --binary "$<TARGET_FILE:stl2step>"
                --census "$<TARGET_FILE:stl2step_census>"
                --dump "$<TARGET_FILE:stl2step_regiondump>")
        set_tests_properties(waiver_audit PROPERTIES
            LABELS "gates"
            TIMEOUT 180
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;$<TARGET_FILE:stl2step_census>;$<TARGET_FILE:stl2step_regiondump>;${_gr_waiver}")
    endif()
    unset(_gr_waiver)
    unset(_gr_hl_script)
endif()
