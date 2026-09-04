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

foreach(_stl2step_gate_dir census regiondump buildtest p2buildtest hybrid_validity)
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

# --- gates-j6watch ---
# J6 freeEdges watch on private Body9/12/18/20. Always registered so CI
# sees the label; SKIP (77) when STL2STEP_PRIVATE_CORPUS is unset or absent.
# Do not bake a private path into this file. API self-test needs no meshes.
set(_j6watch_script "${CMAKE_CURRENT_SOURCE_DIR}/j6_freeedges_watch.py")
if(EXISTS "${_j6watch_script}")
    add_test(NAME j6_freeedges_watch_api
        COMMAND "${Python3_EXECUTABLE}" "${_j6watch_script}" --self-test)
    set_tests_properties(j6_freeedges_watch_api PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step)
        add_test(NAME j6_freeedges_watch
            COMMAND "${Python3_EXECUTABLE}" "${_j6watch_script}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --jobs 4)
        set_tests_properties(j6_freeedges_watch PROPERTIES
            LABELS "gates;private-corpus"
            TIMEOUT 600
            SKIP_RETURN_CODE 77)
    endif()
endif()
unset(_j6watch_script)

# --- cascade-r4 ---
# D1 §5.1 resurrection ban. Greps src/; no engine required.
set(_cascade_nowaiver "${CMAKE_CURRENT_SOURCE_DIR}/cascade_no_waiver.py")
if(EXISTS "${_cascade_nowaiver}")
    add_test(NAME cascade_no_waiver
        COMMAND "${Python3_EXECUTABLE}" "${_cascade_nowaiver}")
    set_tests_properties(cascade_no_waiver PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    add_test(NAME cascade_no_waiver_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_cascade_nowaiver}" --self-test)
    set_tests_properties(cascade_no_waiver_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
endif()
unset(_cascade_nowaiver)

# --- gates-sweep ---
# Standalone 11-file TrueForm sweep (scripts/stress_sweep.py). Not part of
# gates_full and not sharded. SKIPs (exit 77) when the private STL corpus
# is absent so CI stays green.
if(TARGET stl2step)
    get_filename_component(_stl2step_sweep_root
        "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    set(_stl2step_sweep_py "${_stl2step_sweep_root}/scripts/stress_sweep.py")
    if(EXISTS "${_stl2step_sweep_py}")
        add_test(NAME stress_sweep
            COMMAND "${Python3_EXECUTABLE}" "${_stl2step_sweep_py}"
                --binary "$<TARGET_FILE:stl2step>"
                --jobs 4
                --json "${CMAKE_CURRENT_BINARY_DIR}/stress_sweep.json"
                --markdown "${CMAKE_CURRENT_BINARY_DIR}/stress_sweep.md")
        set_tests_properties(stress_sweep PROPERTIES
            LABELS "stress-sweep;private-corpus"
            TIMEOUT 600
            SKIP_RETURN_CODE 77
            REQUIRED_FILES "$<TARGET_FILE:stl2step>")
    endif()
    unset(_stl2step_sweep_py)
    unset(_stl2step_sweep_root)
endif()

# --- gates-partial ---
# DECISION-boundary pickup floors + measured PRG ratchet (D-S3-62).
# STRICT=1: floors (18/30/48, rev=0) are hard. The three PRG cells
# (tightBudget, census.radii, d3f7 R=3) are never-get-worse ceilings
# in tests/gates/baseline/prg-ratchet.json — FAIL if any cell gets
# worse; print current values. Not a pass of the R4 formula budget.
#
# FLIPPED 2026-09-02 at 2971d31: HP ships analytic; EXPECTED-RED
# invert retired. PASS_REGULAR_EXPRESSION removed. Keep
# ENVIRONMENT PARTIAL_RECOVERY_STRICT=1. Do not lower FLOOR_*.
#
# D3F-4/6: floors 18/30/48 (charter combined>=100 STRUCK). D3F-7
# R=3 is ratcheted; other-class over-segmentation stays hard.
set(_prg_script "${CMAKE_CURRENT_SOURCE_DIR}/partial_recovery_gate.py")
set(_prg_ratchet "${CMAKE_CURRENT_SOURCE_DIR}/baseline/prg-ratchet.json")
if(EXISTS "${_prg_script}")
    add_test(NAME partial_recovery_gate_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_prg_script}" --self-test)
    set_tests_properties(partial_recovery_gate_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    add_test(NAME partial_recovery_gate_synthetic
        COMMAND "${Python3_EXECUTABLE}" "${_prg_script}" --synthetic-pass)
    set_tests_properties(partial_recovery_gate_synthetic PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step AND EXISTS "${_prg_ratchet}")
        add_test(NAME partial_recovery_gate
            COMMAND "${Python3_EXECUTABLE}" "${_prg_script}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --ratchet "${_prg_ratchet}"
                    --jobs 4)
        set_tests_properties(partial_recovery_gate PROPERTIES
            LABELS "gates"
            TIMEOUT 300
            ENVIRONMENT "PARTIAL_RECOVERY_STRICT=1"
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_prg_ratchet}")
    endif()
endif()
unset(_prg_script)
unset(_prg_ratchet)

# --- gates-g01 (wave 6 / p3) ---
# Four D3-2 / D3-8 / D3-9 instruments. Append-only. Do not edit run_gates.py.
# seam_stored_both_sides + fillet_exclusion_e0 follow the expected-red protocol
# of partial_recovery_gate (STRICT default ON, PASS_REGULAR_EXPRESSION invert,
# UNEXPECTED stays loud). collapse_fail_ratchet and hybrid_validity_pre_j4 are
# green at dc49e45.

set(_g01_seam "${CMAKE_CURRENT_SOURCE_DIR}/seam_stored_both_sides.py")
if(EXISTS "${_g01_seam}")
    add_test(NAME seam_stored_both_sides_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_g01_seam}" --self-test)
    set_tests_properties(seam_stored_both_sides_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    add_test(NAME seam_stored_both_sides_synthetic
        COMMAND "${Python3_EXECUTABLE}" "${_g01_seam}" --synthetic-pass)
    set_tests_properties(seam_stored_both_sides_synthetic PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step)
        add_test(NAME seam_stored_both_sides
            COMMAND "${Python3_EXECUTABLE}" "${_g01_seam}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --jobs 2)
        set_tests_properties(seam_stored_both_sides PROPERTIES
            LABELS "gates;expected-red"
            TIMEOUT 180
            ENVIRONMENT "SEAM_STORED_BOTH_SIDES_STRICT=1"
            PASS_REGULAR_EXPRESSION "SEAM_STORED_BOTH_SIDES_EXPECTED_RED"
            REQUIRED_FILES "$<TARGET_FILE:stl2step>")
    endif()
endif()
unset(_g01_seam)

set(_g01_coll "${CMAKE_CURRENT_SOURCE_DIR}/collapse_fail_ratchet.py")
set(_g01_coll_ratchet "${CMAKE_CURRENT_SOURCE_DIR}/baseline/collapse-ratchet.json")
if(EXISTS "${_g01_coll}")
    add_test(NAME collapse_fail_ratchet_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_g01_coll}" --self-test)
    set_tests_properties(collapse_fail_ratchet_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    add_test(NAME collapse_fail_ratchet_synthetic
        COMMAND "${Python3_EXECUTABLE}" "${_g01_coll}" --synthetic-pass)
    set_tests_properties(collapse_fail_ratchet_synthetic PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step AND EXISTS "${_g01_coll_ratchet}")
        add_test(NAME collapse_fail_ratchet
            COMMAND "${Python3_EXECUTABLE}" "${_g01_coll}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --ratchet "${_g01_coll_ratchet}"
                    --jobs 2)
        set_tests_properties(collapse_fail_ratchet PROPERTIES
            LABELS "gates"
            TIMEOUT 180
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_g01_coll_ratchet}")
    endif()
endif()
unset(_g01_coll)
unset(_g01_coll_ratchet)

set(_g01_e0 "${CMAKE_CURRENT_SOURCE_DIR}/fillet_exclusion_e0.py")
set(_g01_e0_ratchet "${CMAKE_CURRENT_SOURCE_DIR}/baseline/fillet-e0-ratchet.json")
if(EXISTS "${_g01_e0}")
    add_test(NAME fillet_exclusion_e0_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_g01_e0}" --self-test)
    set_tests_properties(fillet_exclusion_e0_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    add_test(NAME fillet_exclusion_e0_synthetic
        COMMAND "${Python3_EXECUTABLE}" "${_g01_e0}" --synthetic-pass)
    set_tests_properties(fillet_exclusion_e0_synthetic PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step AND TARGET stl2step_regiondump AND EXISTS "${_g01_e0_ratchet}")
        add_test(NAME fillet_exclusion_e0
            COMMAND "${Python3_EXECUTABLE}" "${_g01_e0}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --dump "$<TARGET_FILE:stl2step_regiondump>"
                    --ratchet "${_g01_e0_ratchet}"
                    --jobs 4)
        set_tests_properties(fillet_exclusion_e0 PROPERTIES
            LABELS "gates;expected-red"
            TIMEOUT 180
            ENVIRONMENT "FILLET_EXCLUSION_E0_STRICT=1"
            PASS_REGULAR_EXPRESSION "FILLET_EXCLUSION_E0_EXPECTED_RED"
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;$<TARGET_FILE:stl2step_regiondump>;${_g01_e0_ratchet}"
            DEPENDS corpus_generate)
    endif()
endif()
unset(_g01_e0)
unset(_g01_e0_ratchet)

# --- gates-g01 census allowance ---
# Anti-silent-waiver: censusValidExpected:false set is exactly {Body11}.
set(_g01_allow "${CMAKE_CURRENT_SOURCE_DIR}/corpus_census_allowance.py")
if(EXISTS "${_g01_allow}")
    add_test(NAME corpus_census_allowance
        COMMAND "${Python3_EXECUTABLE}" "${_g01_allow}"
                --corpus "${CMAKE_CURRENT_SOURCE_DIR}/../corpus")
    set_tests_properties(corpus_census_allowance PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
endif()
unset(_g01_allow)

# --- gates-g01 plane census (D-FL-6) ---
# Per-GT-plane face-count ratchet. Never expected-red.
set(_g01_plane "${CMAKE_CURRENT_SOURCE_DIR}/plane_census_gate.py")
set(_g01_plane_ratchet "${CMAKE_CURRENT_SOURCE_DIR}/baseline/plane-census-ratchet.json")
if(EXISTS "${_g01_plane}")
    add_test(NAME plane_census_gate_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_g01_plane}" --self-test)
    set_tests_properties(plane_census_gate_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step AND EXISTS "${_g01_plane_ratchet}")
        add_test(NAME plane_census_gate
            COMMAND "${Python3_EXECUTABLE}" "${_g01_plane}"
                    --binary "$<TARGET_FILE:stl2step>")
        set_tests_properties(plane_census_gate PROPERTIES
            LABELS "gates"
            TIMEOUT 180
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_g01_plane_ratchet}")
    endif()
endif()
unset(_g01_plane)
unset(_g01_plane_ratchet)
# --- 130-BIND edge-class red lines (D-130-2) ---
# unhandled / overTol / overCap == 0 over the corpus + the P2 build fixtures,
# tier-2 mesh-polyline count ratcheted per fixture. Never expected-red as a
# CTest property: a fixture that violates a line on the branch that seeded the
# ratchet is named in edge-class-ratchet.json's expectedRed and stays loud.
set(_edge_class "${CMAKE_CURRENT_SOURCE_DIR}/edge_class_gate.py")
set(_edge_class_ratchet "${CMAKE_CURRENT_SOURCE_DIR}/baseline/edge-class-ratchet.json")
if(EXISTS "${_edge_class}")
    add_test(NAME edge_class_gate_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_edge_class}" --self-test)
    set_tests_properties(edge_class_gate_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 60)
    add_test(NAME edge_class_gate_synthetic
        COMMAND "${Python3_EXECUTABLE}" "${_edge_class}" --synthetic-pass)
    set_tests_properties(edge_class_gate_synthetic PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step AND EXISTS "${_edge_class_ratchet}")
        add_test(NAME edge_class_gate
            COMMAND "${Python3_EXECUTABLE}" "${_edge_class}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --ratchet "${_edge_class_ratchet}"
                    --jobs 4)
        set_tests_properties(edge_class_gate PROPERTIES
            LABELS "gates"
            TIMEOUT 900
            REQUIRED_FILES "$<TARGET_FILE:stl2step>;${_edge_class_ratchet}")
        if(TARGET stl2step_corpus_gen)
            set_tests_properties(edge_class_gate PROPERTIES DEPENDS corpus_generate)
        endif()
    endif()
endif()
unset(_edge_class)
unset(_edge_class_ratchet)

# --- NW-CYL never-worse (D-S3-15, gates-g01-delta-nwcyl) ---
set(_nwcyl_script "${CMAKE_CURRENT_SOURCE_DIR}/nwcyl_gate.py")
set(_nwcyl_baseline "${CMAKE_CURRENT_SOURCE_DIR}/nwcyl_baseline.json")
if(EXISTS "${_nwcyl_script}")
    add_test(NAME nwcyl_gate_selftest
        COMMAND "${Python3_EXECUTABLE}" "${_nwcyl_script}" --self-test)
    set_tests_properties(nwcyl_gate_selftest PROPERTIES
        LABELS "gates"
        TIMEOUT 30)
    if(TARGET stl2step AND TARGET stl2step_census AND TARGET stl2step_regiondump
            AND EXISTS "${_nwcyl_baseline}")
        add_test(NAME nwcyl_gate
            COMMAND "${Python3_EXECUTABLE}" "${_nwcyl_script}"
                    --binary "$<TARGET_FILE:stl2step>"
                    --census "$<TARGET_FILE:stl2step_census>"
                    --dump "$<TARGET_FILE:stl2step_regiondump>"
                    --baseline "${_nwcyl_baseline}"
                    --jobs 4)
        set_tests_properties(nwcyl_gate PROPERTIES
            LABELS "gates"
            TIMEOUT 600
            REQUIRED_FILES
                "$<TARGET_FILE:stl2step>;$<TARGET_FILE:stl2step_census>;$<TARGET_FILE:stl2step_regiondump>;${_nwcyl_baseline}")
        if(TARGET stl2step_corpus_gen)
            set_tests_properties(nwcyl_gate PROPERTIES DEPENDS corpus_generate)
        endif()
    endif()
endif()
unset(_nwcyl_script)
unset(_nwcyl_baseline)

# --- 130-FIXTURES battery gate (expected-red until Wave 2 detectors land) ---
# STRICT default ON, PASS_REGULAR_EXPRESSION invert, UNEXPECTED stays loud.
# Flip: remove PASS_REGULAR_EXPRESSION; keep GATE_130_STRICT=1.
set(_g130 "${CMAKE_CURRENT_SOURCE_DIR}/gate_130.py")
if(EXISTS "${_g130}")
    if(TARGET stl2step)
        set(_g130_cmd
            "${Python3_EXECUTABLE}" "${_g130}"
            --binary "$<TARGET_FILE:stl2step>"
            --corpus "${CMAKE_CURRENT_SOURCE_DIR}/../corpus")
        set(_g130_req "$<TARGET_FILE:stl2step>")
        if(TARGET stl2step_census)
            list(APPEND _g130_cmd --census "$<TARGET_FILE:stl2step_census>")
            set(_g130_req "${_g130_req};$<TARGET_FILE:stl2step_census>")
        endif()
        add_test(NAME gate_130 COMMAND ${_g130_cmd})
        set_tests_properties(gate_130 PROPERTIES
            LABELS "gates;expected-red"
            TIMEOUT 600
            ENVIRONMENT "GATE_130_STRICT=1"
            PASS_REGULAR_EXPRESSION "GATE_130_EXPECTED_RED"
            REQUIRED_FILES "${_g130_req}")
        if(TARGET stl2step_corpus_gen)
            set_tests_properties(gate_130 PROPERTIES DEPENDS corpus_generate)
        endif()
        unset(_g130_cmd)
        unset(_g130_req)
    endif()
endif()
unset(_g130)
