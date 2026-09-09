# D-I-5 / GitHub #8 — single source for the calibrated OCCT test battery.
# Tests calibrated on OCCT 7.9.x skip with a reason on older OCCT; hosted CI pins 7.9.
set(STL2STEP_OCCT_CALIBRATED_MIN_MAJOR 7)
set(STL2STEP_OCCT_CALIBRATED_MIN_MINOR 9)
set(STL2STEP_OCCT_SKIP_EXIT_CODE 77)

# G0.1 twin-run baseline (stl2step 1.0.0); must be resolvable via git cat-file -e.
set(STL2STEP_G01_BASELINE_COMMIT "187ead0d8cf3d3694153cbcff9314d65324fec63")
