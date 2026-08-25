# AC5: stl2step_census /nonexistent.step exits 1 with JSON {"error":...} on stdout.
if(NOT CENSUS)
    message(FATAL_ERROR "CENSUS must be set")
endif()

execute_process(COMMAND "${CENSUS}" "/nonexistent.step"
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 1)
    message(FATAL_ERROR "expected exit 1, got ${rc}; stdout=${out}; stderr=${err}")
endif()
if(NOT out MATCHES "\"error\"")
    message(FATAL_ERROR "expected JSON error on stdout, got: ${out}")
endif()
