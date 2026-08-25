# Run the census CLI twice on the same STEP; stdout must be byte-identical.
if(NOT CENSUS OR NOT FILE)
    message(FATAL_ERROR "CENSUS and FILE must be set")
endif()
if(NOT EXISTS "${FILE}")
    message(FATAL_ERROR "missing STEP fixture: ${FILE}")
endif()

execute_process(COMMAND "${CENSUS}" "${FILE}"
    OUTPUT_VARIABLE out1 RESULT_VARIABLE rc1)
if(NOT rc1 EQUAL 0)
    message(FATAL_ERROR "first census failed on ${FILE} (rc=${rc1})")
endif()
execute_process(COMMAND "${CENSUS}" "${FILE}"
    OUTPUT_VARIABLE out2 RESULT_VARIABLE rc2)
if(NOT rc2 EQUAL 0)
    message(FATAL_ERROR "second census failed on ${FILE} (rc=${rc2})")
endif()
if(NOT out1 STREQUAL out2)
    message(FATAL_ERROR "non-deterministic census for ${FILE}")
endif()
