# Run the dump tool twice on each fixture; stdout must be byte-identical.
if(NOT DUMP OR NOT FILE1)
    message(FATAL_ERROR "DUMP and FILE1 must be set")
endif()

function(dump_twice stl)
    execute_process(COMMAND "${DUMP}" "${stl}"
        OUTPUT_VARIABLE out1 RESULT_VARIABLE rc1)
    if(NOT rc1 EQUAL 0)
        message(FATAL_ERROR "first dump failed on ${stl} (rc=${rc1})")
    endif()
    execute_process(COMMAND "${DUMP}" "${stl}"
        OUTPUT_VARIABLE out2 RESULT_VARIABLE rc2)
    if(NOT rc2 EQUAL 0)
        message(FATAL_ERROR "second dump failed on ${stl} (rc=${rc2})")
    endif()
    if(NOT out1 STREQUAL out2)
        message(FATAL_ERROR "non-deterministic dump for ${stl}")
    endif()
endfunction()

dump_twice("${FILE1}")
if(FILE2)
    dump_twice("${FILE2}")
endif()
if(FILE3)
    dump_twice("${FILE3}")
endif()
