# Byte-determinism: run corpus generator twice into two directories; diff -r must be empty.
if(NOT GEN OR NOT DIR_A OR NOT DIR_B)
    message(FATAL_ERROR "GEN, DIR_A, and DIR_B must be set")
endif()

file(MAKE_DIRECTORY "${DIR_A}")
file(MAKE_DIRECTORY "${DIR_B}")

execute_process(COMMAND "${GEN}" "${DIR_A}" RESULT_VARIABLE rc1)
if(NOT rc1 EQUAL 0)
    message(FATAL_ERROR "first generation failed (rc=${rc1})")
endif()
execute_process(COMMAND "${GEN}" "${DIR_B}" RESULT_VARIABLE rc2)
if(NOT rc2 EQUAL 0)
    message(FATAL_ERROR "second generation failed (rc=${rc2})")
endif()

find_program(_DIFF diff)
if(NOT _DIFF)
    message(FATAL_ERROR "diff not found")
endif()

execute_process(
    COMMAND "${_DIFF}" -qr "${DIR_A}" "${DIR_B}"
    RESULT_VARIABLE diff_rc
    OUTPUT_VARIABLE diff_out
    ERROR_VARIABLE diff_err)
if(NOT diff_rc EQUAL 0)
    message(FATAL_ERROR "corpus not byte-deterministic:\n${diff_out}${diff_err}")
endif()
message(STATUS "corpus determinism OK")
