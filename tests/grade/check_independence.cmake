# Case 10: grader sources stay independent of the engine, harness, and <regex>.
# Forbidden substrings are searched with string(FIND) so this file itself never
# needs a matching-engine invocation.

get_filename_component(_grade_dir "${CMAKE_CURRENT_LIST_DIR}" ABSOLUTE)

set(_forbidden
    "MeshView"
    "stl2step::core"
    "stl2step_harness"
    "stl2step_core"
    "refit"
    "regex")

file(GLOB _grade_files "${_grade_dir}/*")
set(_fail "")
foreach(_f IN LISTS _grade_files)
    if(IS_DIRECTORY "${_f}")
        continue()
    endif()
    file(READ "${_f}" _txt)
    get_filename_component(_name "${_f}" NAME)

    foreach(_bad IN LISTS _forbidden)
        string(FIND "${_txt}" "${_bad}" _pos)
        if(NOT _pos EQUAL -1)
            # This checker names the banned tokens; skip self.
            if(_name STREQUAL "check_independence.cmake")
                continue()
            endif()
            string(APPEND _fail "  ${_name}: contains '${_bad}'\n")
        endif()
    endforeach()

    # Includes from src/ other than stl_quant.hpp.
    string(FIND "${_txt}" "#include" _inc)
    if(NOT _inc EQUAL -1)
        string(FIND "${_txt}" "src/" _srcpos)
        if(NOT _srcpos EQUAL -1)
            string(FIND "${_txt}" "stl_quant.hpp" _ok)
            if(_ok EQUAL -1)
                string(APPEND _fail "  ${_name}: src/ include that is not stl_quant.hpp\n")
            endif()
        endif()
    endif()
endforeach()

if(_fail)
    message(FATAL_ERROR "grade_independence failed:\n${_fail}")
endif()
message(STATUS "grade_independence: clean")
