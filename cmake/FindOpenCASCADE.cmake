# FindOpenCASCADE.cmake - locate OpenCASCADE (OCCT) portably across platforms.
#
# Strategy:
#   1. Prefer OCCT's own exported package (OpenCASCADEConfig.cmake). Homebrew,
#      Debian/Ubuntu (libocct-*-dev), vcpkg and the official installers all ship
#      it, and it knows the exact toolkit layout for its own version.
#   2. Fall back to manual header + per-toolkit library discovery, handling the
#      7.8 rename of the STEP/STL toolkits (TKSTEP*/TKSTL -> TKDESTEP/TKDESTL).
#
# Results:
#   OpenCASCADE_FOUND          - TRUE when everything the engine needs was located
#   OpenCASCADE_INCLUDE_DIRS   - header directory (the one containing Standard.hxx)
#   OpenCASCADE_LIBRARIES      - the toolkits to link
#   OpenCASCADE::OpenCASCADE   - imported INTERFACE target carrying both of the above
#
# Hints (any of): OpenCASCADE_ROOT, CASROOT (env or cache), CMAKE_PREFIX_PATH.

# The toolkits the engine actually references, most-dependent first. Old names
# are listed as fallbacks for the manual path (OCCT < 7.8).
set(_occt_required
    TKDESTL        # STL read/write (was TKSTL)
    TKDESTEP       # STEP read/write (was TKSTEP)
    TKXSBase       # data-exchange base
    TKShHealing    # ShapeFix / sewing repair
    TKMesh         # BRepMesh_IncrementalMesh (--mesh mode)
    TKTopAlgo TKGeomAlgo TKBRep TKG3d TKG2d TKGeomBase TKMath TKernel)

set(_occt_fallback_TKDESTL  TKSTL)
set(_occt_fallback_TKDESTEP TKSTEP)

# ---------------------------------------------------------------------------
# 1. Try the config package first (most reliable; version-aware).
# ---------------------------------------------------------------------------
if(NOT OpenCASCADE_FOUND)
    find_package(OpenCASCADE CONFIG QUIET)
endif()

set(_occt_include "")
if(DEFINED OpenCASCADE_INCLUDE_DIR)
    set(_occt_include "${OpenCASCADE_INCLUDE_DIR}")
elseif(DEFINED OpenCASCADE_INCLUDE_DIRS)
    set(_occt_include "${OpenCASCADE_INCLUDE_DIRS}")
endif()

set(_occt_libs "")
set(_occt_ok TRUE)

if(OpenCASCADE_FOUND AND _occt_include)
    # Config mode succeeded: OCCT exposes each toolkit as an imported target
    # whose bare name matches the toolkit (e.g. TKernel, TKDESTEP).
    foreach(_tk IN LISTS _occt_required)
        if(TARGET ${_tk})
            list(APPEND _occt_libs ${_tk})
        else()
            set(_fb "${_occt_fallback_${_tk}}")
            if(_fb AND TARGET ${_fb})
                list(APPEND _occt_libs ${_fb})
            else()
                # Toolkit target missing - resolve by name in the library dir.
                find_library(_occt_lib_${_tk} NAMES ${_tk} ${_fb}
                             HINTS "${OpenCASCADE_LIBRARY_DIR}" "${OpenCASCADE_BINARY_DIR}")
                if(_occt_lib_${_tk})
                    list(APPEND _occt_libs "${_occt_lib_${_tk}}")
                else()
                    set(_occt_ok FALSE)
                endif()
            endif()
        endif()
    endforeach()
else()
    # -----------------------------------------------------------------------
    # 2. Manual discovery.
    # -----------------------------------------------------------------------
    set(_occt_ok FALSE)  # re-established below only if headers + libs are found

    find_path(OpenCASCADE_INCLUDE_DIR
        NAMES Standard.hxx
        HINTS
            ${OpenCASCADE_ROOT} $ENV{OpenCASCADE_ROOT} $ENV{CASROOT} ${CASROOT}
        PATH_SUFFIXES
            include/opencascade opencascade include inc
        PATHS
            /usr/include/opencascade /usr/local/include/opencascade
            /opt/homebrew/include/opencascade /opt/local/include/opencascade)
    set(_occt_include "${OpenCASCADE_INCLUDE_DIR}")

    if(_occt_include)
        set(_occt_ok TRUE)
        foreach(_tk IN LISTS _occt_required)
            set(_fb "${_occt_fallback_${_tk}}")
            find_library(_occt_lib_${_tk}
                NAMES ${_tk} ${_fb}
                HINTS
                    ${OpenCASCADE_ROOT} $ENV{OpenCASCADE_ROOT} $ENV{CASROOT} ${CASROOT}
                PATH_SUFFIXES lib lib64 win64/vc14/lib win64/gcc/lib
                PATHS
                    /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
                    /opt/homebrew/lib /opt/local/lib)
            if(_occt_lib_${_tk})
                list(APPEND _occt_libs "${_occt_lib_${_tk}}")
            else()
                set(_occt_ok FALSE)
            endif()
        endforeach()
    endif()
endif()

if(_occt_ok AND _occt_include AND _occt_libs)
    set(OpenCASCADE_INCLUDE_DIRS "${_occt_include}")
    set(OpenCASCADE_LIBRARIES    "${_occt_libs}")
else()
    set(OpenCASCADE_INCLUDE_DIRS "")
    set(OpenCASCADE_LIBRARIES    "")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenCASCADE
    REQUIRED_VARS OpenCASCADE_INCLUDE_DIRS OpenCASCADE_LIBRARIES
    VERSION_VAR   OpenCASCADE_VERSION)

if(OpenCASCADE_FOUND AND NOT TARGET OpenCASCADE::OpenCASCADE)
    add_library(OpenCASCADE::OpenCASCADE INTERFACE IMPORTED)
    set_target_properties(OpenCASCADE::OpenCASCADE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${OpenCASCADE_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES      "${OpenCASCADE_LIBRARIES}")
endif()

unset(_occt_required)
unset(_occt_include)
unset(_occt_libs)
unset(_occt_ok)
