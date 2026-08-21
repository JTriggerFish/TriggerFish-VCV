# Repository-level Windows toolchain for scikit-build-core/uv.
#
# MSYS2_ROOT may relocate MSYS2; it defaults to the same location used by
# dev.ps1.  Selecting the compilers in a toolchain file makes the choice before
# CMake performs compiler detection, so an installed MSVC toolchain can never
# be selected for the Python extension.
if(NOT CMAKE_HOST_WIN32)
  return()
endif()

if(DEFINED ENV{MSYS2_ROOT} AND NOT "$ENV{MSYS2_ROOT}" STREQUAL "")
  file(TO_CMAKE_PATH "$ENV{MSYS2_ROOT}" TRIGGERFISH_MSYS2_ROOT)
else()
  set(TRIGGERFISH_MSYS2_ROOT "C:/msys64")
endif()

set(TRIGGERFISH_MINGW_BIN "${TRIGGERFISH_MSYS2_ROOT}/mingw64/bin")
set(TRIGGERFISH_GCC "${TRIGGERFISH_MINGW_BIN}/gcc.exe")
set(TRIGGERFISH_GXX "${TRIGGERFISH_MINGW_BIN}/g++.exe")

if(NOT EXISTS "${TRIGGERFISH_GCC}" OR NOT EXISTS "${TRIGGERFISH_GXX}")
  message(FATAL_ERROR
    "TriggerFish requires MSYS2 MinGW64 for Windows Python builds. "
    "Expected gcc.exe and g++.exe under '${TRIGGERFISH_MINGW_BIN}'. "
    "Install MSYS2 there or set MSYS2_ROOT before invoking uv.")
endif()

# GCC's cc1/assembler/linker subprocesses load MinGW runtime DLLs from this
# directory.  An absolute g++.exe alone is insufficient when uv was launched
# from an ordinary PowerShell whose PATH does not already contain MSYS2.
set(ENV{PATH} "${TRIGGERFISH_MINGW_BIN};$ENV{PATH}")

# scikit-build-core launches the build as a separate subprocess after CMake
# configuration, so the configure-time environment edit above does not reach
# Ninja.  Prefix every compiler invocation with `cmake -E env` as well.
string(REPLACE ";" "\\;" TRIGGERFISH_ESCAPED_BUILD_PATH "$ENV{PATH}")
set(CMAKE_C_COMPILER_LAUNCHER
    "${CMAKE_COMMAND};-E;env;PATH=${TRIGGERFISH_ESCAPED_BUILD_PATH}"
    CACHE STRING "MinGW runtime environment for C compilations" FORCE)
set(CMAKE_CXX_COMPILER_LAUNCHER
    "${CMAKE_COMMAND};-E;env;PATH=${TRIGGERFISH_ESCAPED_BUILD_PATH}"
    CACHE STRING "MinGW runtime environment for C++ compilations" FORCE)
set(CMAKE_C_LINKER_LAUNCHER
    "${CMAKE_COMMAND};-E;env;PATH=${TRIGGERFISH_ESCAPED_BUILD_PATH}"
    CACHE STRING "MinGW runtime environment for C links" FORCE)
set(CMAKE_CXX_LINKER_LAUNCHER
    "${CMAKE_COMMAND};-E;env;PATH=${TRIGGERFISH_ESCAPED_BUILD_PATH}"
    CACHE STRING "MinGW runtime environment for C++ links" FORCE)

set(CMAKE_C_COMPILER "${TRIGGERFISH_GCC}" CACHE FILEPATH
    "TriggerFish MinGW C compiler" FORCE)
set(CMAKE_CXX_COMPILER "${TRIGGERFISH_GXX}" CACHE FILEPATH
    "TriggerFish MinGW C++ compiler" FORCE)
