cmake_minimum_required(VERSION 3.25)

project(TriggerFishPercussionWorkbench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

get_filename_component(TRIGGERFISH_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

set(TRIGGERFISH_CRASH_SOURCES
  ${TRIGGERFISH_ROOT}/src/tfdsp/sampleRate.cpp
  ${TRIGGERFISH_ROOT}/src/tfdsp/percussion/crash_cymbal.cpp
  ${TRIGGERFISH_ROOT}/src/tfdsp/percussion/crash_cymbal_parameters.cpp
  ${CMAKE_CURRENT_LIST_DIR}/wasm/crash_api.cpp
  ${CMAKE_CURRENT_LIST_DIR}/wasm/crash_macros.cpp
)

function(triggerfish_workbench_target target)
  target_sources(${target} PRIVATE ${TRIGGERFISH_CRASH_SOURCES})
  target_include_directories(${target} PRIVATE
    ${TRIGGERFISH_ROOT}/src
    ${TRIGGERFISH_ROOT}/vendor/eigen
    ${CMAKE_CURRENT_LIST_DIR}/wasm
  )
  if(MINGW)
    target_link_options(${target} PRIVATE -static -static-libgcc -static-libstdc++)
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libwinpthread.a
      OUTPUT_VARIABLE MINGW_STATIC_WINPTHREAD
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT EXISTS "${MINGW_STATIC_WINPTHREAD}")
      message(FATAL_ERROR "MinGW static winpthread archive was not found")
    endif()
    target_link_libraries(${target} PRIVATE "${MINGW_STATIC_WINPTHREAD}")
  endif()
endfunction()

if(EMSCRIPTEN)
  add_executable(triggerfish_percussion_wasm)
  triggerfish_workbench_target(triggerfish_percussion_wasm)
  set_target_properties(triggerfish_percussion_wasm PROPERTIES
    OUTPUT_NAME "triggerfish-percussion"
    SUFFIX ".mjs"
  )
  target_compile_options(triggerfish_percussion_wasm PRIVATE -msimd128)
  target_link_options(triggerfish_percussion_wasm PRIVATE
    "-msimd128"
    "--no-entry"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sENVIRONMENT=web,worker,node"
    "-sFILESYSTEM=0"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sSTACK_SIZE=2097152"
    "-sASSERTIONS=1"
    "-sEXPORTED_FUNCTIONS=['_tf_crash_api_version','_tf_crash_create','_tf_crash_destroy','_tf_crash_reset','_tf_crash_trigger','_tf_crash_set_mute','_tf_crash_process','_tf_crash_macro_count','_tf_crash_macro_key','_tf_crash_macro_name','_tf_crash_macro_unit','_tf_crash_macro_scale','_tf_crash_macro_minimum','_tf_crash_macro_maximum','_tf_crash_macro_default','_tf_crash_macro_get','_tf_crash_macro_set','_tf_crash_macro_commit','_malloc','_free']"
    "-sEXPORTED_RUNTIME_METHODS=['HEAPF32','UTF8ToString']"
  )

  # A self-contained copy of the same C++ API for the browser audio rendering
  # thread. AudioWorkletGlobalScope has no fetch API, so the Wasm bytes must be
  # embedded in the module loaded by audioWorklet.addModule().
  add_executable(triggerfish_percussion_audio_wasm)
  triggerfish_workbench_target(triggerfish_percussion_audio_wasm)
  set_target_properties(triggerfish_percussion_audio_wasm PROPERTIES
    OUTPUT_NAME "triggerfish-percussion-audio"
    SUFFIX ".mjs"
  )
  target_compile_options(triggerfish_percussion_audio_wasm PRIVATE -msimd128)
  target_link_options(triggerfish_percussion_audio_wasm PRIVATE
    "-msimd128"
    "--no-entry"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sENVIRONMENT=worklet"
    "-sSINGLE_FILE=1"
    "-sWASM_ASYNC_COMPILATION=0"
    "-sFILESYSTEM=0"
    "-sALLOW_MEMORY_GROWTH=0"
    "-sINITIAL_MEMORY=33554432"
    "-sSTACK_SIZE=2097152"
    "-sASSERTIONS=1"
    "-sEXPORTED_FUNCTIONS=['_tf_crash_create','_tf_crash_destroy','_tf_crash_reset','_tf_crash_trigger','_tf_crash_process','_tf_crash_macro_count','_tf_crash_macro_default','_tf_crash_macro_set','_tf_crash_macro_commit','_malloc','_free']"
    "-sEXPORTED_RUNTIME_METHODS=['HEAPF32']"
  )

  enable_testing()
  find_program(NODE_EXECUTABLE NAMES node nodejs)
  if(NODE_EXECUTABLE)
    add_test(NAME triggerfish_percussion_wasm_smoke
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/tests/wasm_smoke.mjs
        $<TARGET_FILE:triggerfish_percussion_wasm>
    )
    add_test(NAME triggerfish_workbench_analysis_tests
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/tests/analysis_tests.mjs
    )
    add_test(NAME triggerfish_workbench_limiter_tests
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/tests/limiter_tests.mjs
    )
    add_test(NAME triggerfish_workbench_size_meta_tests
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/tests/size_meta_tests.mjs
    )
    add_test(NAME triggerfish_workbench_state_tests
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/tests/state_tests.mjs
    )
    add_test(NAME triggerfish_workbench_settings_tests
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/tests/settings_tests.mjs
    )
    add_custom_target(triggerfish_workbench_site
      COMMAND ${NODE_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/web/build.mjs
        --module $<TARGET_FILE:triggerfish_percussion_wasm>
        --audio-module $<TARGET_FILE:triggerfish_percussion_audio_wasm>
        --output ${CMAKE_BINARY_DIR}/site
      DEPENDS triggerfish_percussion_wasm triggerfish_percussion_audio_wasm
      COMMENT "Assembling the optional percussion workbench site"
    )
  else()
    message(STATUS "Node was not found; the Wasm smoke test is unavailable")
  endif()
else()
  add_executable(triggerfish_workbench_api_tests
    ${CMAKE_CURRENT_LIST_DIR}/tests/crash_api_tests.cpp
  )
  triggerfish_workbench_target(triggerfish_workbench_api_tests)
  add_executable(triggerfish_workbench_signature
    ${CMAKE_CURRENT_LIST_DIR}/tests/crash_signature.cpp
  )
  triggerfish_workbench_target(triggerfish_workbench_signature)
  enable_testing()
  add_test(NAME triggerfish_workbench_api_tests
    COMMAND triggerfish_workbench_api_tests)
endif()
