cmake_minimum_required(VERSION 3.25)

project(TriggerFishDSP VERSION 2.4.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(triggerfish_dsp STATIC
	src/tfseq_compiler.cpp
	src/tfseq.hpp
	src/tfseq_parser.cpp
	src/tfseq_parser.hpp
	src/tfseq_runtime.cpp
	src/tfdsp/sampleRate.cpp
	src/tfdsp/early_reflections.hpp
	src/tfdsp/early_reflections_worker.cpp
	src/tfdsp/early_reflections_worker.hpp
	src/tfdsp/late_reverb.hpp
	src/tfdsp/late_reverb_coefficients.hpp
	src/tfdsp/multiband_decay_filter.hpp
	src/tfdsp/velvet_feedback_matrix.hpp
	src/tfdsp/room_reverb.hpp
	src/tfdsp/scene_pack.hpp
	src/tfdsp/shimmer_velvet_diffuser.hpp
	src/tfdsp/windowed_pitch_shifter.hpp
)
set_target_properties(triggerfish_dsp PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(triggerfish_dsp PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen
	${CMAKE_CURRENT_SOURCE_DIR}/vendor/cpp-peglib
)
find_package(Threads REQUIRED)
target_link_libraries(triggerfish_dsp PUBLIC Threads::Threads)

# Give clangd accurate commands for the Rack-facing sources without making the
# standalone test build depend on the SDK. The target is indexing-only and is
# deliberately excluded from normal builds.
set(RACK_SDK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Rack-SDK" CACHE PATH
	"Path to the VCV Rack SDK used for editor indexing")
if(EXISTS "${RACK_SDK_DIR}/include/rack.hpp")
	add_library(triggerfish_rack_index OBJECT EXCLUDE_FROM_ALL
		src/plugin.cpp
		src/TfSlop.cpp
		src/TfSlop4.cpp
		src/TfVCA.cpp
		src/TfVDPO.cpp
		src/Tf303VoiceCore.cpp
		src/Tf303Oscillator.cpp
		src/Tf4072VoiceCore.cpp
		src/TfWavefoldOscillator.cpp
		src/TfUnisonOscillator.cpp
		src/TfScenePack4.cpp
		src/TfReverb.cpp
		src/TfProgSequencer.cpp
	)
	target_include_directories(triggerfish_rack_index PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/src/models
		${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen
		${CMAKE_CURRENT_SOURCE_DIR}/vendor/cpp-peglib
		${RACK_SDK_DIR}/include
		${RACK_SDK_DIR}/dep/include
	)
	target_compile_definitions(triggerfish_rack_index PRIVATE _USE_MATH_DEFINES)
	# Rack's Windows ABI packs Int24 into exactly three bytes. Clang otherwise
	# applies Microsoft bit-field layout while indexing the MinGW SDK headers.
	if(MINGW)
		target_compile_options(triggerfish_rack_index PRIVATE -mno-ms-bitfields)
	endif()
endif()

include(CTest)
if(BUILD_TESTING)
	add_executable(triggerfish_dsp_tests tests/dsp_tests.cpp)
	target_link_libraries(triggerfish_dsp_tests PRIVATE triggerfish_dsp)
	add_test(NAME triggerfish_dsp_tests COMMAND triggerfish_dsp_tests)

	add_executable(triggerfish_early_reflections_tests
		tests/early_reflections_tests.cpp)
	target_link_libraries(triggerfish_early_reflections_tests PRIVATE
		triggerfish_dsp)
	add_test(NAME triggerfish_early_reflections_tests
		COMMAND triggerfish_early_reflections_tests)

	add_executable(triggerfish_early_reflections_benchmark EXCLUDE_FROM_ALL
		tests/early_reflections_benchmark.cpp)
	target_link_libraries(triggerfish_early_reflections_benchmark PRIVATE
		triggerfish_dsp)

	add_executable(triggerfish_listener_position_analysis EXCLUDE_FROM_ALL
		tools/analyze_listener_position.cpp)
	target_link_libraries(triggerfish_listener_position_analysis PRIVATE
		triggerfish_dsp)

	add_executable(triggerfish_late_reverb_tests tests/late_reverb_tests.cpp)
	target_link_libraries(triggerfish_late_reverb_tests PRIVATE triggerfish_dsp)
	add_test(NAME triggerfish_late_reverb_tests
		COMMAND triggerfish_late_reverb_tests)

	add_executable(triggerfish_room_reverb_tests tests/room_reverb_tests.cpp)
	target_link_libraries(triggerfish_room_reverb_tests PRIVATE triggerfish_dsp)
	add_test(NAME triggerfish_room_reverb_tests
		COMMAND triggerfish_room_reverb_tests)

	add_executable(triggerfish_text_sequencer_tests
		tests/text_sequencer_tests.cpp)
	target_link_libraries(triggerfish_text_sequencer_tests PRIVATE
		triggerfish_dsp)
	add_test(NAME triggerfish_text_sequencer_tests
		COMMAND triggerfish_text_sequencer_tests)
endif()

option(TRIGGERFISH_BUILD_PYTHON "Build the optional Python DSP bindings" OFF)
if(TRIGGERFISH_BUILD_PYTHON)
	find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
	find_package(pybind11 CONFIG REQUIRED)
	pybind11_add_module(_triggerfish_dsp tests/python/bindings.cpp)
	target_link_libraries(_triggerfish_dsp PRIVATE triggerfish_dsp)
	if(MINGW)
		# Python 3.8+ does not search PATH for dependent extension DLLs. Keep
		# the development module self-contained instead of generating an MSVC
		# tree or copying MinGW runtime DLLs into the virtual environment.
		execute_process(
			COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libwinpthread.a
			OUTPUT_VARIABLE MINGW_STATIC_WINPTHREAD
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		if(NOT EXISTS "${MINGW_STATIC_WINPTHREAD}")
			message(FATAL_ERROR "MinGW static winpthread archive was not found")
		endif()
		target_link_options(_triggerfish_dsp PRIVATE
			-static-libgcc
			-static-libstdc++)
		target_link_libraries(_triggerfish_dsp PRIVATE
			"${MINGW_STATIC_WINPTHREAD}")
	endif()
	install(TARGETS _triggerfish_dsp LIBRARY DESTINATION . RUNTIME DESTINATION .)
endif()
