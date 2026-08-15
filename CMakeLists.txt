cmake_minimum_required(VERSION 3.25)

project(TriggerFishDSP VERSION 2.2.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(triggerfish_dsp STATIC
	src/tfdsp/sampleRate.cpp
)
set_target_properties(triggerfish_dsp PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(triggerfish_dsp PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/dep/eigen
)

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
		src/TfDiodeLadderFilter.cpp
		src/Tf303Oscillator.cpp
	)
	target_include_directories(triggerfish_rack_index PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/src/models
		${CMAKE_CURRENT_SOURCE_DIR}/dep/eigen
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
endif()

option(TRIGGERFISH_BUILD_PYTHON "Build the optional Python DSP bindings" OFF)
if(TRIGGERFISH_BUILD_PYTHON)
	find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
	find_package(pybind11 CONFIG REQUIRED)
	pybind11_add_module(_triggerfish_dsp tests/python/bindings.cpp)
	target_link_libraries(_triggerfish_dsp PRIVATE triggerfish_dsp)
	install(TARGETS _triggerfish_dsp LIBRARY DESTINATION . RUNTIME DESTINATION .)
endif()
