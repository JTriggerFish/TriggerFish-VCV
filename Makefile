# If RACK_DIR is not defined when calling the Makefile, default to two directories above
RACK_DIR ?= ../..

# FLAGS will be passed to both the C and C++ compiler
FLAGS += -std=c++17 #-faligned-allocation -faligned-new# -fopt-info-loop-optimized # -fopt-info-vec-missed
CFLAGS +=
CXXFLAGS += -Isrc -Isrc/dsp -Isrc/models

# Careful about linking to shared libraries, since you can't assume much about the user's environment and library search path.
# Static libraries are fine.
LDFLAGS +=

# Add .cpp and .c files to the build
SOURCES += $(wildcard src/*.cpp src/dsp/*.cpp src/models/*.cpp)

# Add files to the ZIP package when running `make dist`
# The compiled plugin is automatically added.
DISTRIBUTABLES += $(wildcard LICENSE*) res

# Include the VCV Rack plugin Makefile framework
include $(RACK_DIR)/plugin.mk
ifdef ARCH_MAC
	#CC=gcc-8
	#CXX=g++-8
endif
