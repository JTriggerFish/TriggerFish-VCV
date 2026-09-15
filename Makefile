# Must follow the format in the Naming section of
# https://vcvrack.com/manual/PluginDevelopmentTutorial.html
SLUG = TriggerFish-Elements

# Must follow the format in the Versioning section of
# https://vcvrack.com/manual/PluginDevelopmentTutorial.html
VERSION = 2.5.4


# FLAGS are passed to both the C and C++ compiler by the Rack SDK.
FLAGS += -std=c++17 -Isrc -Isrc/models -Ivendor/eigen -Ivendor/cpp-peglib

# Careful about linking to shared libraries, since you can't assume much about the user's environment and library search path.
# Static libraries are fine.
LDFLAGS +=

# cpp-peglib uses std::call_once internally. Some Linux toolchains require the
# pthread compile/link flag explicitly even though no parser work runs on the
# audio thread.
ifdef ARCH_LIN
FLAGS += -pthread
LDFLAGS += -pthread
endif

# Add .cpp and .c files to the build
SOURCES += $(wildcard src/*.cpp src/tfdsp/*.cpp src/models/*.cpp)

# Add files to the ZIP package when running `make dist`
# The compiled plugin is automatically added.
DISTRIBUTABLES += $(wildcard LICENSE*) res vendor/cpp-peglib/LICENSE

# Include the VCV Rack plugin Makefile framework
include $(RACK_DIR)/plugin.mk
ifdef ARCH_MAC
	#CC=gcc-8
	#CXX=g++-8
endif
