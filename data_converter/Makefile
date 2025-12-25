# Makefile for waveform processing pipeline

CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
SRCDIR = src
UNAME_S := $(shell uname -s)

# Base include paths
BASE_INCLUDES = -Iinclude -I.

# simdjson configuration (prefer pkg-config, fallback to macOS Homebrew defaults)
JSON_LIBS ?= $(shell pkg-config --libs simdjson 2>/dev/null)
JSON_INCLUDES ?= $(shell pkg-config --cflags simdjson 2>/dev/null)
ifeq ($(JSON_LIBS),)
	JSON_INCLUDES ?= -I/opt/homebrew/include
	JSON_LIBS ?= -L/opt/homebrew/lib -lsimdjson
endif
INCLUDES = $(BASE_INCLUDES) $(JSON_INCLUDES)

# ROOT configuration
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

# HDF5 configuration
HDF5_LIBS ?= $(shell pkg-config --libs hdf5 2>/dev/null)
HDF5_CFLAGS ?= $(shell pkg-config --cflags hdf5 2>/dev/null)
ifeq ($(HDF5_LIBS),)
	HDF5_PREFIX ?= /opt/homebrew/opt/hdf5
	ifeq ($(UNAME_S),Linux)
		HDF5_PREFIX ?= /usr
	endif
	HDF5_CFLAGS = -I$(HDF5_PREFIX)/include
	HDF5_LIBS = -L$(HDF5_PREFIX)/lib -lhdf5
endif

# Targets
TARGETS = convert_to_root analyze_waveforms export_to_hdf5 fast_qa

# Default target
all: $(TARGETS) parallel_analyze.sh qa_comparison

# Stage 1: Convert binary/ASCII to ROOT
convert_to_root: $(SRCDIR)/convert_to_root.cpp include/config/wave_converter_config.h $(SRCDIR)/utils/file_io.cpp include/utils/file_io.h
	@echo "Building convert_to_root..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ROOT_CFLAGS) -o $@ $(SRCDIR)/convert_to_root.cpp $(SRCDIR)/utils/file_io.cpp $(ROOT_LIBS) $(JSON_LIBS)

# Stage 2: Analyze waveforms
analyze_waveforms: $(SRCDIR)/analyze_waveforms.cpp include/config/analysis_config.h $(SRCDIR)/analysis/waveform_math.cpp include/analysis/waveform_math.h $(SRCDIR)/analysis/waveform_plotting.cpp include/analysis/waveform_plotting.h
	@echo "Building analyze_waveforms..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ROOT_CFLAGS) -o $@ $(SRCDIR)/analyze_waveforms.cpp $(SRCDIR)/analysis/waveform_math.cpp $(SRCDIR)/analysis/waveform_plotting.cpp $(ROOT_LIBS) $(JSON_LIBS)

# Stage 3: Export to HDF5
export_to_hdf5: $(SRCDIR)/export_to_hdf5.cpp
	@echo "Building export_to_hdf5..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ROOT_CFLAGS) $(HDF5_CFLAGS) -o $@ $< $(ROOT_LIBS) $(HDF5_LIBS) $(JSON_LIBS)

# Fast QA: Generate quality check plots
fast_qa: $(SRCDIR)/fast_qa.cpp include/config/analysis_config.h
	@echo "Building fast_qa..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ROOT_CFLAGS) -o $@ $(SRCDIR)/fast_qa.cpp $(ROOT_LIBS) $(JSON_LIBS)

# Utility object (optional for reuse)
utils/file_io.o: $(SRCDIR)/utils/file_io.cpp include/utils/file_io.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/utils/file_io.cpp -o $@

# Ensure parallel_analyze.sh is executable
parallel_analyze.sh:
	@if [ -f parallel_analyze.sh ]; then \
		chmod +x parallel_analyze.sh; \
		echo "Made parallel_analyze.sh executable"; \
	fi

# Ensure qa_comparison is executable
qa_comparison:
	@if [ -f src/qa_comparison ]; then \
		chmod +x src/qa_comparison; \
		echo "Made src/qa_comparison executable"; \
	else \
		echo "WARNING: src/qa_comparison not found"; \
	fi

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGETS)
	rm -f *.o

# Clean all generated files (including output data)
cleanall: clean
	@echo "Cleaning all generated files..."
	rm -f *.root *.hdf5 *.h5

# Help target
help:
	@echo "Waveform Processing Pipeline - Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build all executables"
	@echo "  make clean        - Remove executables"
	@echo "  make cleanall     - Remove executables and output files"
	@echo "  make test         - Run test pipeline"
	@echo "  make help         - Show this help"
	@echo ""
	@echo "Prerequisites:"
	@echo "  - ROOT 6.x:   source /opt/root6/root_install/bin/thisroot.sh"
	@echo "  - HDF5:       pkg-config hdf5 or set HDF5_PREFIX (/usr or /opt/homebrew/opt/hdf5)"
	@echo "  - simdjson:   pkg-config simdjson or set JSON_INCLUDES/JSON_LIBS"
	@echo "  - C++17 compiler"
	@echo ""
	@echo "Individual stages:"
	@echo "  make convert_to_root     - Build stage 1 (binary/ASCII to ROOT)"
	@echo "  make analyze_waveforms   - Build stage 2 (waveform analysis)"
	@echo "  make export_to_hdf5      - Build stage 3 (ROOT to HDF5)"
	@echo "  make fast_qa             - Build fast QA tool"
	@echo ""
	@echo "Usage:"
	@echo "  1. Build: make"
	@echo "  2. Run:   ./run_full_pipeline.sh"
	@echo "  or run stages individually:"
	@echo "    ./convert_to_root --config wave_converter_config.json"
	@echo "    ./analyze_waveforms --config analysis_config.json"
	@echo "    ./export_to_hdf5 --mode raw --input waveforms.root --output waveforms.h5"
	@echo "    ./fast_qa --config converter_config.json"
	@echo ""
	@echo "QA Tools:"
	@echo "  ./src/qa_comparison <run_number> [--num-events N]  - Compare ROOT and HDF5 QA plots"

# Test target - run a quick test of the pipeline
test: $(TARGETS)
	@echo "Running test pipeline..."
	@echo "Note: This requires a tiny input dataset configured in converter_config.json"
	@if [ ! -f converter_config.json ]; then \
		echo "ERROR: converter_config.json not found"; \
		exit 1; \
	fi
	@./run_full_pipeline.sh

.PHONY: all clean cleanall help test
