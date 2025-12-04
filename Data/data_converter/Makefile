# Makefile for waveform processing pipeline

CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# ROOT configuration
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

# HDF5 configuration
HDF5_PREFIX = /opt/homebrew/opt/hdf5
HDF5_CFLAGS = -I$(HDF5_PREFIX)/include
HDF5_LIBS = -L$(HDF5_PREFIX)/lib -lhdf5

# Targets
TARGETS = convert_to_root analyze_waveforms export_to_hdf5

# Default target
all: $(TARGETS) parallel_analyze.sh

# Stage 1: Convert binary/ASCII to ROOT
convert_to_root: convert_to_root.cpp wave_converter_config.h
	@echo "Building convert_to_root..."
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) -o $@ $< $(ROOT_LIBS)

# Stage 2: Analyze waveforms
analyze_waveforms: analyze_waveforms.cpp analysis_config.h
	@echo "Building analyze_waveforms..."
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) -o $@ $< $(ROOT_LIBS)

# Stage 3: Export to HDF5
export_to_hdf5: export_to_hdf5.cpp
	@echo "Building export_to_hdf5..."
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) $(HDF5_CFLAGS) -o $@ $< $(ROOT_LIBS) $(HDF5_LIBS)

# Ensure parallel_analyze.sh is executable
parallel_analyze.sh:
	@if [ -f parallel_analyze.sh ]; then \
		chmod +x parallel_analyze.sh; \
		echo "Made parallel_analyze.sh executable"; \
	fi

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGETS)
	rm -f *.o

# Clean all generated files (including output data)
cleanall: clean
	@echo "Cleaning all generated files..."
	rm -f *.root *.h5

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
	@echo "  - HDF5:       Installed at $(HDF5_PREFIX)"
	@echo "  - C++17 compiler"
	@echo ""
	@echo "Individual stages:"
	@echo "  make convert_to_root     - Build stage 1 (binary/ASCII to ROOT)"
	@echo "  make analyze_waveforms   - Build stage 2 (waveform analysis)"
	@echo "  make export_to_hdf5      - Build stage 3 (ROOT to HDF5)"
	@echo ""
	@echo "Usage:"
	@echo "  1. Build: make"
	@echo "  2. Run:   ./run_full_pipeline.sh"
	@echo "  or run stages individually:"
	@echo "    ./convert_to_root --config wave_converter_config.json"
	@echo "    ./analyze_waveforms --config analysis_config.json"
	@echo "    ./export_to_hdf5 --mode raw --input waveforms.root --output waveforms.h5"

# Test target - run a quick test of the pipeline
test: $(TARGETS)
	@echo "Running test pipeline..."
	@echo "Note: This requires valid input files configured in wave_converter_config.json"
	@if [ ! -f wave_converter_config.json ]; then \
		echo "ERROR: wave_converter_config.json not found"; \
		exit 1; \
	fi
	@./run_full_pipeline.sh

.PHONY: all clean cleanall help test
