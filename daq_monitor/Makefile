# Makefile for CAEN DT5742 Real-Time Monitor

CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
SRCDIR = src
INCLUDES = -Iinclude -I. -I/opt/homebrew/include
JSON_LIBS = -L/opt/homebrew/lib -lsimdjson

# Target
TARGET = monitor_realtime

# Source files
SOURCES = $(SRCDIR)/monitor_realtime.cpp \
          $(SRCDIR)/monitor/realtime_monitor.cpp \
          $(SRCDIR)/config/monitor_config.cpp \
          $(SRCDIR)/utils/file_io.cpp

# Header files (for dependency tracking)
HEADERS = include/monitor/realtime_monitor.h \
          include/config/monitor_config.h \
          include/utils/file_io.h \
          include/utils/json_utils.h

# Default target
all: $(TARGET)

# Build monitor
$(TARGET): $(SOURCES) $(HEADERS)
	@echo "Building $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(SOURCES) $(JSON_LIBS)
	@echo "Build successful!"

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGET)
	rm -f *.o
	rm -f monitor.log

# Help target
help:
	@echo "CAEN DT5742 Real-Time Monitor - Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build monitor executable"
	@echo "  make clean        - Remove executable and log files"
	@echo "  make help         - Show this help"
	@echo ""
	@echo "Prerequisites:"
	@echo "  - C++17 compiler (g++ or clang++)"
	@echo "  - simdjson library: brew install simdjson"
	@echo ""
	@echo "Usage:"
	@echo "  1. Build: make"
	@echo "  2. Run:   ./monitor_realtime --file ../Data/AC_LGAD_TEST/wave_0.dat"
	@echo "  or with config:"
	@echo "    ./monitor_realtime --config monitor_config.json"
	@echo ""
	@echo "Command-line options:"
	@echo "  --config FILE     Configuration file (default: monitor_config.json)"
	@echo "  --file FILE       Input binary file (overrides config)"
	@echo "  --no-qa           Disable QA checks (faster)"
	@echo "  --help            Show help message"

.PHONY: all clean help
