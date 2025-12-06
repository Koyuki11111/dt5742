#include "config/monitor_config.h"
#include "monitor/realtime_monitor.h"

#include <csignal>
#include <iostream>
#include <string>

// Global monitor instance for signal handler
RealtimeMonitor* g_monitor = nullptr;

// Signal handler for Ctrl+C
void SignalHandler(int signal) {
  if (signal == SIGINT) {
    std::cout << "\n\nReceived Ctrl+C, stopping monitor...\n";
    if (g_monitor) {
      g_monitor->Stop();
    }
  }
}

void PrintUsage(const char* program_name) {
  std::cout << "\nCAEN DT5742 Real-Time Data Monitor\n\n";
  std::cout << "Usage:\n";
  std::cout << "  " << program_name << " [OPTIONS]\n\n";
  std::cout << "Options:\n";
  std::cout << "  --config FILE       Path to configuration file (default: monitor_config.json)\n";
  std::cout << "  --file FILE         Path to input binary file (overrides config)\n";
  std::cout << "  --no-qa             Disable QA checks (header-only monitoring)\n";
  std::cout << "  --help              Display this help message\n\n";
  std::cout << "Examples:\n";
  std::cout << "  # Monitor with default config\n";
  std::cout << "  " << program_name << "\n\n";
  std::cout << "  # Monitor specific file\n";
  std::cout << "  " << program_name << " --file ../Data/AC_LGAD_TEST/wave_0.dat\n\n";
  std::cout << "  # Monitor without QA checks (faster)\n";
  std::cout << "  " << program_name << " --no-qa\n\n";
  std::cout << "Signals:\n";
  std::cout << "  Ctrl+C              Stop monitoring and print summary\n\n";
}

int main(int argc, char* argv[]) {
  // Default configuration file
  std::string config_file = "monitor_config.json";
  std::string override_file;
  bool disable_qa = false;

  // Parse command-line arguments
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      return 0;
    } else if (arg == "--config" && i + 1 < argc) {
      config_file = argv[++i];
    } else if (arg == "--file" && i + 1 < argc) {
      override_file = argv[++i];
    } else if (arg == "--no-qa") {
      disable_qa = true;
    } else {
      std::cerr << "Unknown option: " << arg << std::endl;
      std::cerr << "Use --help for usage information" << std::endl;
      return 1;
    }
  }

  // Load configuration
  MonitorConfig config = MonitorConfig::LoadFromJson(config_file);

  // Apply command-line overrides
  if (!override_file.empty()) {
    config.input_file = override_file;
  }

  if (disable_qa) {
    config.qa_enabled = false;
  }

  // Validate configuration
  if (config.input_file.empty()) {
    std::cerr << "Error: No input file specified in configuration or command line" << std::endl;
    std::cerr << "Use --file to specify an input file" << std::endl;
    return 1;
  }

  std::cout << "CAEN DT5742 Real-Time Monitor\n";
  std::cout << "Configuration: " << config_file << "\n";
  std::cout << "Input file: " << config.input_file << "\n";
  std::cout << "QA enabled: " << (config.qa_enabled ? "yes" : "no") << "\n";

  // Create monitor
  RealtimeMonitor monitor(config);
  g_monitor = &monitor;

  // Setup signal handler
  std::signal(SIGINT, SignalHandler);

  // Initialize and run
  if (!monitor.Initialize()) {
    std::cerr << "Failed to initialize monitor" << std::endl;
    return 1;
  }

  monitor.Run();

  // Cleanup
  g_monitor = nullptr;

  return 0;
}
