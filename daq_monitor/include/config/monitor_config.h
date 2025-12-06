#pragma once

#include <string>

struct MonitorConfig {
  // File monitoring settings
  std::string input_file;
  int polling_interval_ms = 1000;
  int display_update_interval_ms = 1000;
  int rate_window_seconds = 10;

  // QA settings
  bool qa_enabled = true;
  int qa_sampling_interval = 10;
  int qa_pedestal_samples = 100;
  float qa_baseline_target = 3500.0f;
  float qa_baseline_tolerance = 50.0f;
  float qa_noise_threshold = 10.0f;
  float qa_signal_min = -1000.0f;
  float qa_signal_max = 5000.0f;

  // Logging settings
  bool log_warnings = true;
  std::string log_file = "monitor.log";

  // Load configuration from JSON file
  static MonitorConfig LoadFromJson(const std::string& config_path);
};
