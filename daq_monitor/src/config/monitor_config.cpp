#include "config/monitor_config.h"
#include "utils/json_utils.h"

#include <iostream>

MonitorConfig MonitorConfig::LoadFromJson(const std::string& config_path) {
  MonitorConfig config;

  simdjson::dom::parser parser;
  simdjson::dom::element root;
  std::string error_message;

  if (!ParseJsonFile(config_path, parser, root, &error_message)) {
    std::cerr << "Error loading monitor config: " << error_message << std::endl;
    std::cerr << "Using default configuration" << std::endl;
    return config;
  }

  // Get monitor section
  simdjson::dom::element monitor_section;
  if (!GetObject(root, "monitor", monitor_section)) {
    std::cerr << "Warning: 'monitor' section not found in config, using defaults" << std::endl;
    return config;
  }

  // File monitoring settings
  std::string input_file;
  if (GetString(monitor_section, "input_file", input_file)) {
    config.input_file = input_file;
  }

  double temp_num = 0.0;
  if (GetNumber(monitor_section, "polling_interval_ms", temp_num)) {
    config.polling_interval_ms = static_cast<int>(temp_num);
  }

  if (GetNumber(monitor_section, "display_update_interval_ms", temp_num)) {
    config.display_update_interval_ms = static_cast<int>(temp_num);
  }

  if (GetNumber(monitor_section, "rate_window_seconds", temp_num)) {
    config.rate_window_seconds = static_cast<int>(temp_num);
  }

  // QA settings
  bool qa_enabled;
  if (GetBool(monitor_section, "qa_enabled", qa_enabled)) {
    config.qa_enabled = qa_enabled;
  }

  if (GetNumber(monitor_section, "qa_sampling_interval", temp_num)) {
    config.qa_sampling_interval = static_cast<int>(temp_num);
  }

  if (GetNumber(monitor_section, "qa_pedestal_samples", temp_num)) {
    config.qa_pedestal_samples = static_cast<int>(temp_num);
  }

  if (GetNumber(monitor_section, "qa_baseline_target", temp_num)) {
    config.qa_baseline_target = static_cast<float>(temp_num);
  }

  if (GetNumber(monitor_section, "qa_baseline_tolerance", temp_num)) {
    config.qa_baseline_tolerance = static_cast<float>(temp_num);
  }

  if (GetNumber(monitor_section, "qa_noise_threshold", temp_num)) {
    config.qa_noise_threshold = static_cast<float>(temp_num);
  }

  if (GetNumber(monitor_section, "qa_signal_min", temp_num)) {
    config.qa_signal_min = static_cast<float>(temp_num);
  }

  if (GetNumber(monitor_section, "qa_signal_max", temp_num)) {
    config.qa_signal_max = static_cast<float>(temp_num);
  }

  // Logging settings
  bool log_warnings;
  if (GetBool(monitor_section, "log_warnings", log_warnings)) {
    config.log_warnings = log_warnings;
  }

  std::string log_file;
  if (GetString(monitor_section, "log_file", log_file)) {
    config.log_file = log_file;
  }

  return config;
}
