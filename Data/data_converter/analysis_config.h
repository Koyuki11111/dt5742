#pragma once

#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

struct AnalysisConfig {
  // Common fields (from pipeline_config.json "common" section)
  std::string output_dir = "output";
  int n_channels = 16;
  int max_cores = 8;
  int chunk_size = 100;
  std::string temp_dir = "./temp";
  std::string input_root = "waveforms.root";      // waveforms_root in config
  std::string input_tree = "Waveforms";            // waveforms_tree in config
  std::string output_root = "waveforms_analyzed.root";  // analysis_root in config
  std::string output_tree = "Analysis";            // analysis_tree in config

  // Stage2 specific fields

  // Overall analysis region (per channel, in nanoseconds)
  // Points outside this region will be ignored in all analysis
  std::vector<float> analysis_region_min;  // Start of valid analysis region (ns)
  std::vector<float> analysis_region_max;  // End of valid analysis region (ns)

  // Time regions for analysis (per channel, in nanoseconds)
  std::vector<float> baseline_region_min;  // Baseline calculation region start (ns)
  std::vector<float> baseline_region_max;  // Baseline calculation region end (ns)
  std::vector<float> signal_region_min;
  std::vector<float> signal_region_max;
  std::vector<float> charge_region_min;
  std::vector<float> charge_region_max;

  // Signal polarity (per channel): +1 for positive, -1 for negative
  std::vector<int> signal_polarity;

  // Signal detection threshold (SNR threshold)
  float snr_threshold = 3.0f;

  // CFD thresholds in percent (e.g., 10, 20, 30, 50)
  std::vector<int> cfd_thresholds = {10, 20, 30, 50};

  // Leading edge thresholds in mV
  std::vector<float> le_thresholds = {10.0f, 20.0f, 50.0f};

  // Charge thresholds in percent
  std::vector<int> charge_thresholds = {10, 20, 50};

  // Rise time calculation thresholds
  float rise_time_low = 0.1f;   // 10%
  float rise_time_high = 0.9f;  // 90%

  // Signal quality cuts (per channel)
  std::vector<float> cut_amp_max;

  // Impedance for charge calculation (Ohms)
  float impedance = 50.0f;

  // Waveform plots output options
  bool waveform_plots_enabled = false;
  std::string waveform_plots_dir = "waveform_plots";
  bool waveform_plots_only_signal = true;  // Only save waveforms with detected signal

  // Sensor mapping (per channel)
  std::vector<int> sensor_ids;  // Which sensor each channel belongs to
  std::vector<int> strip_ids;   // Strip number within sensor

  // Constructor with default values
  AnalysisConfig() {
    // Initialize per-channel vectors with default values
    analysis_region_min.assign(16, -100.0f);  // Default: use full waveform
    analysis_region_max.assign(16, 300.0f);
    baseline_region_min.assign(16, -50.0f);
    baseline_region_max.assign(16, -10.0f);
    signal_region_min.assign(16, 0.0f);
    signal_region_max.assign(16, 200.0f);
    charge_region_min.assign(16, 0.0f);
    charge_region_max.assign(16, 200.0f);
    cut_amp_max.assign(16, 1.0f);
    signal_polarity.assign(16, 1);  // Default: positive signals

    // Default sensor mapping: ch0-7 = sensor 1, ch8-15 = sensor 2
    sensor_ids = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2};
    strip_ids = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7};
  }
};

// Simple JSON parser for AnalysisConfig
inline bool ExtractFloatArray(const std::string &text, const std::string &key,
                              std::vector<float> &values) {
  const std::string pattern = "\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    std::string arrayContent = match[1].str();
    values.clear();
    std::regex numPattern("(-?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?)");
    auto numBegin = std::sregex_iterator(arrayContent.begin(), arrayContent.end(), numPattern);
    auto numEnd = std::sregex_iterator();
    for (auto it = numBegin; it != numEnd; ++it) {
      values.push_back(std::stof(it->str()));
    }
    return true;
  }
  return false;
}

inline bool ExtractIntArray(const std::string &text, const std::string &key,
                            std::vector<int> &values) {
  const std::string pattern = "\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    std::string arrayContent = match[1].str();
    values.clear();
    std::regex numPattern("(-?\\d+)");
    auto numBegin = std::sregex_iterator(arrayContent.begin(), arrayContent.end(), numPattern);
    auto numEnd = std::sregex_iterator();
    for (auto it = numBegin; it != numEnd; ++it) {
      values.push_back(std::stoi(it->str()));
    }
    return true;
  }
  return false;
}

inline bool ExtractStringValue(const std::string &text, const std::string &key,
                               std::string &value) {
  const std::string pattern =
      "\"" + key + "\"\\s*:\\s*\"((?:\\\\\"|[^\"])*)\"";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    value = match[1].str();
    return true;
  }
  return false;
}

inline bool ExtractNumberValue(const std::string &text, const std::string &key,
                               double &value) {
  const std::string pattern =
      "\"" + key + "\"\\s*:\\s*(-?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?)";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    value = std::stod(match[1].str());
    return true;
  }
  return false;
}

// Extract a section from JSON (e.g., "common", "stage1", "stage2")
inline bool ExtractSection(const std::string &text, const std::string &section,
                          std::string &sectionContent) {
  // Pattern: "section": { ... }
  // We need to handle nested braces properly
  std::string pattern = "\"" + section + "\"\\s*:\\s*\\{";
  std::regex re(pattern);
  std::smatch match;

  if (!std::regex_search(text, match, re)) {
    return false;
  }

  // Find the starting position of the section content
  size_t start = match.position() + match.length();

  // Count braces to find the end of this section
  int braceCount = 1;
  size_t end = start;
  bool inString = false;
  bool escaped = false;

  while (end < text.size() && braceCount > 0) {
    char c = text[end];

    if (escaped) {
      escaped = false;
      end++;
      continue;
    }

    if (c == '\\') {
      escaped = true;
      end++;
      continue;
    }

    if (c == '"') {
      inString = !inString;
    } else if (!inString) {
      if (c == '{') {
        braceCount++;
      } else if (c == '}') {
        braceCount--;
      }
    }

    end++;
  }

  if (braceCount != 0) {
    return false;  // Mismatched braces
  }

  // Extract the content (excluding the final '}')
  sectionContent = text.substr(start, end - start - 1);
  return true;
}

// Helper function to create directory recursively
inline bool CreateDirectoryIfNeeded(const std::string &path) {
  if (path.empty() || path == ".") {
    return true;  // Current directory always exists
  }

  struct stat info;
  if (stat(path.c_str(), &info) == 0) {
    if (info.st_mode & S_IFDIR) {
      return true;  // Directory already exists
    } else {
      std::cerr << "ERROR: Path exists but is not a directory: " << path << std::endl;
      return false;
    }
  }

  // Try to create directory
  // For nested paths, we need to create parent directories first
  size_t pos = path.find_last_of('/');
  if (pos != std::string::npos && pos > 0) {
    std::string parent = path.substr(0, pos);
    if (!CreateDirectoryIfNeeded(parent)) {
      return false;
    }
  }

  // Create the directory (mkdir with 0755 permissions)
  if (mkdir(path.c_str(), 0755) != 0) {
    std::cerr << "ERROR: Failed to create directory: " << path << std::endl;
    return false;
  }

  return true;
}

// Helper function to build output path with subdirectory
inline std::string BuildOutputPath(const std::string &output_dir,
                                   const std::string &subdir,
                                   const std::string &filename) {
  // If filename is absolute path, use it as-is
  if (!filename.empty() && filename[0] == '/') {
    return filename;
  }

  // If output_dir is empty or ".", just use filename
  if (output_dir.empty() || output_dir == ".") {
    return filename;
  }

  // Build path: output_dir/subdir/filename
  std::string path = output_dir;
  if (path.back() != '/') {
    path += '/';
  }
  path += subdir;
  if (path.back() != '/') {
    path += '/';
  }
  path += filename;

  return path;
}

inline bool LoadAnalysisConfigFromJson(const std::string &path,
                                      AnalysisConfig &cfg,
                                      std::string *errorMessage = nullptr) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    if (errorMessage) {
      *errorMessage = "cannot open config file: " + path;
    }
    return false;
  }

  std::string content((std::istreambuf_iterator<char>(fin)),
                      std::istreambuf_iterator<char>());

  if (content.empty()) {
    if (errorMessage) {
      *errorMessage = "config file is empty: " + path;
    }
    return false;
  }

  // Extract sections
  std::string commonSection, stage2Section;
  bool hasCommon = ExtractSection(content, "common", commonSection);
  bool hasStage2 = ExtractSection(content, "stage2", stage2Section);

  if (!hasCommon && !hasStage2) {
    if (errorMessage) {
      *errorMessage = "config file missing 'common' and 'stage2' sections";
    }
    return false;
  }

  // Parse common section
  if (hasCommon) {
    std::string strValue;
    double numValue = 0.0;

    if (ExtractStringValue(commonSection, "output_dir", strValue)) {
      cfg.output_dir = strValue;
    }
    if (ExtractNumberValue(commonSection, "n_channels", numValue)) {
      cfg.n_channels = static_cast<int>(numValue);
    }
    if (ExtractNumberValue(commonSection, "max_cores", numValue)) {
      cfg.max_cores = static_cast<int>(numValue);
    }
    if (ExtractNumberValue(commonSection, "chunk_size", numValue)) {
      cfg.chunk_size = static_cast<int>(numValue);
    }
    if (ExtractStringValue(commonSection, "temp_dir", strValue)) {
      cfg.temp_dir = strValue;
    }
    if (ExtractStringValue(commonSection, "waveforms_root", strValue)) {
      cfg.input_root = strValue;
    }
    if (ExtractStringValue(commonSection, "waveforms_tree", strValue)) {
      cfg.input_tree = strValue;
    }
    if (ExtractStringValue(commonSection, "analysis_root", strValue)) {
      cfg.output_root = strValue;
    }
    if (ExtractStringValue(commonSection, "analysis_tree", strValue)) {
      cfg.output_tree = strValue;
    }
  }

  // Parse stage2 section
  if (hasStage2) {
    std::string strValue;
    double numValue = 0.0;

    if (ExtractNumberValue(stage2Section, "rise_time_low", numValue)) {
      cfg.rise_time_low = static_cast<float>(numValue);
    }
    if (ExtractNumberValue(stage2Section, "rise_time_high", numValue)) {
      cfg.rise_time_high = static_cast<float>(numValue);
    }
    if (ExtractNumberValue(stage2Section, "impedance", numValue)) {
      cfg.impedance = static_cast<float>(numValue);
    }
    if (ExtractNumberValue(stage2Section, "snr_threshold", numValue)) {
      cfg.snr_threshold = static_cast<float>(numValue);
    }

    ExtractFloatArray(stage2Section, "analysis_region_min", cfg.analysis_region_min);
    ExtractFloatArray(stage2Section, "analysis_region_max", cfg.analysis_region_max);
    ExtractFloatArray(stage2Section, "baseline_region_min", cfg.baseline_region_min);
    ExtractFloatArray(stage2Section, "baseline_region_max", cfg.baseline_region_max);
    ExtractFloatArray(stage2Section, "signal_region_min", cfg.signal_region_min);
    ExtractFloatArray(stage2Section, "signal_region_max", cfg.signal_region_max);
    ExtractFloatArray(stage2Section, "charge_region_min", cfg.charge_region_min);
    ExtractFloatArray(stage2Section, "charge_region_max", cfg.charge_region_max);
    ExtractFloatArray(stage2Section, "cut_amp_max", cfg.cut_amp_max);
    ExtractFloatArray(stage2Section, "le_thresholds", cfg.le_thresholds);

    ExtractIntArray(stage2Section, "cfd_thresholds", cfg.cfd_thresholds);
    ExtractIntArray(stage2Section, "charge_thresholds", cfg.charge_thresholds);
    ExtractIntArray(stage2Section, "signal_polarity", cfg.signal_polarity);

    // Parse waveform plots output options
    std::string boolPattern = "\"waveform_plots_enabled\"\\s*:\\s*(true|false)";
    std::regex boolRe(boolPattern);
    std::smatch boolMatch;
    if (std::regex_search(stage2Section, boolMatch, boolRe)) {
      cfg.waveform_plots_enabled = (boolMatch[1].str() == "true");
    }

    std::string boolPattern2 = "\"waveform_plots_only_signal\"\\s*:\\s*(true|false)";
    std::regex boolRe2(boolPattern2);
    std::smatch boolMatch2;
    if (std::regex_search(stage2Section, boolMatch2, boolRe2)) {
      cfg.waveform_plots_only_signal = (boolMatch2[1].str() == "true");
    }

    if (ExtractStringValue(stage2Section, "waveform_plots_dir", strValue)) {
      cfg.waveform_plots_dir = strValue;
    }

    // Parse sensor mapping
    std::string sensorPattern = "\"sensor_mapping\"\\s*:\\s*\\{([^}]*)\\}";
    std::regex sensorRe(sensorPattern);
    std::smatch sensorMatch;
    if (std::regex_search(stage2Section, sensorMatch, sensorRe)) {
      std::string sensorContent = sensorMatch[1].str();
      ExtractIntArray(sensorContent, "sensor_ids", cfg.sensor_ids);
      ExtractIntArray(sensorContent, "strip_ids", cfg.strip_ids);
    }
  }

  // Ensure per-channel vectors have correct size
  if (cfg.analysis_region_min.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.analysis_region_min.resize(cfg.n_channels, -100.0f);
  }
  if (cfg.analysis_region_max.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.analysis_region_max.resize(cfg.n_channels, 300.0f);
  }
  if (cfg.baseline_region_min.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.baseline_region_min.resize(cfg.n_channels, -50.0f);
  }
  if (cfg.baseline_region_max.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.baseline_region_max.resize(cfg.n_channels, -10.0f);
  }
  if (cfg.signal_region_min.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.signal_region_min.resize(cfg.n_channels, 0.0f);
  }
  if (cfg.signal_region_max.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.signal_region_max.resize(cfg.n_channels, 200.0f);
  }
  if (cfg.charge_region_min.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.charge_region_min.resize(cfg.n_channels, 0.0f);
  }
  if (cfg.charge_region_max.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.charge_region_max.resize(cfg.n_channels, 200.0f);
  }
  if (cfg.cut_amp_max.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.cut_amp_max.resize(cfg.n_channels, 1.0f);
  }
  if (cfg.signal_polarity.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.signal_polarity.resize(cfg.n_channels, 1);
  }

  // Ensure sensor mapping vectors have correct size
  if (cfg.sensor_ids.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.sensor_ids.resize(cfg.n_channels);
    cfg.strip_ids.resize(cfg.n_channels);
    // Set default mapping if not specified
    for (int i = 0; i < cfg.n_channels; ++i) {
      cfg.sensor_ids[i] = (i < 8) ? 1 : 2;
      cfg.strip_ids[i] = i % 8;
    }
  }
  if (cfg.strip_ids.size() < static_cast<size_t>(cfg.n_channels)) {
    cfg.strip_ids.resize(cfg.n_channels);
    for (int i = 0; i < cfg.n_channels; ++i) {
      cfg.strip_ids[i] = i % 8;
    }
  }

  return true;
}
