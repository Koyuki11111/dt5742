#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "utils/filesystem_utils.h"
#include "utils/json_utils.h"
#include "config/common_config.h"

struct AnalysisConfig {
  // Common fields (from converter_config.json "common" section)
  CommonConfig common;
  // Accessors for shared fields
  std::string output_dir() const { return common.output_dir; }
  std::string input_dir()  const { return common.input_dir; }
  std::string daq_name()   const { return common.daq_name; }
  int runnumber()  const { return common.runnumber; }
  int n_channels() const { return common.n_channels; }
  int n_sensors()  const { return common.n_sensors; }
  int max_cores()  const { return common.max_cores; }
  int chunk_size() const { return common.chunk_size; }
  std::string temp_dir() const { return common.temp_dir; }
  std::string input_root() const { return common.waveforms_root; }
  std::string input_tree() const { return common.waveforms_tree; }
  std::string output_root() const { return common.analysis_root; }
  std::string output_tree() const { return common.analysis_tree; }
  void set_n_channels(int v) { common.n_channels = v; }
  void set_n_sensors(int v)  { common.n_sensors  = v; }
  void set_input_root(const std::string &v) { common.waveforms_root = v; }
  void set_output_root(const std::string &v) { common.analysis_root = v; }
  void set_input_tree(const std::string &v) { common.waveforms_tree = v; }
  void set_output_tree(const std::string &v) { common.analysis_tree = v; }

  // waveform_analyzer specific fields

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
  std::vector<int> sensor_cols; // Strip IDs (Unified name)
  std::vector<int> sensor_rows; // Column IDs
  std::vector<std::string> sensor_orientations;  // "vertical" or "horizontal" per sensor
  
  // Constructor with default values
  AnalysisConfig() {
    // Initialize per-channel vectors with default values
    analysis_region_min.assign(common.n_channels, -100.0f);  // Default: use full waveform
    analysis_region_max.assign(common.n_channels, 300.0f);
    baseline_region_min.assign(common.n_channels, -50.0f);
    baseline_region_max.assign(common.n_channels, -10.0f);
    signal_region_min.assign(common.n_channels, 0.0f);
    signal_region_max.assign(common.n_channels, 200.0f);
    charge_region_min.assign(common.n_channels, 0.0f);
    charge_region_max.assign(common.n_channels, 200.0f);
    cut_amp_max.assign(common.n_channels, 1.0f);
    signal_polarity.assign(common.n_channels, 1);  // Default: positive signals

    // Default sensor mapping: ch0-7 = sensor 1, ch8-15 = sensor 2
    sensor_ids      = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2};
    // Default strips 0-7
    sensor_cols     = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7};
    // Default columns (0)
    sensor_rows     = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    sensor_orientations = {"vertical", "vertical"};  // Default: all vertical (local sensors)
  }
};

inline bool LoadAnalysisConfigFromJson(const std::string &path,
                                      AnalysisConfig &cfg,
                                      std::string *errorMessage = nullptr) {
  simdjson::dom::parser parser;
  simdjson::dom::element root;
  if (!ParseJsonFile(path, parser, root, errorMessage)) {
    return false;
  }

  simdjson::dom::element common;
  if (GetObject(root, "common", common)) {
    std::string strValue;
    double numValue = 0.0;

    if (GetString(common, "output_dir", strValue)) {
      cfg.common.output_dir = strValue;
    }    
    if (GetString(common, "input_dir", strValue)) {
      cfg.common.input_dir = strValue;
    }    
    if (GetString(common, "daq_name", strValue)) {
      cfg.common.daq_name = strValue;
    }    
    if (GetNumber(common, "runnumber", numValue)) {
      cfg.common.runnumber = static_cast<int>(numValue);
    }    
    if (GetNumber(common, "n_channels", numValue)) {
      cfg.common.n_channels = static_cast<int>(numValue);
    }
    if (GetNumber(common, "n_sensors", numValue)) {
      cfg.common.n_sensors = static_cast<int>(numValue);
    }
    if (GetNumber(common, "max_cores", numValue)) {
      cfg.common.max_cores = static_cast<int>(numValue);
    }
    if (GetNumber(common, "chunk_size", numValue)) {
      cfg.common.chunk_size = static_cast<int>(numValue);
    }
    if (GetString(common, "temp_dir", strValue)) {
      cfg.common.temp_dir = strValue;
    }
    if (GetString(common, "nsamples_policy", strValue)) {
      cfg.common.nsamples_policy = strValue;
    }
    if (GetString(common, "waveforms_root", strValue)) {
      cfg.common.waveforms_root = strValue;
    }
    if (GetString(common, "waveforms_tree", strValue)) {
      cfg.common.waveforms_tree = strValue;
    }
    if (GetString(common, "analysis_root", strValue)) {
      cfg.common.analysis_root = strValue;
    }
    if (GetString(common, "analysis_tree", strValue)) {
      cfg.common.analysis_tree = strValue;
    }
  }

  simdjson::dom::element waveformAnalyzer;
  if (GetObject(root, "waveform_analyzer", waveformAnalyzer)) {
    std::string strValue;
    double numValue = 0.0;

    if (GetNumber(waveformAnalyzer, "rise_time_low", numValue)) {
      cfg.rise_time_low = static_cast<float>(numValue);
    }
    if (GetNumber(waveformAnalyzer, "rise_time_high", numValue)) {
      cfg.rise_time_high = static_cast<float>(numValue);
    }
    if (GetNumber(waveformAnalyzer, "impedance", numValue)) {
      cfg.impedance = static_cast<float>(numValue);
    }
    if (GetNumber(waveformAnalyzer, "snr_threshold", numValue)) {
      cfg.snr_threshold = static_cast<float>(numValue);
    }

    GetFloatArray(waveformAnalyzer, "analysis_region_min", cfg.analysis_region_min);
    GetFloatArray(waveformAnalyzer, "analysis_region_max", cfg.analysis_region_max);
    GetFloatArray(waveformAnalyzer, "baseline_region_min", cfg.baseline_region_min);
    GetFloatArray(waveformAnalyzer, "baseline_region_max", cfg.baseline_region_max);
    GetFloatArray(waveformAnalyzer, "signal_region_min", cfg.signal_region_min);
    GetFloatArray(waveformAnalyzer, "signal_region_max", cfg.signal_region_max);
    GetFloatArray(waveformAnalyzer, "charge_region_min", cfg.charge_region_min);
    GetFloatArray(waveformAnalyzer, "charge_region_max", cfg.charge_region_max);
    GetFloatArray(waveformAnalyzer, "cut_amp_max", cfg.cut_amp_max);
    GetFloatArray(waveformAnalyzer, "le_thresholds", cfg.le_thresholds);

    GetIntArray(waveformAnalyzer, "cfd_thresholds", cfg.cfd_thresholds);
    GetIntArray(waveformAnalyzer, "charge_thresholds", cfg.charge_thresholds);
    GetIntArray(waveformAnalyzer, "signal_polarity", cfg.signal_polarity);

    bool boolValue = false;
    if (GetBool(waveformAnalyzer, "waveform_plots_enabled", boolValue)) {
      cfg.waveform_plots_enabled = boolValue;
    }
    if (GetBool(waveformAnalyzer, "waveform_plots_only_signal", boolValue)) {
      cfg.waveform_plots_only_signal = boolValue;
    }

    if (GetString(waveformAnalyzer, "waveform_plots_dir", strValue)) {
      cfg.waveform_plots_dir = strValue;
    }

    simdjson::dom::element sensorSection;
    if (GetObject(waveformAnalyzer, "sensor_mapping", sensorSection)) {
      GetIntArray(sensorSection, "sensor_ids", cfg.sensor_ids);
      // Unify: Read strip_ids into sensor_cols
      GetIntArray(sensorSection, "strip_ids", cfg.sensor_cols);
      // Map column_ids to sensor_rows
      GetIntArray(sensorSection, "column_ids", cfg.sensor_rows);
      
      // Parse sensor orientations
      simdjson::dom::array orientationsArray;
      if (sensorSection["sensor_orientations"].get(orientationsArray) == simdjson::SUCCESS) {
        cfg.sensor_orientations.clear();
        for (auto val : orientationsArray) {
          std::string_view sv;
          if (val.get(sv) == simdjson::SUCCESS) {
            cfg.sensor_orientations.push_back(std::string(sv));
          }
        }
      }
    }
  }

  // Ensure per-channel vectors have correct size
  if (cfg.analysis_region_min.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.analysis_region_min.resize(cfg.common.n_channels, -100.0f);
  }
  if (cfg.analysis_region_max.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.analysis_region_max.resize(cfg.common.n_channels, 300.0f);
  }
  if (cfg.baseline_region_min.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.baseline_region_min.resize(cfg.common.n_channels, -50.0f);
  }
  if (cfg.baseline_region_max.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.baseline_region_max.resize(cfg.common.n_channels, -10.0f);
  }
  if (cfg.signal_region_min.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.signal_region_min.resize(cfg.common.n_channels, 0.0f);
  }
  if (cfg.signal_region_max.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.signal_region_max.resize(cfg.common.n_channels, 200.0f);
  }
  if (cfg.charge_region_min.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.charge_region_min.resize(cfg.common.n_channels, 0.0f);
  }
  if (cfg.charge_region_max.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.charge_region_max.resize(cfg.common.n_channels, 200.0f);
  }
  if (cfg.cut_amp_max.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.cut_amp_max.resize(cfg.common.n_channels, 1.0f);
  }
  if (cfg.signal_polarity.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.signal_polarity.resize(cfg.common.n_channels, 1);
  }

  // Ensure sensor mapping vectors have correct size
  if (cfg.sensor_ids.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.sensor_ids.resize(cfg.common.n_channels);
    cfg.sensor_cols.resize(cfg.common.n_channels);
    cfg.sensor_rows.resize(cfg.common.n_channels);
    // Set default mapping if not specified
    for (int i = 0; i < cfg.common.n_channels; ++i) {
      cfg.sensor_ids[i] = (i < 8) ? 1 : 2;
      cfg.sensor_cols[i] = i % 8;
    }
  }
  if (cfg.sensor_cols.size() < static_cast<size_t>(cfg.common.n_channels)) {
    cfg.sensor_cols.resize(cfg.common.n_channels);
    for (int i = 0; i < cfg.common.n_channels; ++i) {
      cfg.sensor_cols[i] = i % 8;
    }
  }

  // Validate sensor orientations (typically 2 entries for local sensors)
  if (cfg.sensor_orientations.empty()) {
    cfg.sensor_orientations.resize(2, "vertical");
  }

  // Validate values and convert to lowercase
  for (size_t i = 0; i < cfg.sensor_orientations.size(); ++i) {
    std::string& orient = cfg.sensor_orientations[i];
    std::transform(orient.begin(), orient.end(), orient.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (orient != "vertical" && orient != "horizontal") {
      std::cerr << "WARNING: Invalid orientation '" << orient
                << "' for sensor " << (i+1) << ", defaulting to 'vertical'" << std::endl;
      orient = "vertical";
    }
  }

  return true;
}
