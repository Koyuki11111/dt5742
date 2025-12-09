#pragma once

#include <string>
#include <iostream>

#include "utils/filesystem_utils.h"
#include "utils/json_utils.h"
#include "config/common_config.h"

struct WaveConverterConfig {
  // Common fields (from converter_config.json "common" section)
  CommonConfig common;

  // Accessors for common fields
  std::string output_dir() const { return common.output_dir; }
  int n_channels() const { return common.n_channels; }
  int max_cores() const { return common.max_cores; }
  int max_events() const { return common.max_events; }
  int chunk_size() const { return common.chunk_size; }
  std::string temp_dir() const { return common.temp_dir; }
  std::string root_file() const { return common.waveforms_root; }
  std::string tree_name() const { return common.waveforms_tree; }
  void set_n_channels(int v) { common.n_channels = v; }
  void set_max_cores(int v) { common.max_cores = v; }
  void set_max_events(int v) { common.max_events = v; }
  void set_chunk_size(int v) { common.chunk_size = v; }
  void set_root_file(const std::string &v) { common.waveforms_root = v; }
  void set_tree_name(const std::string &v) { common.waveforms_tree = v; }
  // waveform_converter specific fields
  std::string input_pattern = "wave_%d.dat";
  std::string input_dir = ".";
  bool input_is_ascii = false;
  std::string special_channel_file = "TR_0_0.dat";
  bool enable_special_override = true;
  int special_channel_index = 3;
  double tsample_ns = 0.2;
  int pedestal_window = 100;
  double ped_target = 3500.0;
  std::string event_policy = "error";
};

inline bool LoadConfigFromJson(const std::string &path,
                               WaveConverterConfig &cfg,
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
    if (GetNumber(common, "n_channels", numValue)) {
      cfg.common.n_channels = static_cast<int>(numValue);
    }
    if (GetNumber(common, "max_cores", numValue)) {
      cfg.common.max_cores = static_cast<int>(numValue);
    }
    if (GetNumber(common, "max_events", numValue)) {
      cfg.common.max_events = static_cast<int>(numValue);
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
  }

  simdjson::dom::element waveformConverter;
  if (GetObject(root, "waveform_converter", waveformConverter)) {
    std::string strValue;
    double numValue = 0.0;
    bool boolValue = false;

    if (GetString(waveformConverter, "input_pattern", strValue)) {
      cfg.input_pattern = strValue;
    }
    if (GetString(waveformConverter, "input_dir", strValue)) {
      cfg.input_dir = strValue;
    }
    if (GetBool(waveformConverter, "input_is_ascii", boolValue)) {
      cfg.input_is_ascii = boolValue;
    }
    if (GetString(waveformConverter, "special_channel_file", strValue)) {
      cfg.special_channel_file = strValue;
    }
    if (GetBool(waveformConverter, "enable_special_override", boolValue)) {
      cfg.enable_special_override = boolValue;
    }
    if (GetNumber(waveformConverter, "special_channel_index", numValue)) {
      cfg.special_channel_index = static_cast<int>(numValue);
    }
    if (GetString(waveformConverter, "event_policy", strValue)) {
      cfg.event_policy = strValue;
    }
    if (GetNumber(waveformConverter, "tsample_ns", numValue)) {
      cfg.tsample_ns = numValue;
    }
    if (GetNumber(waveformConverter, "pedestal_window", numValue)) {
      cfg.pedestal_window = static_cast<int>(numValue);
    }
    if (GetNumber(waveformConverter, "ped_target", numValue)) {
      cfg.ped_target = numValue;
    }
  }

  return true;
}
