#pragma once

#include <string>

// Fields shared across stages.
struct CommonConfig {
  std::string daq_name   = "daq01";
  std::string input_dir  = "output";
  std::string output_dir = "output";
  int runnumber  = 1;
  int n_channels = 16;
  int n_sensors  = 1;
  int max_cores  = 8;
  int chunk_size = 100;
  int max_events = 100;
  std::string temp_dir = "./temp";
  std::string waveforms_root = "waveforms.root";
  std::string waveforms_tree = "Waveforms";
  std::string analysis_root = "waveforms_analyzed.root";
  std::string analysis_tree = "Analysis";
  // Handling for events where channels have different sample counts.
  // Supported values: "strict" (fail on mismatch), "pad" (pad to max samples).
  std::string nsamples_policy = "strict";
};
