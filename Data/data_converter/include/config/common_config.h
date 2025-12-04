#pragma once

#include <string>

// Fields shared across stages.
struct CommonConfig {
  std::string output_dir = "output";
  int n_channels = 16;
  int max_cores = 8;
  int chunk_size = 100;
  std::string temp_dir = "./temp";
  std::string waveforms_root = "waveforms.root";
  std::string waveforms_tree = "Waveforms";
  std::string analysis_root = "waveforms_analyzed.root";
  std::string analysis_tree = "Analysis";
  // Handling for events where channels have different sample counts.
  // Supported values: "strict" (fail on mismatch), "pad" (pad to max samples).
  std::string nsamples_policy = "strict";
};
