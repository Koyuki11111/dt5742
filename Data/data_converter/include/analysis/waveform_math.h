#pragma once

#include <vector>

#include "config/analysis_config.h"

struct WindowIndices {
  int start = 0;
  int end = 0;
};

struct BaselineNoiseMetrics {
  float baseline = 0.0f;
  float rms_noise = 0.0f;
  float noise1_point = 0.0f;
  float amp_min = 0.0f;
  float amp_max = 0.0f;
};

struct PeakMetrics {
  float amplitude = 0.0f;
  int index = 0;
  float time = 0.0f;
};

struct ThresholdCrossing {
  float time = 0.0f;
  float jitter = 0.0f;
  bool found = false;
};

struct WaveformFeatures {
  // Baseline and noise
  float baseline = 0.0f;
  float rmsNoise = 0.0f;
  float noise1Point = 0.0f;
  float ampMinBefore = 0.0f;
  float ampMaxBefore = 0.0f;

  // Signal characteristics
  bool hasSignal = false;
  float ampMax = 0.0f;
  float charge = 0.0f;
  float signalOverNoise = 0.0f;
  float peakTime = 0.0f;

  // Timing
  float riseTime = 0.0f;
  float slewRate = 0.0f;

  // Multi-threshold timing
  std::vector<float> timeCFD;
  std::vector<float> jitterCFD;
  std::vector<float> timeLE;
  std::vector<float> jitterLE;
  std::vector<float> totLE;
  std::vector<float> timeCharge;
};

WindowIndices BuildWindowIndices(const std::vector<float> &time,
                                 float region_min,
                                 float region_max,
                                 int analysis_start,
                                 int analysis_end);

BaselineNoiseMetrics ComputeBaselineAndNoise(const std::vector<float> &amp,
                                             const WindowIndices &baseline_window);

std::vector<float> ApplyBaselineAndPolarity(const std::vector<float> &amp,
                                            float baseline,
                                            int polarity);

PeakMetrics FindPeakInWindow(const std::vector<float> &amp_corr,
                             const std::vector<float> &time,
                             const WindowIndices &signal_window);

float IntegrateChargeWindow(const std::vector<float> &amp_corr,
                            const WindowIndices &charge_window,
                            float dt,
                            float impedance);

std::vector<float> ComputeChargeFractionTimes(const std::vector<float> &amp_corr,
                                              const std::vector<float> &time,
                                              const WindowIndices &charge_window,
                                              float dt,
                                              float impedance,
                                              float total_charge,
                                              const std::vector<int> &thresholds_percent,
                                              float charge_min,
                                              float charge_max);

ThresholdCrossing FindThresholdCrossingBackward(const std::vector<float> &amp_corr,
                                                const std::vector<float> &time,
                                                int from_idx,
                                                int stop_idx,
                                                float threshold,
                                                float rms_noise);

ThresholdCrossing FindTrailingEdgeForward(const std::vector<float> &amp_corr,
                                          const std::vector<float> &time,
                                          int from_idx,
                                          int stop_idx,
                                          float threshold);

WaveformFeatures AnalyzeWaveform(const std::vector<float> &amp,
                                 const std::vector<float> &time,
                                 const AnalysisConfig &cfg,
                                 int channel);
