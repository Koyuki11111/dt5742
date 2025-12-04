#include "analysis/waveform_math.h"

#include <algorithm>
#include <cmath>

namespace {

float Interpolate(float x1, float y1, float x2, float y2, float yTarget) {
  if (std::abs(y2 - y1) < 1e-9f) {
    return x1;
  }
  return x1 + (x2 - x1) / (y2 - y1) * (yTarget - y1);
}

int FindTimeIndex(const std::vector<float> &time, float threshold) {
  for (size_t i = 0; i < time.size(); ++i) {
    if (time[i] >= threshold) {
      return static_cast<int>(i);
    }
  }
  return static_cast<int>(time.size()) - 1;
}

} // namespace

WindowIndices BuildWindowIndices(const std::vector<float> &time,
                                 float region_min,
                                 float region_max,
                                 int analysis_start,
                                 int analysis_end) {
  WindowIndices window;
  if (time.empty()) {
    return window;
  }

  int start_idx = FindTimeIndex(time, region_min);
  int end_idx = FindTimeIndex(time, region_max);

  int nSamples = static_cast<int>(time.size());
  window.start = std::max(0, std::min(std::max(start_idx, analysis_start), nSamples - 1));
  window.end = std::max(window.start, std::min(std::min(end_idx, analysis_end), nSamples - 1));
  return window;
}

BaselineNoiseMetrics ComputeBaselineAndNoise(const std::vector<float> &amp,
                                             const WindowIndices &baseline_window) {
  BaselineNoiseMetrics metrics;
  metrics.amp_min = 100000.0f;
  metrics.amp_max = -100000.0f;

  const int nSamples = static_cast<int>(amp.size());
  if (nSamples == 0) {
    return metrics;
  }

  int start = std::max(0, baseline_window.start);
  int end = std::min(baseline_window.end, nSamples - 1);

  int npointbaseline = 0;
  for (int i = start; i <= end; ++i) {
    float val = amp[i];
    metrics.baseline += val;
    metrics.amp_min = std::min(metrics.amp_min, val);
    metrics.amp_max = std::max(metrics.amp_max, val);
    ++npointbaseline;
  }

  if (npointbaseline > 0) {
    metrics.baseline /= static_cast<float>(npointbaseline);
  }

  int npointnoise = 0;
  int npoint1noise = 0;

  for (int i = start; i <= end; ++i) {
    float val = amp[i];
    metrics.rms_noise += (val - metrics.baseline) * (val - metrics.baseline);
    ++npointnoise;

    float val1 = 0.0f;
    int count1 = 0;
    if (i > 0) {
      val1 += amp[i-1];
      ++count1;
    }
    val1 += val;
    ++count1;
    if (i + 1 < nSamples) {
      val1 += amp[i+1];
      ++count1;
    }
    metrics.noise1_point += (val1 / static_cast<float>(count1)) - metrics.baseline;
    ++npoint1noise;
  }

  if (npointnoise > 0) {
    metrics.rms_noise = std::sqrt(metrics.rms_noise / static_cast<float>(npointnoise));
  }
  if (npoint1noise > 0) {
    metrics.noise1_point /= static_cast<float>(npoint1noise);
  }

  return metrics;
}

std::vector<float> ApplyBaselineAndPolarity(const std::vector<float> &amp,
                                            float baseline,
                                            int polarity) {
  std::vector<float> amp_corr(amp.size());
  for (size_t i = 0; i < amp.size(); ++i) {
    amp_corr[i] = (amp[i] - baseline) * static_cast<float>(polarity);
  }
  return amp_corr;
}

PeakMetrics FindPeakInWindow(const std::vector<float> &amp_corr,
                             const std::vector<float> &time,
                             const WindowIndices &signal_window) {
  PeakMetrics peak;
  if (amp_corr.empty() || time.size() != amp_corr.size()) {
    return peak;
  }

  int nSamples = static_cast<int>(amp_corr.size());
  int start = std::max(0, signal_window.start);
  int end = std::min(signal_window.end, nSamples - 1);

  peak.index = start;
  for (int i = start; i <= end; ++i) {
    float val = amp_corr[i];
    if (val > peak.amplitude) {
      peak.amplitude = val;
      peak.index = i;
    }
  }

  peak.time = time[peak.index];
  return peak;
}

float IntegrateChargeWindow(const std::vector<float> &amp_corr,
                            const WindowIndices &charge_window,
                            float dt,
                            float impedance) {
  float charge = 0.0f;
  int nSamples = static_cast<int>(amp_corr.size());
  if (nSamples == 0) {
    return charge;
  }

  int start = std::max(0, charge_window.start);
  int end = std::min(charge_window.end, nSamples - 1);
  for (int i = start; i < end && i < nSamples; ++i) {
    charge += amp_corr[i] * dt / impedance;
  }

  return charge;
}

std::vector<float> ComputeChargeFractionTimes(const std::vector<float> &amp_corr,
                                              const std::vector<float> &time,
                                              const WindowIndices &charge_window,
                                              float dt,
                                              float impedance,
                                              float total_charge,
                                              const std::vector<int> &thresholds_percent,
                                              float charge_min,
                                              float charge_max) {
  const size_t nChTh = thresholds_percent.size();
  std::vector<float> time_charge(nChTh, 10.0f);
  if (amp_corr.empty() || time.size() != amp_corr.size()) {
    return time_charge;
  }

  int nSamples = static_cast<int>(amp_corr.size());
  int start = std::max(0, charge_window.start);
  int end = std::min(charge_window.end, nSamples - 1);

  std::vector<float> charge_thresholds(nChTh);
  for (size_t i = 0; i < nChTh; ++i) {
    charge_thresholds[i] = total_charge * thresholds_percent[i] / 100.0f;
  }

  float chargetemp = 0.0f;
  size_t chth = 0;
  for (int i = start; i < end && i < nSamples; ++i) {
    chargetemp += amp_corr[i] * dt / impedance;
    if (chth < nChTh && chargetemp > charge_thresholds[chth]) {
      if (i > 0) {
        float prevCharge = chargetemp - amp_corr[i] * dt / impedance;
        float timechargetemp = Interpolate(time[i-1], prevCharge,
                                           time[i], chargetemp,
                                           charge_thresholds[chth]);
        if (timechargetemp >= charge_min && timechargetemp <= charge_max) {
          time_charge[chth] = timechargetemp;
        }
      }
      chth++;
      if (chth == nChTh) {
        break;
      }
    }
  }

  return time_charge;
}

ThresholdCrossing FindThresholdCrossingBackward(const std::vector<float> &amp_corr,
                                                const std::vector<float> &time,
                                                int from_idx,
                                                int stop_idx,
                                                float threshold,
                                                float rms_noise) {
  ThresholdCrossing crossing;
  if (amp_corr.empty() || time.size() != amp_corr.size()) {
    return crossing;
  }

  int nSamples = static_cast<int>(amp_corr.size());
  int start = std::min(from_idx, nSamples - 1);
  int stop = std::max(0, stop_idx);

  for (int i = start; i > stop; --i) {
    if (amp_corr[i] < threshold) {
      if (i + 1 < nSamples) {
        crossing.time = Interpolate(time[i], amp_corr[i],
                                    time[i+1], amp_corr[i+1],
                                    threshold);
        float slew = (amp_corr[i+1] - amp_corr[i]) / (time[i+1] - time[i]);
        if (std::abs(slew) > 1e-9f) {
          crossing.jitter = rms_noise / std::abs(slew);
        }
        crossing.found = true;
      }
      break;
    }
  }

  return crossing;
}

ThresholdCrossing FindTrailingEdgeForward(const std::vector<float> &amp_corr,
                                          const std::vector<float> &time,
                                          int from_idx,
                                          int stop_idx,
                                          float threshold) {
  ThresholdCrossing crossing;
  if (amp_corr.empty() || time.size() != amp_corr.size()) {
    return crossing;
  }

  int nSamples = static_cast<int>(amp_corr.size());
  int start = std::max(1, from_idx);
  int stop = std::min(stop_idx, nSamples - 1);

  for (int i = start; i < stop; ++i) {
    if (amp_corr[i] < threshold) {
      float trailing_time = Interpolate(time[i-1], amp_corr[i-1],
                                        time[i], amp_corr[i],
                                        threshold);
      crossing.time = trailing_time;
      crossing.found = true;
      break;
    }
  }

  return crossing;
}

WaveformFeatures AnalyzeWaveform(const std::vector<float> &amp,
                                 const std::vector<float> &time,
                                 const AnalysisConfig &cfg,
                                 int channel) {
  WaveformFeatures features;
  const int nSamples = static_cast<int>(amp.size());

  if (nSamples == 0 || time.size() != amp.size()) {
    return features;
  }

  const float analysisMin = cfg.analysis_region_min[channel];
  const float analysisMax = cfg.analysis_region_max[channel];
  const float baselineMin = cfg.baseline_region_min[channel];
  const float baselineMax = cfg.baseline_region_max[channel];
  const float signalMin = cfg.signal_region_min[channel];
  const float signalMax = cfg.signal_region_max[channel];
  const float chargeMin = cfg.charge_region_min[channel];
  const float chargeMax = cfg.charge_region_max[channel];
  const int polarity = cfg.signal_polarity[channel];

  float dT = (nSamples > 1) ? (time[1] - time[0]) : 0.2f;

  WindowIndices analysis_window;
  analysis_window.start = 0;
  analysis_window.end = nSamples - 1;
  for (int i = 0; i < nSamples; ++i) {
    if (time[i] >= analysisMin) {
      analysis_window.start = i;
      break;
    }
  }
  for (int i = nSamples - 1; i >= 0; --i) {
    if (time[i] <= analysisMax) {
      analysis_window.end = i;
      break;
    }
  }

  WindowIndices baseline_window = BuildWindowIndices(time, baselineMin, baselineMax,
                                                     analysis_window.start, analysis_window.end);
  BaselineNoiseMetrics baseline_metrics = ComputeBaselineAndNoise(amp, baseline_window);
  features.baseline = baseline_metrics.baseline;
  features.rmsNoise = baseline_metrics.rms_noise;
  features.noise1Point = baseline_metrics.noise1_point;
  features.ampMinBefore = baseline_metrics.amp_min;
  features.ampMaxBefore = baseline_metrics.amp_max;

  std::vector<float> ampCorr = ApplyBaselineAndPolarity(amp, features.baseline, polarity);

  WindowIndices signal_window = BuildWindowIndices(time, signalMin, signalMax,
                                                   analysis_window.start, analysis_window.end);
  PeakMetrics peak = FindPeakInWindow(ampCorr, time, signal_window);
  features.ampMax = peak.amplitude;
  features.peakTime = peak.time;
  int posampmax = peak.index;

  if (features.rmsNoise > 0.0f) {
    features.signalOverNoise = features.ampMax / features.rmsNoise;
    features.hasSignal = (features.signalOverNoise >= cfg.snr_threshold) &&
                         (features.ampMax >= cfg.cut_amp_max[channel]);
  }

  WindowIndices charge_window = BuildWindowIndices(time, chargeMin, chargeMax,
                                                   analysis_window.start, analysis_window.end);
  features.charge = IntegrateChargeWindow(ampCorr, charge_window, dT, cfg.impedance);
  features.timeCharge = ComputeChargeFractionTimes(ampCorr, time, charge_window, dT,
                                                   cfg.impedance, features.charge,
                                                   cfg.charge_thresholds, chargeMin, chargeMax);

  const size_t nCFD = cfg.cfd_thresholds.size();
  features.timeCFD.assign(nCFD, 0.0f);
  features.jitterCFD.assign(nCFD, 0.0f);

  for (int b = static_cast<int>(nCFD) - 1; b >= 0; --b) {
    float threshold = features.ampMax * (cfg.cfd_thresholds[b] / 100.0f);
    ThresholdCrossing crossing = FindThresholdCrossingBackward(ampCorr, time,
                                                               posampmax, signal_window.start,
                                                               threshold, features.rmsNoise);
    if (crossing.found) {
      features.timeCFD[b] = crossing.time;
      features.jitterCFD[b] = crossing.jitter;
    }
  }

  const size_t nLE = cfg.le_thresholds.size();
  features.timeLE.assign(nLE, 20.0f);
  features.jitterLE.assign(nLE, -5.0f);
  features.totLE.assign(nLE, -5.0f);

  for (int b = static_cast<int>(nLE) - 1; b >= 0; --b) {
    float threshold = cfg.le_thresholds[b] / 1000.0f;
    if (features.ampMax <= threshold) {
      continue;
    }

    ThresholdCrossing leading = FindThresholdCrossingBackward(ampCorr, time,
                                                              posampmax, signal_window.start,
                                                              threshold, features.rmsNoise);
    if (leading.found) {
      features.timeLE[b] = leading.time;
      features.jitterLE[b] = leading.jitter;

      ThresholdCrossing trailing = FindTrailingEdgeForward(ampCorr, time,
                                                           posampmax, charge_window.end,
                                                           threshold);
      if (trailing.found) {
        features.totLE[b] = trailing.time - features.timeLE[b];
      }
    }
  }

  float amp90 = features.ampMax * cfg.rise_time_high;
  float amp10 = features.ampMax * cfg.rise_time_low;

  ThresholdCrossing crossing90 = FindThresholdCrossingBackward(ampCorr, time,
                                                               posampmax, signal_window.start,
                                                               amp90, features.rmsNoise);
  ThresholdCrossing crossing10 = FindThresholdCrossingBackward(ampCorr, time,
                                                               posampmax, signal_window.start,
                                                               amp10, features.rmsNoise);

  float time90 = crossing90.found ? crossing90.time : 0.0f;
  float time10 = crossing10.found ? crossing10.time : 0.0f;

  features.riseTime = time90 - time10;
  if (features.riseTime > 0.0f) {
    features.slewRate = (amp90 - amp10) / features.riseTime;
  }

  return features;
}
