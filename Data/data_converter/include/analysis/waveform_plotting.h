#pragma once

#include <vector>

#include "config/analysis_config.h"
#include "analysis/waveform_math.h"

class TFile;

void SaveWaveformPlots(TFile *waveformPlotsFile, int event, int channel,
                       const std::vector<float> &amp,
                       const std::vector<float> &time,
                       const WaveformFeatures &features,
                       const AnalysisConfig &cfg);
