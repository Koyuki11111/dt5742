#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TH2F.h"

#include "analysis_config.h"

namespace {

// Interpolate to find exact crossing point
float Interpolate(float x1, float y1, float x2, float y2, float yTarget) {
  if (std::abs(y2 - y1) < 1e-9f) {
    return x1;
  }
  return x1 + (x2 - x1) / (y2 - y1) * (yTarget - y1);
}

// Find index where time exceeds threshold
int FindTimeIndex(const std::vector<float> &time, float threshold) {
  for (size_t i = 0; i < time.size(); ++i) {
    if (time[i] >= threshold) {
      return static_cast<int>(i);
    }
  }
  return static_cast<int>(time.size()) - 1;
}

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

// Save waveform plots for a single waveform
void SaveWaveformPlots(TFile *waveformPlotsFile, int event, int channel,
                      const std::vector<float> &amp,
                      const std::vector<float> &time,
                      const WaveformFeatures &features,
                      const AnalysisConfig &cfg) {
  if (!waveformPlotsFile || waveformPlotsFile->IsZombie()) {
    return;
  }

  // Create or navigate to event directory inside ROOT file
  char eventDirName[64];
  std::snprintf(eventDirName, sizeof(eventDirName), "event_%06d", event);

  TDirectory *eventDir = waveformPlotsFile->GetDirectory(eventDirName);
  if (!eventDir) {
    eventDir = waveformPlotsFile->mkdir(eventDirName);
  }
  if (!eventDir) {
    std::cerr << "WARNING: Failed to create TDirectory " << eventDirName << std::endl;
    return;
  }

  // Create or navigate to sensor subdirectory within event directory
  int sensorID = cfg.sensor_ids[channel];
  char sensorDirName[128];
  std::snprintf(sensorDirName, sizeof(sensorDirName), "%s/sensor%02d", eventDirName, sensorID);

  TDirectory *sensorDir = waveformPlotsFile->GetDirectory(sensorDirName);
  if (!sensorDir) {
    sensorDir = eventDir->mkdir(Form("sensor%02d", sensorID));
  }
  if (!sensorDir) {
    std::cerr << "WARNING: Failed to create TDirectory " << sensorDirName << std::endl;
    return;
  }
  sensorDir->cd();

  const int nSamples = static_cast<int>(amp.size());
  const int polarity = cfg.signal_polarity[channel];
  const int stripID = cfg.strip_ids[channel];

  // Find analysis region indices
  const float analysisMin = cfg.analysis_region_min[channel];
  const float analysisMax = cfg.analysis_region_max[channel];
  int analysisStartIdx = 0;
  int analysisEndIdx = nSamples - 1;

  for (int i = 0; i < nSamples; ++i) {
    if (time[i] >= analysisMin) {
      analysisStartIdx = i;
      break;
    }
  }
  for (int i = nSamples - 1; i >= 0; --i) {
    if (time[i] <= analysisMax) {
      analysisEndIdx = i;
      break;
    }
  }

  // Calculate number of points in analysis region
  const int nAnalysisPoints = analysisEndIdx - analysisStartIdx + 1;

  // ===== 1. Save raw waveform as TGraph (baseline-subtracted, no polarity flip) =====
  // Only save points within analysis region
  char rawGraphName[64];
  std::snprintf(rawGraphName, sizeof(rawGraphName), "strip%02d_raw", stripID);

  TGraph *grRaw = new TGraph(nAnalysisPoints);
  grRaw->SetName(rawGraphName);
  grRaw->SetTitle(Form("Event %d, Sensor %d, Strip %d (Ch%d) - Raw Waveform (Baseline Subtracted);Time (ns);Amplitude (V)", event, sensorID, stripID, channel));
  int pointIdx = 0;
  for (int i = analysisStartIdx; i <= analysisEndIdx; ++i) {
    grRaw->SetPoint(pointIdx++, time[i], amp[i] - features.baseline);
  }
  grRaw->Write(rawGraphName, TObject::kOverwrite);

  // ===== 2. Save 3-point moving average as TGraph =====
  // Only save points within analysis region
  char avgGraphName[64];
  std::snprintf(avgGraphName, sizeof(avgGraphName), "strip%02d_avg", stripID);

  TGraph *grAvg = new TGraph(nAnalysisPoints);
  grAvg->SetName(avgGraphName);
  grAvg->SetTitle(Form("Event %d, Sensor %d, Strip %d (Ch%d) - 3-Point Moving Average;Time (ns);Amplitude (V)", event, sensorID, stripID, channel));

  pointIdx = 0;
  for (int i = analysisStartIdx; i <= analysisEndIdx; ++i) {
    float avgAmp = 0.0f;
    if (i == 0 || i == analysisStartIdx) {
      // First point in overall or analysis region: average of points i, i+1
      if (i + 1 < nSamples) {
        avgAmp = ((amp[i] - features.baseline) + (amp[i+1] - features.baseline)) / 2.0f;
      } else {
        avgAmp = (amp[i] - features.baseline);
      }
    } else if (i == nSamples - 1 || i == analysisEndIdx) {
      // Last point in overall or analysis region: average of points i-1, i
      avgAmp = ((amp[i-1] - features.baseline) + (amp[i] - features.baseline)) / 2.0f;
    } else {
      // Middle points: 3-point average
      avgAmp = ((amp[i-1] - features.baseline) + (amp[i] - features.baseline) + (amp[i+1] - features.baseline)) / 3.0f;
    }
    grAvg->SetPoint(pointIdx++, time[i], avgAmp);
  }
  grAvg->Write(avgGraphName, TObject::kOverwrite);

  // ===== 3. Save analysis canvas with baseline and CFD markers =====
  // Only show points within analysis region
  char canvasName[64];
  std::snprintf(canvasName, sizeof(canvasName), "strip%02d_analysis", stripID);

  TCanvas *canvas = new TCanvas(canvasName,
                                Form("Event %d, Sensor %d, Strip %d (Ch%d) - Analysis", event, sensorID, stripID, channel),
                                1200, 800);
  canvas->SetGrid();

  // Create baseline-subtracted waveform (NO polarity flip for display)
  // Only include points in analysis region
  TGraph *grAnalysis = new TGraph(nAnalysisPoints);
  const char *polarityStr = (polarity > 0) ? "Positive" : "Negative";
  grAnalysis->SetTitle(Form("Event %d, Sensor %d, Strip %d (Ch%d) - Analysis (%s Signal);Time (ns);Amplitude (V)",
                            event, sensorID, stripID, channel, polarityStr));
  pointIdx = 0;
  for (int i = analysisStartIdx; i <= analysisEndIdx; ++i) {
    float ampSub = amp[i] - features.baseline;
    grAnalysis->SetPoint(pointIdx++, time[i], ampSub);
  }
  grAnalysis->SetLineColor(kBlue);
  grAnalysis->SetLineWidth(2);
  grAnalysis->Draw("AL");

  // Draw baseline (at y=0 after subtraction)
  // Use analysis region time range
  TLine *baselineLine = new TLine(time[analysisStartIdx], 0.0, time[analysisEndIdx], 0.0);
  baselineLine->SetLineColor(kGray+2);
  baselineLine->SetLineStyle(2);
  baselineLine->SetLineWidth(2);
  baselineLine->Draw("same");

  // Create legend
  TLegend *legend = new TLegend(0.65, 0.65, 0.89, 0.89);
  legend->SetTextSize(0.025);
  legend->AddEntry(grAnalysis, "Waveform", "l");
  legend->AddEntry(baselineLine, Form("Baseline = %.3f V", features.baseline), "l");

  // Draw CFD timing markers (with correct polarity consideration)
  const size_t nCFD = cfg.cfd_thresholds.size();
  int colorsCFD[] = {kRed, kGreen+2, kOrange, kViolet};

  for (size_t i = 0; i < nCFD && i < features.timeCFD.size(); ++i) {
    // Threshold is in analysis space (after polarity flip), so flip back for display
    float thresholdDisplay = features.ampMax * (cfg.cfd_thresholds[i] / 100.0f) * polarity;
    float timeCFD = features.timeCFD[i];

    if (timeCFD > 0.0f) {
      // Draw horizontal line at threshold (only in analysis region)
      TLine *threshLine = new TLine(time[analysisStartIdx], thresholdDisplay, time[analysisEndIdx], thresholdDisplay);
      threshLine->SetLineColor(colorsCFD[i % 4]);
      threshLine->SetLineStyle(3);
      threshLine->Draw("same");

      // Draw marker at CFD crossing point
      TMarker *marker = new TMarker(timeCFD, thresholdDisplay, 20);
      marker->SetMarkerColor(colorsCFD[i % 4]);
      marker->SetMarkerSize(1.5);
      marker->Draw("same");

      legend->AddEntry(marker, Form("CFD %d%% @ %.2f ns", cfg.cfd_thresholds[i], timeCFD), "p");
    }
  }

  // Add analysis info text
  TLatex *infoText = new TLatex();
  infoText->SetNDC();
  infoText->SetTextSize(0.025);
  infoText->DrawLatex(0.15, 0.85, Form("Peak Amp: %.3f V (abs)", features.ampMax));
  infoText->DrawLatex(0.15, 0.82, Form("Peak Time: %.2f ns", features.peakTime));
  infoText->DrawLatex(0.15, 0.79, Form("Rise Time: %.2f ns", features.riseTime));
  infoText->DrawLatex(0.15, 0.76, Form("Charge: %.3f pC", features.charge * 1e12));
  infoText->DrawLatex(0.15, 0.73, Form("SNR: %.1f", features.signalOverNoise));
  infoText->DrawLatex(0.15, 0.70, Form("RMS Noise: %.4f V", features.rmsNoise));
  infoText->DrawLatex(0.15, 0.67, Form("Polarity: %s", polarityStr));

  legend->Draw();
  canvas->Write(canvasName, TObject::kOverwrite);

  // Clean up
  delete grRaw;
  delete grAvg;
  delete grAnalysis;
  delete baselineLine;
  delete legend;
  delete infoText;
  delete canvas;

  // Return to main directory
  waveformPlotsFile->cd();
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

  // Calculate time step
  float dT = (nSamples > 1) ? (time[1] - time[0]) : 0.2f;

  // Find valid analysis region indices
  int analysisStartIdx = 0;
  int analysisEndIdx = nSamples - 1;
  for (int i = 0; i < nSamples; ++i) {
    if (time[i] >= analysisMin) {
      analysisStartIdx = i;
      break;
    }
  }
  for (int i = nSamples - 1; i >= 0; --i) {
    if (time[i] <= analysisMax) {
      analysisEndIdx = i;
      break;
    }
  }

  // ========== Step 1: Calculate baseline and noise in configured region ==========
  int npointbaseline = 0;
  features.baseline = 0.0f;
  features.ampMinBefore = 100000.0f;
  features.ampMaxBefore = -100000.0f;

  int posBaselineStart = FindTimeIndex(time, baselineMin);
  int posBaselineEnd = FindTimeIndex(time, baselineMax);

  // Apply analysis region limits to baseline region
  posBaselineStart = std::max(posBaselineStart, analysisStartIdx);
  posBaselineEnd = std::min(posBaselineEnd, analysisEndIdx);

  // Calculate mean baseline in the specified region
  for (int i = posBaselineStart; i <= posBaselineEnd && i < nSamples; ++i) {
    npointbaseline++;
    features.baseline += amp[i];
    if (amp[i] > features.ampMaxBefore) {
      features.ampMaxBefore = amp[i];
    }
    if (amp[i] < features.ampMinBefore) {
      features.ampMinBefore = amp[i];
    }
  }

  if (npointbaseline > 0) {
    features.baseline /= static_cast<float>(npointbaseline);
  }

  features.ampMaxBefore -= features.baseline;
  features.ampMinBefore -= features.baseline;
  features.noise1Point = (nSamples > 10) ? (amp[10] - features.baseline) : 0.0f;

  // Calculate RMS noise in baseline region
  features.rmsNoise = 0.0f;
  for (int i = posBaselineStart; i <= posBaselineEnd && i < nSamples; ++i) {
    float diff = amp[i] - features.baseline;
    features.rmsNoise += diff * diff;
  }
  if (npointbaseline > 1) {
    features.rmsNoise = std::sqrt(features.rmsNoise / (npointbaseline - 1));
  }

  // Create baseline-subtracted waveform (apply polarity)
  std::vector<float> ampSub(nSamples);
  for (int i = 0; i < nSamples; ++i) {
    ampSub[i] = (amp[i] - features.baseline) * polarity;
  }

  // ========== Step 2: Find peak amplitude (considering polarity) ==========
  int posampmax = 0;
  features.ampMax = -100000.0f;

  int posSignalStart = FindTimeIndex(time, signalMin);
  int posSignalEnd = FindTimeIndex(time, signalMax);

  // Apply analysis region limits to signal search
  posSignalStart = std::max(posSignalStart, analysisStartIdx);
  posSignalEnd = std::min(posSignalEnd, analysisEndIdx);

  for (int i = posSignalStart; i < posSignalEnd && i < nSamples; ++i) {
    if (ampSub[i] > features.ampMax) {
      features.ampMax = ampSub[i];
      posampmax = i;
    }
  }

  features.peakTime = time[posampmax];
  features.signalOverNoise = (features.rmsNoise > 0.0f)
                              ? (features.ampMax / features.rmsNoise)
                              : 0.0f;

  // ========== Signal detection: check if SNR exceeds threshold ==========
  features.hasSignal = (features.signalOverNoise > cfg.snr_threshold);

  // ========== Step 3: Calculate charge ==========
  features.charge = 0.0f;
  int posChargeStart = FindTimeIndex(time, chargeMin);
  int posChargeEnd = FindTimeIndex(time, chargeMax);

  // Apply analysis region limits to charge region
  posChargeStart = std::max(posChargeStart, analysisStartIdx);
  posChargeEnd = std::min(posChargeEnd, analysisEndIdx);

  for (int i = posChargeStart; i < posChargeEnd && i < nSamples; ++i) {
    features.charge += ampSub[i];
  }
  features.charge *= dT;
  features.charge /= cfg.impedance;

  // ========== Step 4: Calculate charge-based timing ==========
  const size_t nChTh = cfg.charge_thresholds.size();
  features.timeCharge.assign(nChTh, 10.0f);
  std::vector<float> chargeThresholds(nChTh);

  for (size_t i = 0; i < nChTh; ++i) {
    chargeThresholds[i] = features.charge * cfg.charge_thresholds[i] / 100.0f;
  }

  float chargetemp = 0.0f;
  size_t chth = 0;
  for (int i = posChargeStart; i < posChargeEnd && i < nSamples; ++i) {
    chargetemp += ampSub[i] * dT / cfg.impedance;
    if (chth < nChTh && chargetemp > chargeThresholds[chth]) {
      if (i > 0) {
        float prevCharge = chargetemp - ampSub[i] * dT / cfg.impedance;
        float timechargetemp = Interpolate(time[i-1], prevCharge,
                                          time[i], chargetemp,
                                          chargeThresholds[chth]);
        if (timechargetemp >= chargeMin && timechargetemp <= chargeMax) {
          features.timeCharge[chth] = timechargetemp;
        }
      }
      chth++;
      if (chth == nChTh) {
        break;
      }
    }
  }

  // ========== Step 5: CFD timing (from peak backwards) ==========
  const size_t nCFD = cfg.cfd_thresholds.size();
  features.timeCFD.assign(nCFD, 0.0f);
  features.jitterCFD.assign(nCFD, 0.0f);

  for (int b = static_cast<int>(nCFD) - 1; b >= 0; --b) {
    float threshold = features.ampMax * (cfg.cfd_thresholds[b] / 100.0f);
    for (int i = posampmax; i > posSignalStart; --i) {
      if (ampSub[i] < threshold) {
        if (i + 1 < nSamples) {
          features.timeCFD[b] = Interpolate(time[i], ampSub[i],
                                           time[i+1], ampSub[i+1],
                                           threshold);
          float slewCFD = (ampSub[i+1] - ampSub[i]) / (time[i+1] - time[i]);
          if (std::abs(slewCFD) > 1e-9f) {
            features.jitterCFD[b] = features.rmsNoise / std::abs(slewCFD);
          }
        }
        break;
      }
    }
  }

  // ========== Step 6: Leading Edge timing ==========
  const size_t nLE = cfg.le_thresholds.size();
  features.timeLE.assign(nLE, 20.0f);
  features.jitterLE.assign(nLE, -5.0f);

  for (int b = static_cast<int>(nLE) - 1; b >= 0; --b) {
    float threshold = cfg.le_thresholds[b] / 1000.0f; // mV to V
    if (features.ampMax <= threshold) {
      continue;
    }

    for (int i = posampmax; i > posSignalStart; --i) {
      if (ampSub[i] < threshold) {
        if (i + 1 < nSamples) {
          features.timeLE[b] = Interpolate(time[i], ampSub[i],
                                          time[i+1], ampSub[i+1],
                                          threshold);
          float slewLE = (ampSub[i+1] - ampSub[i]) / (time[i+1] - time[i]);
          if (std::abs(slewLE) > 1e-9f) {
            features.jitterLE[b] = features.rmsNoise / std::abs(slewLE);
          }
        }
        break;
      }
    }
  }

  // ========== Step 7: Time over Threshold (trailing edge) ==========
  features.totLE.assign(nLE, -5.0f);

  for (int b = static_cast<int>(nLE) - 1; b >= 0; --b) {
    float threshold = cfg.le_thresholds[b] / 1000.0f;
    if (features.ampMax <= threshold) {
      continue;
    }

    for (int i = posampmax; i < posChargeEnd && i < nSamples - 1; ++i) {
      if (ampSub[i] < threshold) {
        float trailingTime = Interpolate(time[i-1], ampSub[i-1],
                                        time[i], ampSub[i],
                                        threshold);
        features.totLE[b] = trailingTime - features.timeLE[b];
        break;
      }
    }
  }

  // ========== Step 8: Rise time (10% to 90%) ==========
  float amp90 = features.ampMax * cfg.rise_time_high;
  float amp10 = features.ampMax * cfg.rise_time_low;
  float time90 = 0.0f;
  float time10 = 0.0f;

  // Find 90% point
  for (int i = posampmax; i > posSignalStart; --i) {
    if (ampSub[i] < amp90) {
      if (i + 1 < nSamples) {
        time90 = Interpolate(time[i], ampSub[i],
                            time[i+1], ampSub[i+1],
                            amp90);
      }
      break;
    }
  }

  // Find 10% point
  for (int i = posampmax; i > posSignalStart; --i) {
    if (ampSub[i] < amp10) {
      if (i + 1 < nSamples) {
        time10 = Interpolate(time[i], ampSub[i],
                            time[i+1], ampSub[i+1],
                            amp10);
      }
      break;
    }
  }

  features.riseTime = time90 - time10;
  if (features.riseTime > 0.0f) {
    features.slewRate = (amp90 - amp10) / features.riseTime;
  }

  return features;
}

bool RunAnalysis(const AnalysisConfig &cfg, Long64_t eventStart = -1, Long64_t eventEnd = -1) {
  // Create waveform plots ROOT file if enabled
  TFile *waveformPlotsFile = nullptr;
  if (cfg.waveform_plots_enabled) {
    // Build waveform plots output path: output_dir/waveform_plots/waveform_plots_dir.root
    std::string waveformPlotsFileName = BuildOutputPath(cfg.output_dir, "waveform_plots", cfg.waveform_plots_dir + ".root");

    // Create directory if needed
    size_t lastSlash = waveformPlotsFileName.find_last_of('/');
    if (lastSlash != std::string::npos) {
      std::string dirPath = waveformPlotsFileName.substr(0, lastSlash);
      if (!CreateDirectoryIfNeeded(dirPath)) {
        std::cerr << "WARNING: Failed to create waveform plots output directory: " << dirPath << std::endl;
      }
    }

    waveformPlotsFile = TFile::Open(waveformPlotsFileName.c_str(), "RECREATE");
    if (!waveformPlotsFile || waveformPlotsFile->IsZombie()) {
      std::cerr << "WARNING: Failed to create waveform plots output file "
                << waveformPlotsFileName << std::endl;
      std::cerr << "         Continuing without waveform plots output..." << std::endl;
      waveformPlotsFile = nullptr;
    } else {
      std::cout << "Waveform plots output enabled. Saving to: " << waveformPlotsFileName << std::endl;
      if (cfg.waveform_plots_only_signal) {
        std::cout << "  Only saving waveforms with detected signals (SNR > "
                  << cfg.snr_threshold << ")" << std::endl;
      }
    }
  }

  // Build input path: output_dir/root/input_root
  std::string inputPath = BuildOutputPath(cfg.output_dir, "root", cfg.input_root);

  // Open input ROOT file
  TFile *inputFile = TFile::Open(inputPath.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open input ROOT file " << inputPath << std::endl;
    return false;
  }
  std::cout << "Reading input file: " << inputPath << std::endl;

  TTree *inputTree = dynamic_cast<TTree *>(inputFile->Get(cfg.input_tree.c_str()));
  if (!inputTree) {
    std::cerr << "ERROR: cannot find tree " << cfg.input_tree << std::endl;
    inputFile->Close();
    return false;
  }

  // Get number of entries
  Long64_t totalEntries = inputTree->GetEntries();
  if (totalEntries == 0) {
    std::cerr << "ERROR: input tree has no entries" << std::endl;
    inputFile->Close();
    return false;
  }

  // Determine event range to process
  Long64_t startEntry = (eventStart >= 0) ? eventStart : 0;
  Long64_t endEntry = (eventEnd >= 0) ? eventEnd : totalEntries;

  // Validate range
  if (startEntry < 0) startEntry = 0;
  if (endEntry > totalEntries) endEntry = totalEntries;
  if (startEntry >= endEntry) {
    std::cerr << "ERROR: invalid event range [" << startEntry << ", " << endEntry << ")" << std::endl;
    inputFile->Close();
    return false;
  }

  Long64_t nEntries = endEntry - startEntry;
  std::cout << "Processing event range [" << startEntry << ", " << endEntry << ") - "
            << nEntries << " events" << std::endl;

  // Set up input branches
  int eventIdx = 0;
  int nChannels = 0;
  int nsamples = 0;
  std::vector<float> *timeAxis = nullptr;
  std::vector<std::vector<float> *> chPed(cfg.n_channels, nullptr);

  inputTree->SetBranchAddress("event", &eventIdx);
  inputTree->SetBranchAddress("n_channels", &nChannels);
  inputTree->SetBranchAddress("nsamples", &nsamples);
  inputTree->SetBranchAddress("time_ns", &timeAxis);

  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    char bname[32];
    std::snprintf(bname, sizeof(bname), "ch%02d_ped", ch);
    if (inputTree->GetBranch(bname)) {
      inputTree->SetBranchAddress(bname, &chPed[ch]);
    }
  }

  // Build output path: output_dir/root/output_root
  std::string outputPath = BuildOutputPath(cfg.output_dir, "root", cfg.output_root);

  // Create directory if needed
  size_t lastSlash = outputPath.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = outputPath.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      inputFile->Close();
      return false;
    }
  }

  // Create output ROOT file
  TFile *outputFile = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "ERROR: cannot create output ROOT file " << outputPath << std::endl;
    inputFile->Close();
    return false;
  }
  std::cout << "Creating output file: " << outputPath << std::endl;

  TTree *outputTree = new TTree(cfg.output_tree.c_str(), "Analyzed waveform features");

  // Create output branches - per channel vectors
  int event = 0;
  std::vector<int> sensorID(cfg.n_channels);
  std::vector<int> stripID(cfg.n_channels);
  std::vector<bool> hasSignal(cfg.n_channels);
  std::vector<float> baseline(cfg.n_channels);
  std::vector<float> rmsNoise(cfg.n_channels);
  std::vector<float> noise1Point(cfg.n_channels);
  std::vector<float> ampMinBefore(cfg.n_channels);
  std::vector<float> ampMaxBefore(cfg.n_channels);
  std::vector<float> ampMax(cfg.n_channels);
  std::vector<float> charge(cfg.n_channels);
  std::vector<float> signalOverNoise(cfg.n_channels);
  std::vector<float> peakTime(cfg.n_channels);
  std::vector<float> riseTime(cfg.n_channels);
  std::vector<float> slewRate(cfg.n_channels);

  // Initialize sensor and strip IDs from config
  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    sensorID[ch] = cfg.sensor_ids[ch];
    stripID[ch] = cfg.strip_ids[ch];
  }

  outputTree->Branch("event", &event);
  outputTree->Branch("sensorID", &sensorID);
  outputTree->Branch("stripID", &stripID);
  outputTree->Branch("hasSignal", &hasSignal);
  outputTree->Branch("baseline", &baseline);
  outputTree->Branch("rmsNoise", &rmsNoise);
  outputTree->Branch("noise1Point", &noise1Point);
  outputTree->Branch("ampMinBefore", &ampMinBefore);
  outputTree->Branch("ampMaxBefore", &ampMaxBefore);
  outputTree->Branch("ampMax", &ampMax);
  outputTree->Branch("charge", &charge);
  outputTree->Branch("signalOverNoise", &signalOverNoise);
  outputTree->Branch("peakTime", &peakTime);
  outputTree->Branch("riseTime", &riseTime);
  outputTree->Branch("slewRate", &slewRate);

  // Multi-threshold timing branches (per channel, per threshold)
  const size_t nCFD = cfg.cfd_thresholds.size();
  const size_t nLE = cfg.le_thresholds.size();
  const size_t nCharge = cfg.charge_thresholds.size();

  std::vector<std::vector<float>> timeCFD(cfg.n_channels, std::vector<float>(nCFD));
  std::vector<std::vector<float>> jitterCFD(cfg.n_channels, std::vector<float>(nCFD));
  std::vector<std::vector<float>> timeLE(cfg.n_channels, std::vector<float>(nLE));
  std::vector<std::vector<float>> jitterLE(cfg.n_channels, std::vector<float>(nLE));
  std::vector<std::vector<float>> totLE(cfg.n_channels, std::vector<float>(nLE));
  std::vector<std::vector<float>> timeCharge(cfg.n_channels, std::vector<float>(nCharge));

  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    for (size_t i = 0; i < nCFD; ++i) {
      outputTree->Branch(Form("ch%02d_timeCFD_%dpc", ch, cfg.cfd_thresholds[i]),
                        &timeCFD[ch][i]);
      outputTree->Branch(Form("ch%02d_jitterCFD_%dpc", ch, cfg.cfd_thresholds[i]),
                        &jitterCFD[ch][i]);
    }
    for (size_t i = 0; i < nLE; ++i) {
      outputTree->Branch(Form("ch%02d_timeLE_%.1fmV", ch, cfg.le_thresholds[i]),
                        &timeLE[ch][i]);
      outputTree->Branch(Form("ch%02d_jitterLE_%.1fmV", ch, cfg.le_thresholds[i]),
                        &jitterLE[ch][i]);
      outputTree->Branch(Form("ch%02d_totLE_%.1fmV", ch, cfg.le_thresholds[i]),
                        &totLE[ch][i]);
    }
    for (size_t i = 0; i < nCharge; ++i) {
      outputTree->Branch(Form("ch%02d_timeCharge_%dpc", ch, cfg.charge_thresholds[i]),
                        &timeCharge[ch][i]);
    }
  }

  // Process all events in the specified range
  std::cout << "Analyzing " << nEntries << " events..." << std::endl;

  // Calculate progress reporting interval (report at least 10 times)
  Long64_t reportInterval = (nEntries < 10) ? 1 : nEntries / 10;

  for (Long64_t i = 0; i < nEntries; ++i) {
    Long64_t entry = startEntry + i;

    if (i % reportInterval == 0 || i == nEntries - 1) {
      std::cout << "Processing entry " << entry << " (" << i << " / " << nEntries
                << " = " << (100 * i / nEntries) << "%)" << std::endl;
    }

    inputTree->GetEntry(entry);
    event = eventIdx;

    if (!timeAxis || timeAxis->empty()) {
      std::cerr << "WARNING: empty time axis at entry " << entry << std::endl;
      continue;
    }

    // Analyze each channel
    for (int ch = 0; ch < cfg.n_channels; ++ch) {
      if (!chPed[ch] || chPed[ch]->empty()) {
        continue;
      }

      WaveformFeatures features = AnalyzeWaveform(*chPed[ch], *timeAxis, cfg, ch);

      hasSignal[ch] = features.hasSignal;
      baseline[ch] = features.baseline;
      rmsNoise[ch] = features.rmsNoise;
      noise1Point[ch] = features.noise1Point;
      ampMinBefore[ch] = features.ampMinBefore;
      ampMaxBefore[ch] = features.ampMaxBefore;
      ampMax[ch] = features.ampMax;
      charge[ch] = features.charge;
      signalOverNoise[ch] = features.signalOverNoise;
      peakTime[ch] = features.peakTime;
      riseTime[ch] = features.riseTime;
      slewRate[ch] = features.slewRate;

      for (size_t i = 0; i < nCFD && i < features.timeCFD.size(); ++i) {
        timeCFD[ch][i] = features.timeCFD[i];
        jitterCFD[ch][i] = features.jitterCFD[i];
      }
      for (size_t i = 0; i < nLE && i < features.timeLE.size(); ++i) {
        timeLE[ch][i] = features.timeLE[i];
        jitterLE[ch][i] = features.jitterLE[i];
        totLE[ch][i] = features.totLE[i];
      }
      for (size_t i = 0; i < nCharge && i < features.timeCharge.size(); ++i) {
        timeCharge[ch][i] = features.timeCharge[i];
      }

      // Save waveform plots if enabled
      if (waveformPlotsFile) {
        // Check if we should save this waveform
        bool shouldSave = !cfg.waveform_plots_only_signal || features.hasSignal;
        if (shouldSave) {
          SaveWaveformPlots(waveformPlotsFile, eventIdx, ch, *chPed[ch], *timeAxis, features, cfg);
        }
      }
    }

    // Create 2D histograms for each sensor showing signal amplitudes
    if (waveformPlotsFile) {
      // Determine which sensors are present and how many strips each has
      std::map<int, std::vector<int>> sensorStrips;  // sensor ID -> list of strip IDs
      std::map<int, std::vector<int>> sensorChannels; // sensor ID -> list of channel indices

      for (int ch = 0; ch < cfg.n_channels; ++ch) {
        int sensorID = cfg.sensor_ids[ch];
        int stripID = cfg.strip_ids[ch];
        sensorStrips[sensorID].push_back(stripID);
        sensorChannels[sensorID].push_back(ch);
      }

      // Create event directory if not exists
      char eventDirName[64];
      std::snprintf(eventDirName, sizeof(eventDirName), "event_%06d", eventIdx);
      TDirectory *eventDir = waveformPlotsFile->GetDirectory(eventDirName);
      if (!eventDir) {
        eventDir = waveformPlotsFile->mkdir(eventDirName);
      }

      // Create histogram for each sensor
      for (const auto &sensorPair : sensorStrips) {
        int sensorID = sensorPair.first;
        const std::vector<int> &strips = sensorPair.second;
        const std::vector<int> &channels = sensorChannels[sensorID];

        // Find strip range
        int minStrip = *std::min_element(strips.begin(), strips.end());
        int maxStrip = *std::max_element(strips.begin(), strips.end());
        int nStrips = maxStrip - minStrip + 1;

        // Create histogram: 1 bin in X, nStrips bins in Y
        TH2F *hist = new TH2F(Form("sensor%02d_amplitude_map", sensorID),
                              Form("Event %d - Sensor %02d Amplitude Map;X;Strip;Amplitude (V)",
                                   eventIdx, sensorID),
                              1, 0, 1,  // X axis: single bin
                              nStrips, minStrip, maxStrip + 1);  // Y axis: one bin per strip

        // Fill histogram with amplitude values
        for (size_t i = 0; i < channels.size(); ++i) {
          int ch = channels[i];
          int stripID = strips[i];
          float amplitude = ampMax[ch];  // Use maximum amplitude
          hist->Fill(0.5, stripID, amplitude);  // X=0.5 (center of bin), Y=stripID
        }

        // Save to sensor directory within event
        TDirectory *sensorDir = eventDir->GetDirectory(Form("sensor%02d", sensorID));
        if (!sensorDir) {
          sensorDir = eventDir->mkdir(Form("sensor%02d", sensorID));
        }
        sensorDir->cd();
        hist->Write(hist->GetName(), TObject::kOverwrite);
        delete hist;
      }

      waveformPlotsFile->cd();
    }

    outputTree->Fill();
  }

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
  inputFile->Close();

  // Close waveform plots file if it was created
  if (waveformPlotsFile) {
    waveformPlotsFile->cd();
    waveformPlotsFile->Close();
    delete waveformPlotsFile;
    std::string waveformPlotsFullPath = BuildOutputPath(cfg.output_dir, "waveform_plots", cfg.waveform_plots_dir + ".root");
    std::cout << "Waveform plots output saved to " << waveformPlotsFullPath << std::endl;
  }

  std::string outputFullPath = BuildOutputPath(cfg.output_dir, "root", cfg.output_root);
  std::cout << "Analysis complete. Output written to " << outputFullPath << std::endl;
  return true;
}

void PrintUsage(const char *prog) {
  std::cout << "Analyze waveforms: Extract timing and amplitude features from ROOT file\n"
            << "Usage: " << prog << " [options]\n"
            << "Options:\n"
            << "  --config PATH          Load analysis settings from JSON file\n"
            << "  --input FILE           Override input ROOT file\n"
            << "  --output FILE          Override output ROOT file\n"
            << "  --event-range START:END  Process only events in range [START, END)\n"
            << "  --waveform-plots       Enable waveform plots output (saves detailed waveform plots)\n"
            << "  --waveform-plots-file NAME  Set waveform plots output ROOT file name (default: waveform_plots.root)\n"
            << "  --waveform-plots-all   Save all waveforms (default: only with signal)\n"
            << "  -h, --help             Show this help message\n";
}

} // namespace

int main(int argc, char **argv) {
  AnalysisConfig cfg;

  // Try to load default config
  std::string defaultPath = "pipeline_config.json";
  std::string err;
  if (LoadAnalysisConfigFromJson(defaultPath, cfg, &err)) {
    std::cout << "Loaded configuration from " << defaultPath << std::endl;
  }

  // Event range parameters
  Long64_t eventStart = -1;
  Long64_t eventEnd = -1;

  // Parse command line
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      return 0;
    } else if (arg == "--config") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --config requires a value" << std::endl;
        return 1;
      }
      if (!LoadAnalysisConfigFromJson(argv[++i], cfg, &err)) {
        std::cerr << "ERROR: " << err << std::endl;
        return 1;
      }
      std::cout << "Loaded configuration from " << argv[i] << std::endl;
    } else if (arg == "--input") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --input requires a value" << std::endl;
        return 1;
      }
      cfg.input_root = argv[++i];
    } else if (arg == "--output") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --output requires a value" << std::endl;
        return 1;
      }
      cfg.output_root = argv[++i];
    } else if (arg == "--event-range") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --event-range requires a value (START:END)" << std::endl;
        return 1;
      }
      std::string range = argv[++i];
      size_t colonPos = range.find(':');
      if (colonPos == std::string::npos) {
        std::cerr << "ERROR: --event-range format must be START:END" << std::endl;
        return 1;
      }
      try {
        eventStart = std::stoll(range.substr(0, colonPos));
        eventEnd = std::stoll(range.substr(colonPos + 1));
        std::cout << "Event range: [" << eventStart << ", " << eventEnd << ")" << std::endl;
      } catch (...) {
        std::cerr << "ERROR: invalid event range format" << std::endl;
        return 1;
      }
    } else if (arg == "--waveform-plots") {
      cfg.waveform_plots_enabled = true;
      std::cout << "Waveform plots output enabled" << std::endl;
    } else if (arg == "--waveform-plots-file") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --waveform-plots-file requires a value" << std::endl;
        return 1;
      }
      cfg.waveform_plots_dir = argv[++i];
      // Remove .root extension if provided
      if (cfg.waveform_plots_dir.size() > 5 &&
          cfg.waveform_plots_dir.substr(cfg.waveform_plots_dir.size() - 5) == ".root") {
        cfg.waveform_plots_dir = cfg.waveform_plots_dir.substr(0, cfg.waveform_plots_dir.size() - 5);
      }
    } else if (arg == "--waveform-plots-all") {
      cfg.waveform_plots_only_signal = false;
      std::cout << "Will save all waveforms (not just signals)" << std::endl;
    } else {
      std::cerr << "ERROR: unknown option " << arg << std::endl;
      PrintUsage(argv[0]);
      return 1;
    }
  }

  try {
    if (!RunAnalysis(cfg, eventStart, eventEnd)) {
      return 2;
    }
  } catch (const std::exception &ex) {
    std::cerr << "Unhandled exception: " << ex.what() << std::endl;
    return 3;
  }

  return 0;
}
