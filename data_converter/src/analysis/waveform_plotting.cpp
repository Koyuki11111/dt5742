#include "analysis/waveform_plotting.h"

#include <cstdio>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"

void SaveWaveformPlots(TFile *waveformPlotsFile, int event, int channel,
                       const std::vector<float> &amp,
                       const std::vector<float> &time,
                       const WaveformFeatures &features,
                       const AnalysisConfig &cfg) {
  if (!waveformPlotsFile || waveformPlotsFile->IsZombie()) {
    return;
  }

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
  // In the unified naming scheme, sensor_cols holds variable that acts as strip ID
  const int stripID = cfg.sensor_cols[channel];

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

  const int nAnalysisPoints = analysisEndIdx - analysisStartIdx + 1;

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

  char avgGraphName[64];
  std::snprintf(avgGraphName, sizeof(avgGraphName), "strip%02d_avg", stripID);

  TGraph *grAvg = new TGraph(nAnalysisPoints);
  grAvg->SetName(avgGraphName);
  grAvg->SetTitle(Form("Event %d, Sensor %d, Strip %d (Ch%d) - 3-Point Moving Average;Time (ns);Amplitude (V)", event, sensorID, stripID, channel));

  pointIdx = 0;
  for (int i = analysisStartIdx; i <= analysisEndIdx; ++i) {
    float avgAmp = 0.0f;
    if (i == 0 || i == analysisStartIdx) {
      if (i + 1 < nSamples) {
        avgAmp = ((amp[i] - features.baseline) + (amp[i+1] - features.baseline)) / 2.0f;
      } else {
        avgAmp = (amp[i] - features.baseline);
      }
    } else if (i == nSamples - 1 || i == analysisEndIdx) {
      avgAmp = ((amp[i-1] - features.baseline) + (amp[i] - features.baseline)) / 2.0f;
    } else {
      avgAmp = ((amp[i-1] - features.baseline) + (amp[i] - features.baseline) + (amp[i+1] - features.baseline)) / 3.0f;
    }
    grAvg->SetPoint(pointIdx++, time[i], avgAmp);
  }
  grAvg->Write(avgGraphName, TObject::kOverwrite);

  char canvasName[64];
  std::snprintf(canvasName, sizeof(canvasName), "strip%02d_analysis", stripID);

  TCanvas *canvas = new TCanvas(canvasName,
                                Form("Event %d, Sensor %d, Strip %d (Ch%d) - Analysis", event, sensorID, stripID, channel),
                                1200, 800);
  canvas->SetGrid();

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

  TLine *baselineLine = new TLine(time[analysisStartIdx], 0.0,
                                  time[analysisEndIdx], 0.0);
  baselineLine->SetLineStyle(2);
  baselineLine->SetLineColor(kRed);
  baselineLine->Draw("SAME");

  float displayPeakTime = features.peakTime;
  float displayAmp = features.ampMax * ((polarity > 0) ? 1.0f : -1.0f);
  TMarker *peakMarker = new TMarker(displayPeakTime, displayAmp, 20);
  peakMarker->SetMarkerColor(kBlack);
  peakMarker->SetMarkerSize(1.2);
  peakMarker->Draw("SAME");

  TLegend *legend = new TLegend(0.65, 0.70, 0.90, 0.90);
  legend->AddEntry(grAnalysis, "Waveform", "l");
  legend->AddEntry(baselineLine, "Baseline", "l");

  for (size_t i = 0; i < features.timeCFD.size(); ++i) {
    float tCFD = features.timeCFD[i];
    if (tCFD <= 0.0f) continue;

    float frac = (i < cfg.cfd_thresholds.size()) ? cfg.cfd_thresholds[i] : 0;
    TLine *cfdLine = new TLine(tCFD, 0.0, tCFD, displayAmp);
    cfdLine->SetLineColor(kGreen + static_cast<int>(i % 5));
    cfdLine->SetLineStyle(7);
    cfdLine->Draw("SAME");
    legend->AddEntry(cfdLine, Form("CFD %.0f%%", frac), "l");
  }

  for (size_t i = 0; i < features.timeLE.size(); ++i) {
    float tLE = features.timeLE[i];
    if (tLE <= 0.0f) continue;

    float threshold = (i < cfg.le_thresholds.size()) ? cfg.le_thresholds[i] : 0;
    TLine *leLine = new TLine(tLE, 0.0, tLE, displayAmp);
    leLine->SetLineColor(kMagenta + static_cast<int>(i % 3));
    leLine->SetLineStyle(9);
    leLine->Draw("SAME");
    legend->AddEntry(leLine, Form("LE %.1f mV", threshold), "l");
  }

  for (size_t i = 0; i < features.timeCharge.size(); ++i) {
    float tCharge = features.timeCharge[i];
    if (tCharge <= 0.0f) continue;

    float threshold = (i < cfg.charge_thresholds.size()) ? cfg.charge_thresholds[i] : 0;
    TMarker *chargeMarker = new TMarker(tCharge, displayAmp * 0.8f, 21);
    chargeMarker->SetMarkerColor(kOrange + static_cast<int>(i % 3));
    chargeMarker->SetMarkerSize(1.0);
    chargeMarker->Draw("SAME");
    legend->AddEntry(chargeMarker, Form("Charge %.0f%%", threshold), "p");
  }

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

  delete grRaw;
  delete grAvg;
  delete grAnalysis;
  delete baselineLine;
  delete legend;
  delete infoText;
  delete canvas;

  waveformPlotsFile->cd();
}
