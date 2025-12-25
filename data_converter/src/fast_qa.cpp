#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include "config/analysis_config.h"
#include "utils/filesystem_utils.h"

using namespace std;

namespace {

std::string to6digits(int n) {
  std::ostringstream oss;
  oss << std::setw(6) << std::setfill('0') << n;
  return oss.str();
}

bool EnsureParentDirectory(const std::string &path) {
  size_t lastSlash = path.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = path.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      return false;
    }
  }
  return true;
}

bool RunFastQA(const AnalysisConfig &cfg) {
  gStyle->SetOptStat(0);  // Disable statistics box

  string outname_base = cfg.output_dir() + '/';
  outname_base += to6digits(cfg.runnumber()) + '/';
  outname_base += cfg.daq_name() + "/output/";

  // Build input path: output_dir/root/analysis_root
  std::string inputPath = BuildOutputPath(outname_base, "root", cfg.output_root());

  // Open analyzed ROOT file
  TFile *inputFile = TFile::Open(inputPath.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open analyzed ROOT file " << inputPath << std::endl;
    return false;
  }
  std::cout << "Reading analyzed file: " << inputPath << std::endl;

  TTree *tree = dynamic_cast<TTree *>(inputFile->Get(cfg.output_tree().c_str()));
  if (!tree) {
    std::cerr << "ERROR: cannot find tree " << cfg.output_tree() << std::endl;
    inputFile->Close();
    return false;
  }

  // Set up branch addresses
  int event = 0;
  int nChannels = 0;
  std::vector<int> *sensorID = nullptr;
  std::vector<int> *sensorCol = nullptr;
  std::vector<int> *sensorRow = nullptr;
  std::vector<bool> *isHorizontal = nullptr;
  std::vector<float> *ampMax = nullptr;
  std::vector<float> *baseline = nullptr;

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("nChannels", &nChannels);
  tree->SetBranchAddress("sensorID", &sensorID);
  tree->SetBranchAddress("sensorCol", &sensorCol);
  tree->SetBranchAddress("sensorRow", &sensorRow);
  tree->SetBranchAddress("isHorizontal", &isHorizontal);
  tree->SetBranchAddress("ampMax", &ampMax);
  tree->SetBranchAddress("baseline", &baseline);

  // Create output file
  std::string qualityCheckFileName = BuildOutputPath(outname_base, "quality_check",
                                                     "quality_check.root");
  if (!EnsureParentDirectory(qualityCheckFileName)) {
    std::cerr << "ERROR: failed to create quality_check output directory for "
              << qualityCheckFileName << std::endl;
    inputFile->Close();
    return false;
  }

  TFile *outputFile = TFile::Open(qualityCheckFileName.c_str(), "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "ERROR: cannot create quality_check output file "
              << qualityCheckFileName << std::endl;
    inputFile->Close();
    return false;
  }
  std::cout << "Creating quality check file: " << qualityCheckFileName << std::endl;

  // Create histograms for accumulation across all events
  std::map<int, TH1F*> ampMaxHists;
  std::map<int, TH1F*> baselineHists;

  int nCh = cfg.n_channels();
  for (int ch = 0; ch < nCh; ++ch) {
    ampMaxHists[ch] = new TH1F(Form("ampMax_ch%02d", ch),
                                Form("Channel %d Amplitude;Amplitude (ADC);Events", ch),
                                500, 0, 5000);  // ADC range: 0-5000
    baselineHists[ch] = new TH1F(Form("baseline_ch%02d", ch),
                                  Form("Channel %d Baseline;Baseline (ADC);Events", ch),
                                  200, 3400, 3600);
  }

  // Process all events
  Long64_t nEntries = tree->GetEntries();
  std::cout << "Processing " << nEntries << " events for quality check..." << std::endl;

  Long64_t reportInterval = (nEntries < 10) ? 1 : nEntries / 10;

  // Create event directories in output file
  TDirectory *eventsDir = outputFile->mkdir("events");
  if (!eventsDir) {
    std::cerr << "ERROR: failed to create events directory" << std::endl;
    outputFile->Close();
    inputFile->Close();
    return false;
  }

  for (Long64_t i = 0; i < nEntries; ++i) {
    if (i % reportInterval == 0 || i == nEntries - 1) {
      std::cout << "Processing entry " << i << " / " << nEntries
                << " (" << (100 * i / nEntries) << "%)" << std::endl;
    }

    tree->GetEntry(i);

    // Debug: print first event's ampMax values
    if (i == 0) {
      std::cout << "DEBUG: First event ampMax values:" << std::endl;
      if (ampMax && ampMax->size() > 0) {
        for (int ch = 0; ch < std::min(nCh, static_cast<int>(ampMax->size())); ++ch) {
          std::cout << "  ch" << ch << ": " << ampMax->at(ch) << " ADC" << std::endl;
        }
      } else {
        std::cout << "  WARNING: ampMax is empty or null!" << std::endl;
      }
    }

    // Fill histograms
    for (int ch = 0; ch < nCh; ++ch) {
      if (ch < static_cast<int>(ampMax->size())) {
        ampMaxHists[ch]->Fill(ampMax->at(ch));
      }
      if (ch < static_cast<int>(baseline->size())) {
        baselineHists[ch]->Fill(baseline->at(ch));
      }
    }

    // Create amplitude maps per sensor for this event
    std::map<int, std::vector<int>> sensorChannels;
    for (int ch = 0; ch < nCh; ++ch) {
      if (ch < static_cast<int>(sensorID->size())) {
        int sid = sensorID->at(ch);
        sensorChannels[sid].push_back(ch);
      }
    }

    // Determine max sensor ID for canvas sizing
    int maxSensorID = 0;
    for (const auto &sensorPair : sensorChannels) {
      if (sensorPair.first > maxSensorID) {
        maxSensorID = sensorPair.first;
      }
    }

    int numPads = std::max(4, maxSensorID + 1);

    // Create canvas with all sensors
    char canvasName[64];
    std::snprintf(canvasName, sizeof(canvasName), "event_%06d_quality_check", event);
    TCanvas *canvas = new TCanvas(canvasName,
                                  Form("Event %d - All Sensors Quality Check", event),
                                  600 * numPads, 800);
    canvas->Divide(numPads, 1);

    // Store histograms to keep them alive until canvas is written
    std::vector<TH2F*> eventHistograms;

    // Create amplitude map for each sensor
    for (const auto &sensorPair : sensorChannels) {
      int sid = sensorPair.first;
      const std::vector<int> &channels = sensorPair.second;

      if (channels.empty()) continue;

      bool isHoriz = isHorizontal->at(channels[0]);

      TH2F *hist = nullptr;
      if (!isHoriz) {
        // Vertical: X=Column (sensor_row), Y=Strip (sensor_col)
        // Bin edges at -0.5, 0.5, 1.5 so that integer values 0, 1 are at bin centers
        hist = new TH2F(Form("sensor%02d_amplitude_map", sid),
                        Form("Event %d - Sensor %02d Amplitude Map;Column;Strip;Amplitude (ADC)",
                             event, sid),
                        2, -0.5, 1.5,  // X axis: column (bin centers at 0, 1)
                        5, -0.5, 4.5);  // Y axis: strip (bin centers at 0, 1, 2, 3, 4)
      } else {
        // Horizontal: X=Strip (sensor_col), Y=Column (sensor_row)
        hist = new TH2F(Form("sensor%02d_amplitude_map", sid),
                        Form("Event %d - Sensor %02d Amplitude Map;Strip;Column;Amplitude (ADC)",
                             event, sid),
                        5, -0.5, 4.5,  // X axis: strip (bin centers at 0, 1, 2, 3, 4)
                        2, -0.5, 1.5);  // Y axis: column (bin centers at 0, 1)
      }

      // Fill histogram with amplitude values
      for (int ch : channels) {
        if (ch >= static_cast<int>(sensorCol->size()) ||
            ch >= static_cast<int>(sensorRow->size()) ||
            ch >= static_cast<int>(ampMax->size())) {
          continue;
        }

        int strip = sensorCol->at(ch);
        int col = sensorRow->at(ch);
        float amplitude = ampMax->at(ch);

        // Debug first event
        if (i == 0 && ch < 4) {
          std::cout << "DEBUG: Event " << event << " ch" << ch
                    << " -> sensor=" << sid << " col=" << col << " strip=" << strip
                    << " amp=" << amplitude << " ADC" << std::endl;
        }

        if (!isHoriz) {
          hist->Fill(col, strip, amplitude);
        } else {
          hist->Fill(strip, col, amplitude);
        }
      }

      // Debug: check histogram content
      if (i == 0) {
        std::cout << "DEBUG: Sensor " << sid << " histogram entries=" << hist->GetEntries()
                  << " max=" << hist->GetMaximum() << std::endl;
      }

      // Draw in corresponding pad
      canvas->cd(sid + 1);
      hist->SetStats(0);

      // Set fixed Z-axis range for consistent comparison across events
      hist->SetMinimum(0);
      hist->SetMaximum(300);  // Fixed maximum for all events

      hist->Draw("COLZ TEXT");

      // Store histogram to keep it alive
      eventHistograms.push_back(hist);
    }

    // Save canvas to events directory
    eventsDir->cd();
    canvas->Write(canvasName, TObject::kOverwrite);
    delete canvas;

    // Now safe to delete histograms
    for (auto *h : eventHistograms) {
      delete h;
    }
  }

  // Return to main directory
  outputFile->cd();

  std::cout << "Creating summary histograms..." << std::endl;

  // Create overlay canvas for amplitude
  TCanvas *cAmpMax = new TCanvas("ampMax_all_channels",
                                 "All Channels Amplitude Distribution",
                                 1200, 800);
  cAmpMax->SetGrid();
  cAmpMax->SetLogy();  // Log scale for better visibility

  TLegend *legAmp = new TLegend(0.75, 0.4, 0.95, 0.9);
  legAmp->SetTextSize(0.025);

  // Define color palette
  const int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2,
                        kViolet, kSpring-5, kAzure+7, kPink-3, kTeal-5, kYellow+2,
                        kRed+2, kBlue+2, kGreen-2, kMagenta+2};

  double maxY = 0;
  for (int ch = 0; ch < nCh; ++ch) {
    double thisMax = ampMaxHists[ch]->GetMaximum();
    if (thisMax > maxY) maxY = thisMax;
  }

  for (int ch = 0; ch < nCh; ++ch) {
    ampMaxHists[ch]->SetLineColor(colors[ch % 16]);
    ampMaxHists[ch]->SetLineWidth(2);
    ampMaxHists[ch]->GetYaxis()->SetRangeUser(0.5, maxY * 1.5);
    ampMaxHists[ch]->GetXaxis()->SetTitle("Amplitude (ADC)");
    ampMaxHists[ch]->GetYaxis()->SetTitle("Events");
    ampMaxHists[ch]->Draw(ch == 0 ? "HIST" : "HIST SAME");
    legAmp->AddEntry(ampMaxHists[ch],
                     Form("%s_ch%02d", cfg.daq_name().c_str(), ch), "l");
  }
  legAmp->Draw();
  cAmpMax->Write();

  // Create overlay canvas for baseline
  TCanvas *cBaseline = new TCanvas("baseline_all_channels",
                                   "All Channels Baseline Distribution",
                                   1200, 800);
  cBaseline->SetGrid();

  TLegend *legBase = new TLegend(0.15, 0.15, 0.45, 0.55);
  legBase->SetTextSize(0.022);
  legBase->SetNColumns(2);

  maxY = 0;
  for (int ch = 0; ch < nCh; ++ch) {
    double thisMax = baselineHists[ch]->GetMaximum();
    if (thisMax > maxY) maxY = thisMax;
  }

  for (int ch = 0; ch < nCh; ++ch) {
    baselineHists[ch]->SetLineColor(colors[ch % 16]);
    baselineHists[ch]->SetLineWidth(2);
    baselineHists[ch]->GetYaxis()->SetRangeUser(0, maxY * 1.2);
    baselineHists[ch]->GetXaxis()->SetTitle("Baseline (ADC)");
    baselineHists[ch]->GetYaxis()->SetTitle("Events");
    baselineHists[ch]->Draw(ch == 0 ? "HIST" : "HIST SAME");

    double mean = baselineHists[ch]->GetMean();
    double sigma = baselineHists[ch]->GetStdDev();
    legBase->AddEntry(baselineHists[ch],
                      Form("%s_ch%02d: %.1f#pm%.1f",
                           cfg.daq_name().c_str(), ch, mean, sigma), "l");
  }
  legBase->Draw();
  cBaseline->Write();

  // Write individual histograms
  std::cout << "Writing individual histograms..." << std::endl;
  for (int ch = 0; ch < nCh; ++ch) {
    ampMaxHists[ch]->Write();
    baselineHists[ch]->Write();
  }

  // Clean up
  for (auto &pair : ampMaxHists) delete pair.second;
  for (auto &pair : baselineHists) delete pair.second;
  delete cAmpMax;
  delete cBaseline;

  outputFile->Close();
  inputFile->Close();

  std::cout << "Quality check complete. Output written to " << qualityCheckFileName << std::endl;
  return true;
}

void PrintUsage(const char *prog) {
  std::cout << "Fast QA: Generate quality check plots from analyzed waveforms\n"
            << "Usage: " << prog << " [options]\n"
            << "Options:\n"
            << "  --config PATH    Load settings from JSON file\n"
            << "  -h, --help       Show this help message\n";
}

} // namespace

int main(int argc, char **argv) {
  AnalysisConfig cfg;

  // Try to load default config
  std::string defaultPath = "converter_config.json";
  std::string err;
  if (LoadAnalysisConfigFromJson(defaultPath, cfg, &err)) {
    std::cout << "Loaded configuration from " << defaultPath << std::endl;
  }

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
    } else {
      std::cerr << "ERROR: unknown option " << arg << std::endl;
      PrintUsage(argv[0]);
      return 1;
    }
  }

  try {
    if (!RunFastQA(cfg)) {
      return 2;
    }
  } catch (const std::exception &ex) {
    std::cerr << "Unhandled exception: " << ex.what() << std::endl;
    return 3;
  }

  return 0;
}
