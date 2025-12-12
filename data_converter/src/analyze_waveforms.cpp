#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "utils/file_io.h"
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

#include "config/analysis_config.h"
#include "analysis/waveform_math.h"
#include "analysis/waveform_plotting.h"

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

enum class NsamplesPolicy { kStrict, kPad };

NsamplesPolicy ResolveNsamplesPolicy(const std::string &policyText) {
  std::string lowered = policyText;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lowered == "pad") {
    return NsamplesPolicy::kPad;
  }
  if (lowered != "strict") {
    std::cerr << "WARNING: unknown nsamples_policy '" << policyText
              << "', defaulting to 'strict'" << std::endl;
  }
  return NsamplesPolicy::kStrict;
}

// Helper to check if sensor should be displayed horizontally
bool IsSensorHorizontal(int sensorID, const AnalysisConfig& cfg) {
  // Find unique sensor IDs in current config and map to local index
  std::set<int> uniqueSensors(cfg.sensor_ids.begin(), cfg.sensor_ids.end());
  std::vector<int> sortedSensors(uniqueSensors.begin(), uniqueSensors.end());
  std::sort(sortedSensors.begin(), sortedSensors.end());
  
  // Find index of this sensorID in the sorted unique list
  auto it = std::find(sortedSensors.begin(), sortedSensors.end(), sensorID);
  if (it == sortedSensors.end()) {
    return false;  // Sensor not found
  }

  int localIndex = std::distance(sortedSensors.begin(), it);
  if (localIndex < 0 || localIndex >= static_cast<int>(cfg.sensor_orientations.size())) {
    return false;  // Out of bounds, default to vertical
  }

  return cfg.sensor_orientations[localIndex] == "horizontal";
}

bool RunAnalysis(const AnalysisConfig &cfg, Long64_t eventStart = -1, Long64_t eventEnd = -1) {
  // Maximum file size for waveform plots: 4 GB
  const Long64_t MAX_PLOTS_FILE_SIZE = 4LL * 1024 * 1024 * 1024;  // 4 GB in bytes
  
  // Create waveform plots ROOT file if enabled
  TFile *waveformPlotsFile = nullptr;
  int waveformPlotsFileCounter = 0;

  string outname_base = cfg.output_dir()+'/';
  outname_base += to6digits(cfg.runnumber())+'/';
  outname_base += cfg.daq_name()+"/output/";
  
  auto openWaveformPlotsFile = [&](TFile *&outFile, int fileNum) {
    if (!cfg.waveform_plots_enabled) {
      return;
    }

    std::string baseFileName = cfg.waveform_plots_dir;
    std::string waveformPlotsFileName;
    
    if (fileNum == 0) {
      waveformPlotsFileName = BuildOutputPath(outname_base, "waveform_plots",
                                               baseFileName + ".root");
    } else {
      char suffix[32];
      std::snprintf(suffix, sizeof(suffix), "_%03d.root", fileNum);
      waveformPlotsFileName = BuildOutputPath(outname_base, "waveform_plots",
                                               baseFileName + suffix);
    }

    if (!EnsureParentDirectory(waveformPlotsFileName)) {
      std::cerr << "WARNING: Failed to create waveform plots output directory for "
                << waveformPlotsFileName << std::endl;
      return;
    }
    outFile = TFile::Open(waveformPlotsFileName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
      std::cerr << "WARNING: Failed to create waveform plots output file "
                << waveformPlotsFileName << std::endl;
      std::cerr << "         Continuing without waveform plots output..." << std::endl;
      outFile = nullptr;
      return;
    }
    std::cout << "Waveform plots output enabled. Saving to: " << waveformPlotsFileName << std::endl;
    if (cfg.waveform_plots_only_signal && fileNum == 0) {
      std::cout << "  Only saving waveforms with detected signals (SNR > "
                << cfg.snr_threshold << ")" << std::endl;
    }
  };
  
  auto checkAndRotateWaveformPlotsFile = [&](TFile *&outFile) {
    if (!outFile || outFile->IsZombie()) {
      return;
    }

    Long64_t currentSize = outFile->GetSize();
    if (currentSize >= MAX_PLOTS_FILE_SIZE) {
      std::cout << "Waveform plots file size reached " << (currentSize / (1024.0 * 1024.0 * 1024.0))
                << " GB. Rotating to new file..." << std::endl;

      // Close current file
      std::string currentFileName = outFile->GetName();
      outFile->Close();
      delete outFile;
      outFile = nullptr;

      std::cout << "Saved waveform plots to: " << currentFileName << std::endl;

      // Open new file with incremented counter
      waveformPlotsFileCounter++;
      openWaveformPlotsFile(outFile, waveformPlotsFileCounter);
    }
  };

  openWaveformPlotsFile(waveformPlotsFile, waveformPlotsFileCounter);

  // Create quality check ROOT file
  TFile *qualityCheckFile = nullptr;
  if (cfg.waveform_plots_enabled) {
    // Use same naming scheme as waveform_plots_dir for quality check
    // If waveform_plots_dir is "waveform_plots", use "quality_check"
    // If waveform_plots_dir is "waveform_plots_chunk_0", use "quality_check_chunk_0"
    std::string qualityCheckBaseName = cfg.waveform_plots_dir;

    // Replace "waveform_plots" with "quality_check" in the base name
    size_t pos = qualityCheckBaseName.find("waveform_plots");
    if (pos != std::string::npos) {
      qualityCheckBaseName.replace(pos, std::string("waveform_plots").length(), "quality_check");
    } else {
      // Fallback: just use "quality_check" prefix
      qualityCheckBaseName = "quality_check_" + qualityCheckBaseName;
    }
    
    std::string qualityCheckFileName = BuildOutputPath(outname_base, "quality_check",
                                                       qualityCheckBaseName + ".root");
    if (!EnsureParentDirectory(qualityCheckFileName)) {
      std::cerr << "WARNING: Failed to create quality_check output directory for "
                << qualityCheckFileName << std::endl;
    } else {
      qualityCheckFile = TFile::Open(qualityCheckFileName.c_str(), "RECREATE");
      if (!qualityCheckFile || qualityCheckFile->IsZombie()) {
        std::cerr << "WARNING: Failed to create quality_check output file "
                  << qualityCheckFileName << std::endl;
        std::cerr << "         Continuing without quality_check output..." << std::endl;
        qualityCheckFile = nullptr;
      } else {
        std::cout << "Quality check output enabled. Saving to: " << qualityCheckFileName << std::endl;
      }
    }
  }

  // Build input path: output_dir/root/input_root
  std::string inputPath = BuildOutputPath(outname_base, "root", cfg.input_root());

  // Open input ROOT file
  TFile *inputFile = TFile::Open(inputPath.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open input ROOT file " << inputPath << std::endl;
    return false;
  }
  std::cout << "Reading input file: " << inputPath << std::endl;

  TTree *inputTree = dynamic_cast<TTree *>(inputFile->Get(cfg.input_tree().c_str()));
  if (!inputTree) {
    std::cerr << "ERROR: cannot find tree " << cfg.input_tree() << std::endl;
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

  const NsamplesPolicy policy = ResolveNsamplesPolicy(cfg.common.nsamples_policy);

  // Set up input branches
  int eventIdx = 0;
  int nChannels = 0;
  int nsamples = 0;
  std::vector<float> *timeAxis = nullptr;
  std::vector<std::vector<float> *> chPed(cfg.n_channels(), nullptr);
  std::vector<int> *nsamplesPerChannel = nullptr;

  inputTree->SetBranchAddress("event", &eventIdx);
  inputTree->SetBranchAddress("n_channels", &nChannels);
  inputTree->SetBranchAddress("nsamples", &nsamples);
  inputTree->SetBranchAddress("time_ns", &timeAxis);
  if (inputTree->GetBranch("nsamples_per_channel")) {
    inputTree->SetBranchAddress("nsamples_per_channel", &nsamplesPerChannel);
  }

  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
    char bname[32];
    std::snprintf(bname, sizeof(bname), "ch%02d_ped", ch);
    if (inputTree->GetBranch(bname)) {
      inputTree->SetBranchAddress(bname, &chPed[ch]);
    }
  }

  // Build output path: output_dir/root/output_root
  std::string outputPath = BuildOutputPath(outname_base, "root", cfg.output_root());
  
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

  TTree *outputTree = new TTree(cfg.output_tree().c_str(), "Analyzed waveform features");

  // Create output branches - per channel vectors
  int event = 0;
  std::vector<int> sensorID(cfg.n_channels());
  std::vector<int> sensorCol(cfg.n_channels());
  std::vector<int> sensorRow(cfg.n_channels());
  std::vector<int> stripID(cfg.n_channels());
  std::vector<bool> isHorizontal(cfg.n_channels());
  std::vector<bool> hasSignal(cfg.n_channels());
  std::vector<float> baseline(cfg.n_channels());
  std::vector<float> rmsNoise(cfg.n_channels());
  std::vector<float> noise1Point(cfg.n_channels());
  std::vector<float> ampMinBefore(cfg.n_channels());
  std::vector<float> ampMaxBefore(cfg.n_channels());
  std::vector<float> ampMax(cfg.n_channels());
  std::vector<float> charge(cfg.n_channels());
  std::vector<float> signalOverNoise(cfg.n_channels());
  std::vector<float> peakTime(cfg.n_channels());
  std::vector<float> riseTime(cfg.n_channels());
  std::vector<float> slewRate(cfg.n_channels());
  std::vector<float> jitterRMS(cfg.n_channels());
  
  // Initialize sensor and strip IDs from config
  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
    sensorID[ch]  = cfg.sensor_ids[ch];
    sensorCol[ch] = cfg.sensor_cols[ch];
    sensorRow[ch] = cfg.sensor_rows[ch];
    stripID[ch]   = cfg.strip_ids[ch];
    // Check sensor orientation
    isHorizontal[ch] = IsSensorHorizontal(sensorID[ch], cfg);
  }
  
  auto defineScalarBranches = [&]() {
    outputTree->Branch("nChannels", &nChannels);
    outputTree->Branch("event", &event);
    outputTree->Branch("sensorID", &sensorID);
    outputTree->Branch("sensorCol", &sensorCol);
    outputTree->Branch("sensorRow", &sensorRow);
    outputTree->Branch("stripID", &stripID);
    outputTree->Branch("isHorizontal", &isHorizontal);
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
    outputTree->Branch("jitterRMS", &jitterRMS);
  };

  // Multi-threshold timing branches (per channel, per threshold)
  const size_t nCFD = cfg.cfd_thresholds.size();
  const size_t nLE = cfg.le_thresholds.size();
  const size_t nCharge = cfg.charge_thresholds.size();

  std::vector<std::vector<float>> timeCFD(cfg.n_channels(), std::vector<float>(nCFD));
  std::vector<std::vector<float>> jitterCFD(cfg.n_channels(), std::vector<float>(nCFD));
  std::vector<std::vector<float>> timeLE(cfg.n_channels(), std::vector<float>(nLE));
  std::vector<std::vector<float>> jitterLE(cfg.n_channels(), std::vector<float>(nLE));
  std::vector<std::vector<float>> totLE(cfg.n_channels(), std::vector<float>(nLE));
  std::vector<std::vector<float>> timeCharge(cfg.n_channels(), std::vector<float>(nCharge));

  auto defineTimingBranches = [&]() {
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
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
  };

  defineScalarBranches();
  defineTimingBranches();

  // Process all events in the specified range
  std::cout << "Analyzing " << nEntries << " events..." << std::endl;

  // Calculate progress reporting interval (report at least 10 times)
  Long64_t reportInterval = (nEntries < 10) ? 1 : nEntries / 10;
  bool nsamplesError = false;
  bool loggedNsamplesTrim = false;
  std::vector<float> trimmedAmpBuf;
  std::vector<float> trimmedTimeBuf;

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

    std::vector<int> effectiveSamples(cfg.n_channels(), 0);
    int minSamples = std::numeric_limits<int>::max();
    int maxSamples = 0;
    bool haveSamples = false;

    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (!chPed[ch] || chPed[ch]->empty()) {
        continue;
      }

      int chSamples = nsamples;
      if (nsamplesPerChannel &&
          ch < static_cast<int>(nsamplesPerChannel->size())) {
        chSamples = nsamplesPerChannel->at(ch);
      }

      chSamples =
          std::min(chSamples, static_cast<int>(chPed[ch]->size()));
      chSamples =
          std::min(chSamples, static_cast<int>(timeAxis->size()));

      if (chSamples <= 0) {
        continue;
      }

      haveSamples = true;
      effectiveSamples[ch] = chSamples;
      minSamples = std::min(minSamples, chSamples);
      maxSamples = std::max(maxSamples, chSamples);
    }

    if (!haveSamples) {
      continue;
    }

    const bool hasMismatch = maxSamples != minSamples;
    if (hasMismatch && policy == NsamplesPolicy::kStrict) {
      std::cerr << "ERROR: nsamples mismatch at entry " << entry
                << " (min " << minSamples << ", max " << maxSamples << ")"
                << std::endl;
      nsamplesError = true;
      break;
    }

    if (hasMismatch && !loggedNsamplesTrim &&
        policy == NsamplesPolicy::kPad) {
      std::cout << "INFO: nsamples mismatch detected at entry " << entry
                << ", trimming analysis to per-channel sample counts" << std::endl;
      loggedNsamplesTrim = true;
    }

    // Analyze each channel
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (!chPed[ch] || chPed[ch]->empty()) {
        continue;
      }

      const int samplesToUse = effectiveSamples[ch];
      if (samplesToUse <= 0) {
        continue;
      }

      const bool needsTrim =
          samplesToUse != static_cast<int>(chPed[ch]->size()) ||
          static_cast<size_t>(samplesToUse) != timeAxis->size();

      const std::vector<float> *ampPtr = chPed[ch];
      const std::vector<float> *timePtr = timeAxis;

      if (needsTrim) {
        trimmedAmpBuf.assign(chPed[ch]->begin(),
                             chPed[ch]->begin() + samplesToUse);
        trimmedTimeBuf.assign(timeAxis->begin(),
                              timeAxis->begin() + samplesToUse);
        ampPtr = &trimmedAmpBuf;
        timePtr = &trimmedTimeBuf;
      }

      WaveformFeatures features = AnalyzeWaveform(*ampPtr, *timePtr, cfg, ch);

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
      jitterRMS[ch] = features.jitterRMS;
      
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
    if (waveformPlotsFile || qualityCheckFile) {
      // Determine which sensors are present and how many strips each has
      std::map<int, std::vector<int>> sensorStrips;   // sensor ID -> list of strip IDs
      std::map<int, std::vector<int>> sensorChannels; // sensor ID -> list of channel indices
      std::map<int, std::vector<int>> sensorCols;     // sensor ID -> list of colum indices
      std::map<int, std::vector<int>> sensorRows;     // sensor ID -> list of row indices

      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        int sensorID = cfg.sensor_ids[ch];
	int sensorCol= cfg.sensor_cols[ch];
	int sensorRow= cfg.sensor_rows[ch];
	int stripID = cfg.strip_ids[ch];
        sensorStrips[sensorID].push_back(stripID);
        sensorChannels[sensorID].push_back(ch);
	sensorCols[sensorID].push_back(sensorCol);
	sensorRows[sensorID].push_back(sensorRow);
      }

      // Create event directory if not exists (for waveformPlotsFile)
      char eventDirName[64];
      std::snprintf(eventDirName, sizeof(eventDirName), "event_%06d", eventIdx);
      TDirectory *eventDir = nullptr;
      if (waveformPlotsFile) {
        eventDir = waveformPlotsFile->GetDirectory(eventDirName);
        if (!eventDir) {
          eventDir = waveformPlotsFile->mkdir(eventDirName);
        }
      }

      // Store histograms for quality check canvas
      std::map<int, TH2F*> sensorHistograms;

      // Create histogram for each sensor
      for (const auto &sensorPair : sensorChannels) {
        int sensorID = sensorPair.first;
        const std::vector<int> &strips = sensorPair.second;
        const std::vector<int> &channels = sensorChannels[sensorID];
	
        // Find strip range
        int minStrip = *std::min_element(strips.begin(), strips.end());
        int maxStrip = *std::max_element(strips.begin(), strips.end());
        //int nStrips = maxStrip - minStrip + 1;

	// Check sensor orientation
        bool isHorizontal = IsSensorHorizontal(sensorID, cfg);
	
        TH2F *hist = nullptr;

        if (isHorizontal) {	  
          // Vertical: X=1 bin, Y=strips (current behavior)
          hist = new TH2F(Form("sensor%02d_amplitude_map", sensorID),
                          Form("Event %d - Sensor %02d Amplitude Map;X;Strip;Amplitude (V)",
                               eventIdx, sensorID),
                          2, 0, 2,  // X axis: single bin
                          64, 0, 64);  // Y axis: strips
	} else {
	  // Horizontal: X=strips, Y=1 bin (rotated 90 degrees)
          hist = new TH2F(Form("sensor%02d_amplitude_map", sensorID),
                          Form("Event %d - Sensor %02d Amplitude Map;Strip;Y;Amplitude (V)",
                               eventIdx, sensorID),
                          64, 0, 64,  // X axis: strips
                          2, 0, 2);  // Y axis: single bin
        }
	
	std::vector<float> weight_col_pos;
	
        // Fill histogram with amplitude values
        for (size_t i = 0; i < channels.size(); ++i) {
          int ch = channels[i];
          //int stripID = strips[i];	  
	  float amplitude = ampMax[ch];
	  if (isHorizontal)
	    hist->Fill(sensorRow[i], sensorCol[i], amplitude);
	  else
	    hist->Fill(sensorCol[i], sensorRow[i], amplitude);
        }
	
        // Save to sensor directory within event (for waveformPlotsFile)
        if (waveformPlotsFile && eventDir) {
          TDirectory *sensorDir = eventDir->GetDirectory(Form("sensor%02d", sensorID));
          if (!sensorDir) {
            sensorDir = eventDir->mkdir(Form("sensor%02d", sensorID));
          }
          sensorDir->cd();	  
	  hist->Write(hist->GetName(), TObject::kOverwrite);
	}
        // Store histogram for quality check canvas	
	sensorHistograms[sensorID] = hist;
      }
      // Create quality check canvas with all 4 sensors
      if (qualityCheckFile) {
        qualityCheckFile->cd();
        // Create canvas with 4 columns (1x4 layout)
        char canvasName[64];
        std::snprintf(canvasName, sizeof(canvasName), "event_%06d_quality_check", eventIdx);
	TCanvas *canvas = new TCanvas(canvasName,
                                      Form("Event %d - All Sensors Quality Check", eventIdx),
                                      2400, 800);  // Wide canvas for 4 columns
        canvas->Divide(3, 1);  // 4 columns, 1 row
	
        // Draw each sensor in its own pad (sensors 0-4(maximum sensor number))
        for (int sensorID = 0; sensorID < 4; ++sensorID) {
          canvas->cd(sensorID+1);  // Move to pad sensorID	  
	  auto it = sensorHistograms.find(sensorID);	  
	  if (it != sensorHistograms.end()) {
	    // Clone the histogram to avoid deletion issues
	    TH2F *histClone = (TH2F*)it->second->Clone(Form("sensor%02d_qc_clone", sensorID));
            histClone->SetStats(0);  // Hide statistics box
            histClone->Draw("COLZ");  // Draw with color scale
          } else {
            // Create empty histogram if sensor not found
            TH2F *emptyHist = new TH2F(Form("sensor%02d_empty", sensorID),
                                       Form("Event %d - Sensor %02d (No Data);X;Strip;Amplitude (V)",
                                            eventIdx, sensorID),
                                       1, 0, 1, 8, 0, 8);
            emptyHist->SetStats(0);
            emptyHist->Draw("COLZ");
          }
        }
	
        // Save canvas to quality check file
        canvas->Write(canvasName, TObject::kOverwrite);
        delete canvas;  // Canvas deletion will handle cloned histograms
      }

      // Clean up sensor histograms
      for (auto &pair : sensorHistograms) {
        delete pair.second;
      }

      if (waveformPlotsFile) {
        waveformPlotsFile->cd();
      }
    }

    // Check if waveform plots file needs rotation (after saving all plots for this event)
    if (waveformPlotsFile) {
      checkAndRotateWaveformPlotsFile(waveformPlotsFile);
    }

    outputTree->Fill();
  }

  if (nsamplesError) {
    if (waveformPlotsFile) {
      waveformPlotsFile->cd();
      waveformPlotsFile->Close();
      delete waveformPlotsFile;
      waveformPlotsFile = nullptr;
    }
    if (qualityCheckFile) {
      qualityCheckFile->cd();
      qualityCheckFile->Close();
      delete qualityCheckFile;
      qualityCheckFile = nullptr;
    }
    outputFile->Close();
    inputFile->Close();
    return false;
  }

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
  inputFile->Close();

  // Close waveform plots file if it was created
  if (waveformPlotsFile) {
    std::string finalFileName = waveformPlotsFile->GetName();
    waveformPlotsFile->cd();
    waveformPlotsFile->Close();
    delete waveformPlotsFile;
    std::cout << "Waveform plots output saved to " << finalFileName << std::endl;
    if (waveformPlotsFileCounter > 0) {
      std::cout << "  Total files created: " << (waveformPlotsFileCounter + 1)
                << " (split due to 4GB size limit)" << std::endl;
    }
  }

  // Close quality check file if it was created
  if (qualityCheckFile) {
    std::string finalFileName = qualityCheckFile->GetName();
    qualityCheckFile->cd();
    qualityCheckFile->Close();
    delete qualityCheckFile;
    std::cout << "Quality check output saved to " << finalFileName << std::endl;
  }

  std::string outputFullPath = BuildOutputPath(outname_base, "root", cfg.output_root());
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
  std::string defaultPath = "converter_config.json";
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
      cfg.set_input_root(argv[++i]);
    } else if (arg == "--output") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --output requires a value" << std::endl;
        return 1;
      }
      cfg.set_output_root(argv[++i]);
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
