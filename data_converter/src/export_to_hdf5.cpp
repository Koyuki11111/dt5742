#include <cstdint>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <regex>

#include "TFile.h"
#include "TTree.h"

#include "utils/filesystem_utils.h"
#include "hdf5.h"
#include "utils/json_utils.h"

namespace {

// Structure to hold DAQ configuration and paths
struct DaqConfig {
  std::string configPath;
  std::string rootFilePath;
  std::string daqName;
  int nChannels;
  std::vector<int> sensorIds;
  std::vector<int> columnIds;
  std::vector<int> stripIds;
};

// Helper function to extract paths and sensor mapping from config JSON
inline bool ExtractDaqConfig(const std::string &configPath,
                              DaqConfig &daqConfig,
                              const std::string &mode,
                              int defaultColumn = 1) {
  simdjson::dom::parser parser;
  simdjson::dom::element root;
  std::string err;
  if (!ParseJsonFile(configPath, parser, root, &err)) {
    std::cerr << "ERROR: " << err << std::endl;
    return false;
  }

  daqConfig.configPath = configPath;

  // Extract common section
  simdjson::dom::element common;
  if (!GetObject(root, "common", common)) {
    std::cerr << "ERROR: common section not found in config" << std::endl;
    return false;
  }

  std::string outputDir, daqName, rootFileName;
  double runnumber = 0;
  double nChannels = 16;

  GetString(common, "output_dir", outputDir);
  GetString(common, "daq_name", daqName);
  GetNumber(common, "runnumber", runnumber);
  GetNumber(common, "n_channels", nChannels);

  daqConfig.nChannels = static_cast<int>(nChannels);

  if (mode == "raw") {
    GetString(common, "waveforms_root", rootFileName);
  } else {
    GetString(common, "analysis_root", rootFileName);
  }

  daqConfig.daqName = daqName;

  // Build ROOT file path: outputDir/runnumber/daqName/output/root/rootFileName
  char runDirBuf[32];
  std::snprintf(runDirBuf, sizeof(runDirBuf), "%06d", static_cast<int>(runnumber));
  daqConfig.rootFilePath = outputDir + "/" + std::string(runDirBuf) + "/" +
                           daqName + "/output/root/" + rootFileName;

  // Extract waveform_analyzer section
  simdjson::dom::element waveformAnalyzer;
  if (!GetObject(root, "waveform_analyzer", waveformAnalyzer)) {
    std::cerr << "ERROR: waveform_analyzer section not found in config" << std::endl;
    return false;
  }

  simdjson::dom::element sensorMapping;
  if (!GetObject(waveformAnalyzer, "sensor_mapping", sensorMapping)) {
    std::cerr << "ERROR: sensor_mapping not found in config" << std::endl;
    return false;
  }

  daqConfig.sensorIds.clear();
  if (!GetIntArray(sensorMapping, "sensor_ids", daqConfig.sensorIds) ||
      daqConfig.sensorIds.empty()) {
    std::cerr << "ERROR: sensor_ids not found in sensor_mapping" << std::endl;
    return false;
  }

  daqConfig.columnIds.clear();
  if (!GetIntArray(sensorMapping, "column_ids", daqConfig.columnIds) ||
      daqConfig.columnIds.empty()) {
    // If not provided, fill with defaultColumn
    daqConfig.columnIds.assign(daqConfig.sensorIds.size(), defaultColumn);
  }

  daqConfig.stripIds.clear();
  if (!GetIntArray(sensorMapping, "strip_ids", daqConfig.stripIds) ||
      daqConfig.stripIds.empty()) {
    // If not provided, fill with channel index (0, 1, 2, ...)
    daqConfig.stripIds.resize(daqConfig.sensorIds.size());
    for (size_t i = 0; i < daqConfig.sensorIds.size(); ++i) {
      daqConfig.stripIds[i] = static_cast<int>(i);
    }
  }

  return true;
}

// Legacy helper for backwards compatibility
inline bool ExtractSensorIds(const std::string &configPath,
                             std::vector<int> &sensorIds,
                             std::vector<int> *columnIds = nullptr,
                             std::vector<int> *stripIds = nullptr,
                             int defaultColumn = 1) {
  DaqConfig daqConfig;
  if (!ExtractDaqConfig(configPath, daqConfig, "analysis", defaultColumn)) {
    return false;
  }

  sensorIds = daqConfig.sensorIds;
  if (columnIds) *columnIds = daqConfig.columnIds;
  if (stripIds) *stripIds = daqConfig.stripIds;
  return true;
}

#pragma pack(push, 1)
struct WaveformMeta {
  uint32_t event;
  uint16_t channel;
  uint16_t nsamples;
  uint32_t board_id;
  uint32_t event_counter;
  float pedestal;
};

struct AnalysisFeatureMeta {
  uint32_t event;
  uint16_t channel;
  uint16_t sensor_id;
  uint16_t column_id;
  uint16_t strip_id;
  float baseline;
  float rmsNoise;
  float noise1Point;
  float ampMinBefore;
  float ampMaxBefore;
  float ampMax;
  float charge;
  float signalOverNoise;
  float peakTime;
  float riseTime;
  float slewRate;
};
#pragma pack(pop)

bool ExportRawWaveforms(const std::string &rootFile,
                       const std::string &treeName,
                       const std::string &hdf5File,
                       int nChannels,
                       int sensorFilter = -1,
                       const std::vector<int> *sensorIds = nullptr) {
  TFile *fin = TFile::Open(rootFile.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open ROOT file " << rootFile << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(fin->Get(treeName.c_str()));
  if (!tree) {
    std::cerr << "ERROR: tree " << treeName << " not found" << std::endl;
    fin->Close();
    return false;
  }

  int eventIdx = 0;
  int nsamples = 0;
  int nChannelsBranch = 0;
  float samplingNs = 0.0f;
  float pedTarget = 0.0f;

  std::vector<float> *timeAxis = nullptr;
  std::vector<float> *pedestals = nullptr;
  std::vector<uint32_t> *boardIds = nullptr;
  std::vector<uint32_t> *eventCounters = nullptr;
  std::vector<int> *nsamplesPerChannel = nullptr;

  tree->SetBranchAddress("event", &eventIdx);
  tree->SetBranchAddress("nsamples", &nsamples);
  tree->SetBranchAddress("n_channels", &nChannelsBranch);
  tree->SetBranchAddress("sampling_ns", &samplingNs);
  tree->SetBranchAddress("ped_target", &pedTarget);
  tree->SetBranchAddress("time_ns", &timeAxis);
  tree->SetBranchAddress("pedestals", &pedestals);
  tree->SetBranchAddress("board_ids", &boardIds);
  tree->SetBranchAddress("event_counters", &eventCounters);
  if (tree->GetBranch("nsamples_per_channel")) {
    tree->SetBranchAddress("nsamples_per_channel", &nsamplesPerChannel);
  }

  const int maxChannels = nChannels;
  std::vector<std::vector<float> *> chPedPtrs(maxChannels, nullptr);

  for (int ch = 0; ch < maxChannels; ++ch) {
    char bname[32];
    std::snprintf(bname, sizeof(bname), "ch%02d_ped", ch);
    if (tree->GetBranch(bname)) {
      tree->SetBranchAddress(bname, &chPedPtrs[ch]);
    }
  }

  const Long64_t nEntries = tree->GetEntries();
  if (nEntries <= 0) {
    std::cerr << "WARNING: tree contains no entries, skipping HDF5 export"
              << std::endl;
    fin->Close();
    return false;
  }

  std::vector<WaveformMeta> metadata;
  metadata.reserve(static_cast<size_t>(nEntries) * maxChannels);

  std::vector<std::vector<float>> waveformRows;
  waveformRows.reserve(static_cast<size_t>(nEntries) * maxChannels);
  std::vector<float> rowPadValues;
  rowPadValues.reserve(static_cast<size_t>(nEntries) * maxChannels);
  size_t maxSamplesPerRow = 0;
  bool loggedNsamplesTrim = false;

  std::vector<float> timeAxisCopy;

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);

    if (!timeAxisCopy.size() && timeAxis) {
      timeAxisCopy.assign(timeAxis->begin(), timeAxis->end());
    }

    if (!pedestals || !boardIds || !eventCounters) {
      std::cerr << "ERROR: missing per-channel vectors in tree entry "
                << entry << std::endl;
      fin->Close();
      return false;
    }

    for (int ch = 0; ch < maxChannels; ++ch) {
      // Filter by sensor if requested
      if (sensorFilter >= 0 && sensorIds && ch < static_cast<int>(sensorIds->size())) {
        if ((*sensorIds)[ch] != sensorFilter) {
          continue;  // Skip this channel, it's not from the requested sensor
        }
      }

      auto *vecPtr = chPedPtrs[ch];
      if (!vecPtr) {
        continue;
      }

      int chSamples = nsamples;
      if (nsamplesPerChannel &&
          ch < static_cast<int>(nsamplesPerChannel->size())) {
        chSamples = nsamplesPerChannel->at(ch);
      }

      chSamples =
          std::min(chSamples, static_cast<int>(vecPtr->size()));
      if (timeAxis) {
        chSamples =
            std::min(chSamples, static_cast<int>(timeAxis->size()));
      }

      if (chSamples <= 0) {
        continue;
      }

      if (!loggedNsamplesTrim &&
          chSamples != static_cast<int>(vecPtr->size())) {
        std::cout << "INFO: trimming waveform samples at entry " << entry
                  << " ch" << ch << " to " << chSamples
                  << " for HDF5 export" << std::endl;
        loggedNsamplesTrim = true;
      }

      WaveformMeta meta{};
      meta.event = static_cast<uint32_t>(eventIdx);
      meta.channel = static_cast<uint16_t>(ch);
      meta.nsamples = static_cast<uint16_t>(chSamples);
      meta.board_id =
          (ch < static_cast<int>(boardIds->size())) ? (*boardIds)[ch] : 0;
      meta.event_counter =
          (ch < static_cast<int>(eventCounters->size()))
              ? (*eventCounters)[ch]
              : 0;
      meta.pedestal =
          (ch < static_cast<int>(pedestals->size())) ? (*pedestals)[ch] : 0.0f;

      metadata.push_back(meta);
      maxSamplesPerRow =
          std::max(maxSamplesPerRow, static_cast<size_t>(chSamples));

      std::vector<float> row(vecPtr->begin(),
                             vecPtr->begin() + chSamples);
      waveformRows.push_back(std::move(row));
      rowPadValues.push_back(pedTarget);
    }
  }

  fin->Close();

  if (metadata.empty()) {
    std::cerr << "WARNING: no waveform metadata filled, aborting HDF5 export"
              << std::endl;
    return false;
  }

  if (maxSamplesPerRow == 0) {
    std::cerr << "ERROR: no waveform samples found to export" << std::endl;
    return false;
  }

  if (waveformRows.size() != metadata.size()) {
    std::cerr << "ERROR: internal mismatch between metadata rows and waveforms ("
              << metadata.size() << " vs " << waveformRows.size() << ")"
              << std::endl;
    return false;
  }

  const size_t rows = metadata.size();
  const size_t samplesPerRow = maxSamplesPerRow;

  std::vector<float> waveforms;
  waveforms.reserve(rows * samplesPerRow);
  for (size_t i = 0; i < waveformRows.size(); ++i) {
    const auto &row = waveformRows[i];
    const float padValue =
        (i < rowPadValues.size()) ? rowPadValues[i] : pedTarget;
    waveforms.insert(waveforms.end(), row.begin(), row.end());
    if (row.size() < samplesPerRow) {
      waveforms.insert(waveforms.end(), samplesPerRow - row.size(),
                       padValue);
    }
  }

  if (waveforms.size() != rows * samplesPerRow) {
    std::cerr << "ERROR: waveform buffer size mismatch (" << waveforms.size()
              << " vs " << rows * samplesPerRow << ")" << std::endl;
    return false;
  }

  hid_t file =
      H5Fcreate(hdf5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    std::cerr << "ERROR: cannot create HDF5 file " << hdf5File << std::endl;
    return false;
  }

  // Metadata dataset
  hsize_t metaDim = rows;
  hid_t metaSpace = H5Screate_simple(1, &metaDim, nullptr);
  hid_t metaType = H5Tcreate(H5T_COMPOUND, sizeof(WaveformMeta));
  H5Tinsert(metaType, "event", HOFFSET(WaveformMeta, event), H5T_NATIVE_UINT32);
  H5Tinsert(metaType, "channel", HOFFSET(WaveformMeta, channel),
            H5T_NATIVE_UINT16);
  H5Tinsert(metaType, "nsamples", HOFFSET(WaveformMeta, nsamples),
            H5T_NATIVE_UINT16);
  H5Tinsert(metaType, "board_id", HOFFSET(WaveformMeta, board_id),
            H5T_NATIVE_UINT32);
  H5Tinsert(metaType, "event_counter",
            HOFFSET(WaveformMeta, event_counter), H5T_NATIVE_UINT32);
  H5Tinsert(metaType, "pedestal", HOFFSET(WaveformMeta, pedestal),
            H5T_NATIVE_FLOAT);

  hid_t metaSet =
      H5Dcreate(file, "Metadata", metaType, metaSpace, H5P_DEFAULT, H5P_DEFAULT,
                H5P_DEFAULT);
  if (metaSet < 0) {
    std::cerr << "ERROR: cannot create Metadata dataset" << std::endl;
    H5Tclose(metaType);
    H5Sclose(metaSpace);
    H5Fclose(file);
    return false;
  }

  H5Dwrite(metaSet, metaType, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata.data());

  // Waveform dataset (rows x samples)
  hsize_t waveDims[2] = {rows, samplesPerRow};
  hid_t waveSpace = H5Screate_simple(2, waveDims, nullptr);
  hid_t waveSet =
      H5Dcreate(file, "Waveforms", H5T_NATIVE_FLOAT, waveSpace, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
  if (waveSet < 0) {
    std::cerr << "ERROR: cannot create Waveforms dataset" << std::endl;
    H5Dclose(metaSet);
    H5Tclose(metaType);
    H5Sclose(metaSpace);
    H5Sclose(waveSpace);
    H5Fclose(file);
    return false;
  }

  H5Dwrite(waveSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           waveforms.data());

  // Time axis dataset
  if (!timeAxisCopy.empty()) {
    hsize_t timeDim = timeAxisCopy.size();
    hid_t timeSpace = H5Screate_simple(1, &timeDim, nullptr);
    hid_t timeSet =
        H5Dcreate(file, "TimeAxis_ns", H5T_NATIVE_FLOAT, timeSpace, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
    if (timeSet >= 0) {
      H5Dwrite(timeSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               timeAxisCopy.data());
      H5Dclose(timeSet);
    } else {
      std::cerr << "WARNING: failed to create TimeAxis_ns dataset" << std::endl;
    }
    H5Sclose(timeSpace);
  }

  // File-level attributes
  hid_t attrSpace = H5Screate(H5S_SCALAR);
  if (attrSpace >= 0) {
    hid_t attrSampling =
        H5Acreate2(file, "sampling_ns", H5T_NATIVE_FLOAT, attrSpace,
                   H5P_DEFAULT, H5P_DEFAULT);
    if (attrSampling >= 0) {
      H5Awrite(attrSampling, H5T_NATIVE_FLOAT, &samplingNs);
      H5Aclose(attrSampling);
    }

    hid_t attrPedTarget =
        H5Acreate2(file, "ped_target", H5T_NATIVE_FLOAT, attrSpace,
                   H5P_DEFAULT, H5P_DEFAULT);
    if (attrPedTarget >= 0) {
      H5Awrite(attrPedTarget, H5T_NATIVE_FLOAT, &pedTarget);
      H5Aclose(attrPedTarget);
    }
    H5Sclose(attrSpace);
  }

  H5Dclose(waveSet);
  H5Sclose(waveSpace);
  H5Dclose(metaSet);
  H5Tclose(metaType);
  H5Sclose(metaSpace);
  H5Fclose(file);

  std::cout << "HDF5 raw waveforms written to " << hdf5File << std::endl;
  return true;
}

bool ExportAnalysisFeatures(const std::string &rootFile,
                            const std::string &treeName,
                            const std::string &hdf5File,
                            int nChannels,
                            int sensorFilter = -1,
                            const std::vector<int> *sensorIds = nullptr,
                            const std::vector<int> *columnIds = nullptr,
                            const std::vector<int> *stripIds = nullptr,
                            bool append = false) {
  TFile *fin = TFile::Open(rootFile.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open ROOT file " << rootFile << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(fin->Get(treeName.c_str()));
  if (!tree) {
    std::cerr << "ERROR: tree " << treeName << " not found" << std::endl;
    fin->Close();
    return false;
  }

  int event = 0;
  std::vector<float> *baseline = nullptr;
  std::vector<float> *rmsNoise = nullptr;
  std::vector<float> *noise1Point = nullptr;
  std::vector<float> *ampMinBefore = nullptr;
  std::vector<float> *ampMaxBefore = nullptr;
  std::vector<float> *ampMax = nullptr;
  std::vector<float> *charge = nullptr;
  std::vector<float> *signalOverNoise = nullptr;
  std::vector<float> *peakTime = nullptr;
  std::vector<float> *riseTime = nullptr;
  std::vector<float> *slewRate = nullptr;

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("baseline", &baseline);
  tree->SetBranchAddress("rmsNoise", &rmsNoise);
  tree->SetBranchAddress("noise1Point", &noise1Point);
  tree->SetBranchAddress("ampMinBefore", &ampMinBefore);
  tree->SetBranchAddress("ampMaxBefore", &ampMaxBefore);
  tree->SetBranchAddress("ampMax", &ampMax);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("signalOverNoise", &signalOverNoise);
  tree->SetBranchAddress("peakTime", &peakTime);
  tree->SetBranchAddress("riseTime", &riseTime);
  tree->SetBranchAddress("slewRate", &slewRate);

  const Long64_t nEntries = tree->GetEntries();
  if (nEntries <= 0) {
    std::cerr << "WARNING: tree contains no entries" << std::endl;
    fin->Close();
    return false;
  }

  std::vector<AnalysisFeatureMeta> features;
  features.reserve(static_cast<size_t>(nEntries) * nChannels);

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);

    if (!baseline || !ampMax) {
      continue;
    }

    for (int ch = 0; ch < nChannels; ++ch) {
      // Filter by sensor if requested
      if (sensorFilter >= 0 && sensorIds && ch < static_cast<int>(sensorIds->size())) {
        if ((*sensorIds)[ch] != sensorFilter) {
          continue;  // Skip this channel, it's not from the requested sensor
        }
      }

      AnalysisFeatureMeta meta{};
      meta.event = static_cast<uint32_t>(event);
      meta.channel = static_cast<uint16_t>(ch);

      // Set sensor mapping info
      if (sensorIds && ch < static_cast<int>(sensorIds->size())) {
        meta.sensor_id = static_cast<uint16_t>((*sensorIds)[ch]);
      } else {
        meta.sensor_id = 0;
      }

      if (columnIds && ch < static_cast<int>(columnIds->size())) {
        meta.column_id = static_cast<uint16_t>((*columnIds)[ch]);
      } else {
        meta.column_id = 1;
      }

      if (stripIds && ch < static_cast<int>(stripIds->size())) {
        meta.strip_id = static_cast<uint16_t>((*stripIds)[ch]);
      } else {
        meta.strip_id = static_cast<uint16_t>(ch);
      }

      meta.baseline = (ch < static_cast<int>(baseline->size())) ? (*baseline)[ch] : 0.0f;
      meta.rmsNoise = (ch < static_cast<int>(rmsNoise->size())) ? (*rmsNoise)[ch] : 0.0f;
      meta.noise1Point = (ch < static_cast<int>(noise1Point->size())) ? (*noise1Point)[ch] : 0.0f;
      meta.ampMinBefore = (ch < static_cast<int>(ampMinBefore->size())) ? (*ampMinBefore)[ch] : 0.0f;
      meta.ampMaxBefore = (ch < static_cast<int>(ampMaxBefore->size())) ? (*ampMaxBefore)[ch] : 0.0f;
      meta.ampMax = (ch < static_cast<int>(ampMax->size())) ? (*ampMax)[ch] : 0.0f;
      meta.charge = (ch < static_cast<int>(charge->size())) ? (*charge)[ch] : 0.0f;
      meta.signalOverNoise = (ch < static_cast<int>(signalOverNoise->size())) ? (*signalOverNoise)[ch] : 0.0f;
      meta.peakTime = (ch < static_cast<int>(peakTime->size())) ? (*peakTime)[ch] : 0.0f;
      meta.riseTime = (ch < static_cast<int>(riseTime->size())) ? (*riseTime)[ch] : 0.0f;
      meta.slewRate = (ch < static_cast<int>(slewRate->size())) ? (*slewRate)[ch] : 0.0f;

      features.push_back(meta);
    }
  }

  fin->Close();

  if (features.empty()) {
    std::cerr << "WARNING: no features extracted" << std::endl;
    return false;
  }

  hid_t file = -1;
  if (append) {
    file = H5Fopen(hdf5File.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
      std::cerr << "ERROR: cannot open HDF5 file for appending " << hdf5File << std::endl;
      return false;
    }
  } else {
    file = H5Fcreate(hdf5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  if (file < 0) {
    std::cerr << "ERROR: cannot create HDF5 file " << hdf5File << std::endl;
    return false;
  }

  hsize_t dim = features.size();
  hid_t space = H5Screate_simple(1, &dim, nullptr);
  hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(AnalysisFeatureMeta));

  H5Tinsert(type, "event", HOFFSET(AnalysisFeatureMeta, event), H5T_NATIVE_UINT32);
  H5Tinsert(type, "channel", HOFFSET(AnalysisFeatureMeta, channel), H5T_NATIVE_UINT16);
  H5Tinsert(type, "sensor_id", HOFFSET(AnalysisFeatureMeta, sensor_id), H5T_NATIVE_UINT16);
  H5Tinsert(type, "column_id", HOFFSET(AnalysisFeatureMeta, column_id), H5T_NATIVE_UINT16);
  H5Tinsert(type, "strip_id", HOFFSET(AnalysisFeatureMeta, strip_id), H5T_NATIVE_UINT16);
  H5Tinsert(type, "baseline", HOFFSET(AnalysisFeatureMeta, baseline), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "rmsNoise", HOFFSET(AnalysisFeatureMeta, rmsNoise), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "noise1Point", HOFFSET(AnalysisFeatureMeta, noise1Point), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "ampMinBefore", HOFFSET(AnalysisFeatureMeta, ampMinBefore), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "ampMaxBefore", HOFFSET(AnalysisFeatureMeta, ampMaxBefore), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "ampMax", HOFFSET(AnalysisFeatureMeta, ampMax), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "charge", HOFFSET(AnalysisFeatureMeta, charge), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "signalOverNoise", HOFFSET(AnalysisFeatureMeta, signalOverNoise), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "peakTime", HOFFSET(AnalysisFeatureMeta, peakTime), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "riseTime", HOFFSET(AnalysisFeatureMeta, riseTime), H5T_NATIVE_FLOAT);
  H5Tinsert(type, "slewRate", HOFFSET(AnalysisFeatureMeta, slewRate), H5T_NATIVE_FLOAT);

  hid_t dset = H5Dcreate(file, "AnalysisFeatures", type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    std::cerr << "ERROR: cannot create dataset" << std::endl;
    H5Tclose(type);
    H5Sclose(space);
    H5Fclose(file);
    return false;
  }

  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, features.data());

  H5Dclose(dset);
  H5Tclose(type);
  H5Sclose(space);
  H5Fclose(file);

  std::cout << "HDF5 analysis features written to " << hdf5File << std::endl;
  return true;
}

bool ExportCorryHits(const std::string &rootFile,
                     const std::string &treeName,
                     const std::string &hdf5File,
                     int nChannels,
                     int sensorFilter = -1,
                     const std::vector<int> *sensorIds = nullptr,
                     const std::vector<int> *columnIds = nullptr,
                     const std::vector<int> *stripIds = nullptr,
                     int defaultColumn = 1,
                     bool onlyCorryFields = true) {
  TFile *fin = TFile::Open(rootFile.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open ROOT file " << rootFile << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(fin->Get(treeName.c_str()));
  if (!tree) {
    std::cerr << "ERROR: tree " << treeName << " not found" << std::endl;
    fin->Close();
    return false;
  }

  int event = 0;
  std::vector<float> *charge = nullptr;
  std::vector<float> *peakTime = nullptr;

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("peakTime", &peakTime);

  const Long64_t nEntries = tree->GetEntries();
  if (nEntries <= 0) {
    std::cerr << "WARNING: tree contains no entries" << std::endl;
    fin->Close();
    return false;
  }

  #pragma pack(push, 1)
  struct HitRow {
    uint16_t column;
    uint16_t row;
    uint8_t raw;
    double charge;
    double timestamp;
    uint32_t trigger_number;
  };
  #pragma pack(pop)

  std::vector<HitRow> hits;
  hits.reserve(static_cast<size_t>(nEntries) * nChannels);

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);

    for (int ch = 0; ch < nChannels; ++ch) {
      if (sensorFilter >= 0 && sensorIds && ch < static_cast<int>(sensorIds->size())) {
        if ((*sensorIds)[ch] != sensorFilter) {
          continue;
        }
      }

      HitRow hit{};
      // Column: default or per-channel mapping if provided
      if (columnIds && ch < static_cast<int>(columnIds->size())) {
        hit.column = static_cast<uint16_t>((*columnIds)[ch]);
      } else {
        hit.column = static_cast<uint16_t>(defaultColumn);
      }
      // Row: use strip_ids if available, otherwise use channel index
      if (stripIds && ch < static_cast<int>(stripIds->size())) {
        hit.row = static_cast<uint16_t>((*stripIds)[ch]);
      } else {
        hit.row = static_cast<uint16_t>(ch);
      }
      hit.raw = 0u;
      hit.charge =
          (charge && ch < static_cast<int>(charge->size())) ? static_cast<double>((*charge)[ch]) : 0.0;
      hit.timestamp =
          (peakTime && ch < static_cast<int>(peakTime->size())) ? static_cast<double>((*peakTime)[ch]) : 0.0;
      hit.trigger_number = static_cast<uint32_t>(event);

      hits.push_back(hit);
    }
  }

  fin->Close();

  if (hits.empty()) {
    std::cerr << "WARNING: no hits extracted for Corryvreckan format" << std::endl;
    return false;
  }

  hid_t file =
      H5Fcreate(hdf5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    std::cerr << "ERROR: cannot create HDF5 file " << hdf5File << std::endl;
    return false;
  }

  hsize_t dim = hits.size();
  hid_t space = H5Screate_simple(1, &dim, nullptr);

  hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(HitRow));
  H5Tinsert(type, "column", HOFFSET(HitRow, column), H5T_NATIVE_UINT16);
  H5Tinsert(type, "row", HOFFSET(HitRow, row), H5T_NATIVE_UINT16);
  H5Tinsert(type, "raw", HOFFSET(HitRow, raw), H5T_NATIVE_UINT8);
  H5Tinsert(type, "charge", HOFFSET(HitRow, charge), H5T_NATIVE_DOUBLE);
  H5Tinsert(type, "timestamp", HOFFSET(HitRow, timestamp), H5T_NATIVE_DOUBLE);
  H5Tinsert(type, "trigger_number", HOFFSET(HitRow, trigger_number), H5T_NATIVE_UINT32);

  hid_t dset = H5Dcreate(file, "Hits", type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    std::cerr << "ERROR: cannot create Hits dataset" << std::endl;
    H5Tclose(type);
    H5Sclose(space);
    H5Fclose(file);
    return false;
  }

  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, hits.data());

  // Mark whether only Corryvreckan fields are stored
  hid_t attrSpace = H5Screate(H5S_SCALAR);
  if (attrSpace >= 0) {
    unsigned char flag = onlyCorryFields ? 1 : 0;
    hid_t attr = H5Acreate2(file, "corry_only_fields", H5T_NATIVE_UCHAR, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
    if (attr >= 0) {
      H5Awrite(attr, H5T_NATIVE_UCHAR, &flag);
      H5Aclose(attr);
    }
    H5Sclose(attrSpace);
  }

  H5Dclose(dset);
  H5Tclose(type);
  H5Sclose(space);
  H5Fclose(file);

  std::cout << "HDF5 Corryvreckan Hits written to " << hdf5File << std::endl;
  return true;
}

// Export analysis features from multiple DAQ configs, merging data by sensor
bool ExportAnalysisFeaturesMultiDAQ(const std::vector<DaqConfig> &daqConfigs,
                                     const std::string &treeName,
                                     const std::string &outputDir,
                                     const std::string &baseOutputName,
                                     bool splitBySensor = true) {
  if (daqConfigs.empty()) {
    std::cerr << "ERROR: no DAQ configs provided" << std::endl;
    return false;
  }

  // Collect all unique sensor IDs across all DAQs
  std::set<int> uniqueSensorIds;
  for (const auto &cfg : daqConfigs) {
    for (int sid : cfg.sensorIds) {
      uniqueSensorIds.insert(sid);
    }
  }

  std::cout << "Found " << uniqueSensorIds.size() << " unique sensors across "
            << daqConfigs.size() << " DAQs" << std::endl;

  // Process each sensor
  for (int sensorId : uniqueSensorIds) {
    std::cout << "\nProcessing sensor " << sensorId << "..." << std::endl;

    #pragma pack(push, 1)
    struct HitRow {
      uint16_t column;
      uint16_t row;
      uint8_t raw;
      double charge;
      double timestamp;
      uint32_t trigger_number;
    };
    #pragma pack(pop)

    std::vector<HitRow> allHits;

    // Collect data from all DAQs for this sensor
    for (const auto &daqCfg : daqConfigs) {
      std::cout << "  Reading " << daqCfg.daqName << ": " << daqCfg.rootFilePath << std::endl;

      TFile *fin = TFile::Open(daqCfg.rootFilePath.c_str(), "READ");
      if (!fin || fin->IsZombie()) {
        std::cerr << "  WARNING: cannot open ROOT file " << daqCfg.rootFilePath
                  << ", skipping" << std::endl;
        continue;
      }

      TTree *tree = dynamic_cast<TTree *>(fin->Get(treeName.c_str()));
      if (!tree) {
        std::cerr << "  WARNING: tree " << treeName << " not found, skipping" << std::endl;
        fin->Close();
        continue;
      }

      int event = 0;
      std::vector<float> *ampMax = nullptr;
      std::vector<float> *peakTime = nullptr;

      tree->SetBranchAddress("event", &event);
      tree->SetBranchAddress("ampMax", &ampMax);
      tree->SetBranchAddress("peakTime", &peakTime);

      const Long64_t nEntries = tree->GetEntries();
      size_t channelsAdded = 0;

      for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        if (!ampMax) {
          continue;
        }

        for (int ch = 0; ch < daqCfg.nChannels; ++ch) {
          // Filter by sensor
          if (ch >= static_cast<int>(daqCfg.sensorIds.size())) {
            continue;
          }
          if (daqCfg.sensorIds[ch] != sensorId) {
            continue;  // This channel belongs to a different sensor
          }

          HitRow hit{};

          // Column: from columnIds mapping
          if (ch < static_cast<int>(daqCfg.columnIds.size())) {
            hit.column = static_cast<uint16_t>(daqCfg.columnIds[ch]);
          } else {
            hit.column = 1;
          }

          // Row: from stripIds mapping
          if (ch < static_cast<int>(daqCfg.stripIds.size())) {
            hit.row = static_cast<uint16_t>(daqCfg.stripIds[ch]);
          } else {
            hit.row = static_cast<uint16_t>(ch);
          }

          hit.raw = 0u;
          // Use ampMax as a proxy for charge (ADC units)
          hit.charge = (ch < static_cast<int>(ampMax->size())) ? static_cast<double>((*ampMax)[ch]) : 0.0;
          hit.timestamp = (ch < static_cast<int>(peakTime->size())) ? static_cast<double>((*peakTime)[ch]) : 0.0;
          hit.trigger_number = static_cast<uint32_t>(event);

          allHits.push_back(hit);
          if (entry == 0) channelsAdded++;
        }
      }

      fin->Close();
      std::cout << "  Added " << channelsAdded << " channels from " << daqCfg.daqName << std::endl;
    }

    // Sort hits by trigger_number (event), then by column, then by row
    // This ensures Corryvreckan reads data event-by-event with all columns interleaved
    std::cout << "  Sorting " << allHits.size() << " hits by event, column, row..." << std::endl;
    std::sort(allHits.begin(), allHits.end(), [](const HitRow &a, const HitRow &b) {
      if (a.trigger_number != b.trigger_number) return a.trigger_number < b.trigger_number;
      if (a.column != b.column) return a.column < b.column;
      return a.row < b.row;
    });

    if (allHits.empty()) {
      std::cerr << "  WARNING: no hits for sensor " << sensorId << ", skipping" << std::endl;
      continue;
    }

    // Create output HDF5 file for this sensor
    std::string hdf5File;
    if (splitBySensor) {
      hdf5File = outputDir + "/sensor" + std::to_string(sensorId) + "_" + baseOutputName;
    } else {
      hdf5File = outputDir + "/" + baseOutputName;
    }

    // Create output directory if needed
    size_t lastSlash = hdf5File.find_last_of('/');
    if (lastSlash != std::string::npos) {
      std::string dirPath = hdf5File.substr(0, lastSlash);
      if (!CreateDirectoryIfNeeded(dirPath)) {
        std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
        return false;
      }
    }

    hid_t file = H5Fcreate(hdf5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
      std::cerr << "ERROR: cannot create HDF5 file " << hdf5File << std::endl;
      return false;
    }

    hsize_t dim = allHits.size();
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(HitRow));

    H5Tinsert(type, "column", HOFFSET(HitRow, column), H5T_NATIVE_UINT16);
    H5Tinsert(type, "row", HOFFSET(HitRow, row), H5T_NATIVE_UINT16);
    H5Tinsert(type, "raw", HOFFSET(HitRow, raw), H5T_NATIVE_UINT8);
    H5Tinsert(type, "charge", HOFFSET(HitRow, charge), H5T_NATIVE_DOUBLE);
    H5Tinsert(type, "timestamp", HOFFSET(HitRow, timestamp), H5T_NATIVE_DOUBLE);
    H5Tinsert(type, "trigger_number", HOFFSET(HitRow, trigger_number), H5T_NATIVE_UINT32);

    hid_t dset = H5Dcreate(file, "Hits", type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
      std::cerr << "ERROR: cannot create Hits dataset" << std::endl;
      H5Tclose(type);
      H5Sclose(space);
      H5Fclose(file);
      return false;
    }

    H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, allHits.data());

    H5Dclose(dset);
    H5Tclose(type);
    H5Sclose(space);
    H5Fclose(file);

    std::cout << "  Wrote " << allHits.size() << " hits to " << hdf5File << std::endl;
  }

  std::cout << "\nMulti-DAQ export completed successfully" << std::endl;
  return true;
}

void PrintUsage(const char *prog) {
  std::cout << "Export ROOT data to HDF5 format\n"
            << "Usage: " << prog << " [options]\n"
            << "\n"
            << "=== Multi-DAQ Mode (NEW) ===\n"
            << "  --config FILE       DAQ config file (can be specified multiple times)\n"
            << "  --mode MODE         Export mode: 'analysis' only for multi-DAQ\n"
            << "  --tree NAME         Input tree name (default: 'Analysis')\n"
            << "  --output-dir DIR    Output directory for HDF5 files (required)\n"
            << "  --output-name NAME  Base output filename (default: 'merged_analysis.h5')\n"
            << "  --split-by-sensor   Split output by sensor (default: true)\n"
            << "\n"
            << "=== Single-DAQ Mode (Legacy) ===\n"
            << "  --mode MODE         Export mode: 'raw', 'analysis', or 'corry' (required)\n"
            << "  --input FILE        Input ROOT file (required)\n"
            << "  --tree NAME         Input tree name (required)\n"
            << "  --output FILE       Output HDF5 file (required)\n"
            << "  --channels N        Number of channels (default: 16)\n"
            << "  --output-dir DIR    Output directory (default: 'output')\n"
            << "  --sensor-id ID      Export only channels from this sensor ID\n"
            << "  --sensor-mapping FILE  Load sensor mapping from analysis config JSON\n"
            << "  --use-sensor-mapping BOOL  Enable/disable applying mapping (default: true)\n"
            << "  --corry-only-fields BOOL   If true, store only fields used by Corryvreckan (default: true)\n"
            << "  --column-id ID      Default column value for corry mode (default: 1)\n"
            << "\n"
            << "=== Common Options ===\n"
            << "  -h, --help          Show this help message\n"
            << "\n"
            << "=== Examples ===\n"
            << "Multi-DAQ merging by sensor:\n"
            << "  " << prog << " --config converter_config_daq00.json --config converter_config_daq01.json \\\n"
            << "            --mode analysis --output-dir /data/000139/merged/hdf5\n"
            << "\n"
            << "Single-DAQ legacy mode:\n"
            << "  " << prog << " --mode analysis --input waveforms_analyzed.root --tree Analysis \\\n"
            << "            --output output.h5 --channels 16\n";
}

} // namespace

int main(int argc, char **argv) {
  // Multi-DAQ mode variables
  std::vector<std::string> configFiles;
  std::string outputName = "merged_analysis.h5";
  bool splitBySensor = true;

  // Single-DAQ mode variables (legacy)
  std::string mode;
  std::string inputRoot;
  std::string treeName = "Analysis";  // Default tree name
  std::string outputHdf5;
  std::string outputDir;
  std::string sensorMappingFile;
  bool useSensorMapping = true;
  bool corryOnlyFields = true;
  int nChannels = 16;
  int sensorFilter = -1;  // -1 means no filtering
  int defaultColumnId = 1;
  std::vector<int> sensorIds;
  std::vector<int> columnIds;
  std::vector<int> stripIds;

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
      configFiles.push_back(argv[++i]);
    } else if (arg == "--output-name") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --output-name requires a value" << std::endl;
        return 1;
      }
      outputName = argv[++i];
    } else if (arg == "--split-by-sensor") {
      splitBySensor = true;
    } else if (arg == "--mode") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --mode requires a value" << std::endl;
        return 1;
      }
      mode = argv[++i];
    } else if (arg == "--input") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --input requires a value" << std::endl;
        return 1;
      }
      inputRoot = argv[++i];
    } else if (arg == "--tree") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --tree requires a value" << std::endl;
        return 1;
      }
      treeName = argv[++i];
    } else if (arg == "--output") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --output requires a value" << std::endl;
        return 1;
      }
      outputHdf5 = argv[++i];
    } else if (arg == "--output-dir") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --output-dir requires a value" << std::endl;
        return 1;
      }
      outputDir = argv[++i];
    } else if (arg == "--channels") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --channels requires a value" << std::endl;
        return 1;
      }
      try {
        nChannels = std::stoi(argv[++i]);
      } catch (...) {
        std::cerr << "ERROR: invalid number for --channels" << std::endl;
        return 1;
      }
    } else if (arg == "--sensor-id") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --sensor-id requires a value" << std::endl;
        return 1;
      }
      try {
        sensorFilter = std::stoi(argv[++i]);
      } catch (...) {
        std::cerr << "ERROR: invalid number for --sensor-id" << std::endl;
        return 1;
      }
    } else if (arg == "--sensor-mapping") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --sensor-mapping requires a value" << std::endl;
        return 1;
      }
      sensorMappingFile = argv[++i];
    } else if (arg == "--column-id") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --column-id requires a value" << std::endl;
        return 1;
      }
      try {
        defaultColumnId = std::stoi(argv[++i]);
      } catch (...) {
        std::cerr << "ERROR: invalid number for --column-id" << std::endl;
        return 1;
      }
    } else if (arg == "--use-sensor-mapping") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --use-sensor-mapping requires a value (true/false)" << std::endl;
        return 1;
      }
      std::string val = argv[++i];
      std::transform(val.begin(), val.end(), val.begin(), ::tolower);
      if (val == "true" || val == "1" || val == "yes") {
        useSensorMapping = true;
      } else if (val == "false" || val == "0" || val == "no") {
        useSensorMapping = false;
      } else {
        std::cerr << "ERROR: invalid value for --use-sensor-mapping (use true/false)" << std::endl;
        return 1;
      }
    } else if (arg == "--corry-only-fields") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --corry-only-fields requires a value (true/false)" << std::endl;
        return 1;
      }
      std::string val = argv[++i];
      std::transform(val.begin(), val.end(), val.begin(), ::tolower);
      if (val == "true" || val == "1" || val == "yes") {
        corryOnlyFields = true;
      } else if (val == "false" || val == "0" || val == "no") {
        corryOnlyFields = false;
      } else {
        std::cerr << "ERROR: invalid value for --corry-only-fields (use true/false)" << std::endl;
        return 1;
      }
    } else {
      std::cerr << "ERROR: unknown option " << arg << std::endl;
      PrintUsage(argv[0]);
      return 1;
    }
  }

  // Multi-DAQ mode: if --config is specified, use the new multi-DAQ logic
  if (!configFiles.empty()) {
    std::cout << "=== Multi-DAQ Mode ===" << std::endl;
    std::cout << "Config files: " << configFiles.size() << std::endl;
    for (const auto &cfg : configFiles) {
      std::cout << "  - " << cfg << std::endl;
    }

    if (mode.empty()) {
      mode = "analysis";  // Default to analysis mode
    }
    if (mode != "analysis") {
      std::cerr << "ERROR: multi-DAQ mode only supports 'analysis' mode" << std::endl;
      return 1;
    }
    if (outputDir.empty()) {
      std::cerr << "ERROR: --output-dir is required for multi-DAQ mode" << std::endl;
      return 1;
    }

    // Load all DAQ configs
    std::vector<DaqConfig> daqConfigs;
    for (const auto &cfgFile : configFiles) {
      DaqConfig daqCfg;
      if (!ExtractDaqConfig(cfgFile, daqCfg, mode, defaultColumnId)) {
        std::cerr << "ERROR: failed to load config from " << cfgFile << std::endl;
        return 1;
      }
      daqConfigs.push_back(daqCfg);
      std::cout << "Loaded " << daqCfg.daqName << " from " << cfgFile << std::endl;
      std::cout << "  ROOT file: " << daqCfg.rootFilePath << std::endl;
      std::cout << "  Channels: " << daqCfg.nChannels << std::endl;
    }

    // Run multi-DAQ export
    bool ok = ExportAnalysisFeaturesMultiDAQ(daqConfigs, treeName, outputDir, outputName, splitBySensor);
    return ok ? 0 : 2;
  }

  // Single-DAQ mode (legacy)
  std::cout << "=== Single-DAQ Mode (Legacy) ===" << std::endl;

  if (mode.empty() || inputRoot.empty() || treeName.empty() || outputHdf5.empty()) {
    std::cerr << "ERROR: missing required arguments for single-DAQ mode" << std::endl;
    PrintUsage(argv[0]);
    return 1;
  }

  if (outputDir.empty()) {
    outputDir = "output";  // Default for legacy mode
  }

  // Load sensor mapping if requested
  const std::vector<int> *sensorIdsPtr = nullptr;
  const std::vector<int> *columnIdsPtr = nullptr;
  const std::vector<int> *stripIdsPtr = nullptr;
  if (sensorFilter >= 0) {
    if (!useSensorMapping) {
      std::cerr << "ERROR: --sensor-id requires mapping, but --use-sensor-mapping=false" << std::endl;
      return 1;
    }
    if (sensorMappingFile.empty()) {
      std::cerr << "ERROR: --sensor-id requires --sensor-mapping" << std::endl;
      return 1;
    }
    if (!ExtractSensorIds(sensorMappingFile, sensorIds, &columnIds, &stripIds, defaultColumnId)) {
      std::cerr << "ERROR: failed to load sensor IDs from " << sensorMappingFile << std::endl;
      return 1;
    }
    sensorIdsPtr = &sensorIds;
    columnIdsPtr = &columnIds;
    stripIdsPtr = &stripIds;
    std::cout << "Filtering for sensor ID " << sensorFilter << std::endl;
  } else if (!sensorMappingFile.empty() && useSensorMapping) {
    // Allow mapping without filtering
    if (!ExtractSensorIds(sensorMappingFile, sensorIds, &columnIds, &stripIds, defaultColumnId)) {
      std::cerr << "ERROR: failed to load sensor IDs from " << sensorMappingFile << std::endl;
      return 1;
    }
    sensorIdsPtr = &sensorIds;
    columnIdsPtr = &columnIds;
    stripIdsPtr = &stripIds;
  }

  // Build full paths with directory structure
  std::string inputPath = BuildPath(outputDir, "root", inputRoot);
  std::string outputPath = BuildPath(outputDir, "hdf5", outputHdf5);

  // Create output directory if needed
  size_t lastSlash = outputPath.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = outputPath.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      return 1;
    }
  }

  try {
    bool ok = false;
    if (mode == "raw") {
      ok = ExportRawWaveforms(inputPath, treeName, outputPath, nChannels, sensorFilter, sensorIdsPtr);
    } else if (mode == "analysis") {
      ok = ExportAnalysisFeatures(inputPath, treeName, outputPath, nChannels, sensorFilter, sensorIdsPtr, columnIdsPtr, stripIdsPtr);
    } else if (mode == "corry") {
      ok = ExportCorryHits(inputPath,
                           treeName,
                           outputPath,
                           nChannels,
                           sensorFilter,
                           sensorIdsPtr,
                           columnIdsPtr,
                           stripIdsPtr,
                           defaultColumnId,
                           corryOnlyFields);
      if (ok && !corryOnlyFields) {
        // Append analysis features for richer files if requested
        bool appended = ExportAnalysisFeatures(
            inputPath, treeName, outputPath, nChannels, sensorFilter, sensorIdsPtr, columnIdsPtr, stripIdsPtr, true /*append*/);
        if (!appended) {
          std::cerr << "ERROR: failed to append AnalysisFeatures dataset" << std::endl;
          return 1;
        }
      }
    } else {
      std::cerr << "ERROR: unknown mode '" << mode << "'. Use 'raw', 'analysis', or 'corry'" << std::endl;
      return 1;
    }

    if (!ok) {
      return 2;
    }
  } catch (const std::exception &ex) {
    std::cerr << "Unhandled exception: " << ex.what() << std::endl;
    return 3;
  }

  return 0;
}
