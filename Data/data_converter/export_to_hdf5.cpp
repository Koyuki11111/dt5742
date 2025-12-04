#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <regex>
#include <sys/stat.h>
#include <sys/types.h>

#include "TFile.h"
#include "TTree.h"

#include "hdf5.h"

namespace {

// Helper function to create directory recursively
inline bool CreateDirectoryIfNeeded(const std::string &path) {
  if (path.empty() || path == ".") {
    return true;  // Current directory always exists
  }

  struct stat info;
  if (stat(path.c_str(), &info) == 0) {
    if (info.st_mode & S_IFDIR) {
      return true;  // Directory already exists
    } else {
      std::cerr << "ERROR: Path exists but is not a directory: " << path << std::endl;
      return false;
    }
  }

  // Try to create directory
  // For nested paths, we need to create parent directories first
  size_t pos = path.find_last_of('/');
  if (pos != std::string::npos && pos > 0) {
    std::string parent = path.substr(0, pos);
    if (!CreateDirectoryIfNeeded(parent)) {
      return false;
    }
  }

  // Create the directory (mkdir with 0755 permissions)
  if (mkdir(path.c_str(), 0755) != 0) {
    std::cerr << "ERROR: Failed to create directory: " << path << std::endl;
    return false;
  }

  return true;
}

// Helper function to build path with subdirectory
inline std::string BuildPath(const std::string &output_dir,
                             const std::string &subdir,
                             const std::string &filename) {
  // If filename is absolute path, use it as-is
  if (!filename.empty() && filename[0] == '/') {
    return filename;
  }

  // If output_dir is empty or ".", just use filename
  if (output_dir.empty() || output_dir == ".") {
    return filename;
  }

  // Build path: output_dir/subdir/filename
  std::string path = output_dir;
  if (path.back() != '/') {
    path += '/';
  }
  path += subdir;
  if (path.back() != '/') {
    path += '/';
  }
  path += filename;

  return path;
}

// Helper function to extract sensor_ids from analysis config JSON
inline bool ExtractSensorIds(const std::string &configPath, std::vector<int> &sensorIds) {
  std::ifstream fin(configPath);
  if (!fin.is_open()) {
    std::cerr << "ERROR: cannot open config file: " << configPath << std::endl;
    return false;
  }

  std::string content((std::istreambuf_iterator<char>(fin)),
                      std::istreambuf_iterator<char>());
  fin.close();

  // Extract sensor_mapping block
  std::string sensorPattern = "\"sensor_mapping\"\\s*:\\s*\\{([^}]*)\\}";
  std::regex sensorRe(sensorPattern);
  std::smatch sensorMatch;
  if (!std::regex_search(content, sensorMatch, sensorRe)) {
    std::cerr << "ERROR: sensor_mapping not found in config" << std::endl;
    return false;
  }

  std::string sensorContent = sensorMatch[1].str();

  // Extract sensor_ids array
  std::string arrayPattern = "\"sensor_ids\"\\s*:\\s*\\[([^\\]]*)\\]";
  std::regex arrayRe(arrayPattern);
  std::smatch arrayMatch;
  if (!std::regex_search(sensorContent, arrayMatch, arrayRe)) {
    std::cerr << "ERROR: sensor_ids not found in sensor_mapping" << std::endl;
    return false;
  }

  std::string arrayContent = arrayMatch[1].str();
  sensorIds.clear();

  // Parse integers
  std::regex numPattern("(-?\\d+)");
  auto numBegin = std::sregex_iterator(arrayContent.begin(), arrayContent.end(), numPattern);
  auto numEnd = std::sregex_iterator();
  for (auto it = numBegin; it != numEnd; ++it) {
    sensorIds.push_back(std::stoi(it->str()));
  }

  return !sensorIds.empty();
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

  tree->SetBranchAddress("event", &eventIdx);
  tree->SetBranchAddress("nsamples", &nsamples);
  tree->SetBranchAddress("n_channels", &nChannelsBranch);
  tree->SetBranchAddress("sampling_ns", &samplingNs);
  tree->SetBranchAddress("ped_target", &pedTarget);
  tree->SetBranchAddress("time_ns", &timeAxis);
  tree->SetBranchAddress("pedestals", &pedestals);
  tree->SetBranchAddress("board_ids", &boardIds);
  tree->SetBranchAddress("event_counters", &eventCounters);

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

  std::vector<float> waveforms;
  waveforms.reserve(static_cast<size_t>(nEntries) * maxChannels * 1024);

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

      if (static_cast<int>(vecPtr->size()) != nsamples) {
        std::cerr << "WARNING: nsamples mismatch entry " << entry << " ch" << ch
                  << " branch size " << vecPtr->size()
                  << " header " << nsamples << std::endl;
      }

      WaveformMeta meta{};
      meta.event = static_cast<uint32_t>(eventIdx);
      meta.channel = static_cast<uint16_t>(ch);
      meta.nsamples = static_cast<uint16_t>(vecPtr->size());
      meta.board_id =
          (ch < static_cast<int>(boardIds->size())) ? (*boardIds)[ch] : 0;
      meta.event_counter =
          (ch < static_cast<int>(eventCounters->size()))
              ? (*eventCounters)[ch]
              : 0;
      meta.pedestal =
          (ch < static_cast<int>(pedestals->size())) ? (*pedestals)[ch] : 0.0f;

      metadata.push_back(meta);
      waveforms.insert(waveforms.end(), vecPtr->begin(), vecPtr->end());
    }
  }

  fin->Close();

  if (metadata.empty()) {
    std::cerr << "WARNING: no waveform metadata filled, aborting HDF5 export"
              << std::endl;
    return false;
  }

  const size_t rows = metadata.size();
  const size_t samplesPerRow =
      metadata.front().nsamples > 0 ? metadata.front().nsamples : nsamples;

  if (rows * samplesPerRow != waveforms.size()) {
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

  hid_t file =
      H5Fcreate(hdf5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    std::cerr << "ERROR: cannot create HDF5 file " << hdf5File << std::endl;
    return false;
  }

  hsize_t dim = features.size();
  hid_t space = H5Screate_simple(1, &dim, nullptr);
  hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(AnalysisFeatureMeta));

  H5Tinsert(type, "event", HOFFSET(AnalysisFeatureMeta, event), H5T_NATIVE_UINT32);
  H5Tinsert(type, "channel", HOFFSET(AnalysisFeatureMeta, channel), H5T_NATIVE_UINT16);
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

void PrintUsage(const char *prog) {
  std::cout << "Export ROOT data to HDF5 format\n"
            << "Usage: " << prog << " [options]\n"
            << "Options:\n"
            << "  --mode MODE         Export mode: 'raw' or 'analysis' (required)\n"
            << "  --input FILE        Input ROOT file (required)\n"
            << "  --tree NAME         Input tree name (required)\n"
            << "  --output FILE       Output HDF5 file (required)\n"
            << "  --channels N        Number of channels (default: 16)\n"
            << "  --output-dir DIR    Output directory (default: 'output')\n"
            << "  --sensor-id ID      Export only channels from this sensor ID\n"
            << "  --sensor-mapping FILE  Load sensor mapping from analysis config JSON\n"
            << "  -h, --help          Show this help message\n"
            << "\n"
            << "Note: If output-dir is specified, files will be organized:\n"
            << "  - Input ROOT: output-dir/root/input\n"
            << "  - Output HDF5: output-dir/hdf5/output\n"
            << "\n"
            << "Sensor filtering:\n"
            << "  --sensor-id filters channels by sensor. Requires --sensor-mapping\n"
            << "  to load sensor_ids from analysis_config.json\n";
}

} // namespace

int main(int argc, char **argv) {
  std::string mode;
  std::string inputRoot;
  std::string treeName;
  std::string outputHdf5;
  std::string outputDir = "output";  // Default output directory
  std::string sensorMappingFile;
  int nChannels = 16;
  int sensorFilter = -1;  // -1 means no filtering
  std::vector<int> sensorIds;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      return 0;
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
    } else {
      std::cerr << "ERROR: unknown option " << arg << std::endl;
      PrintUsage(argv[0]);
      return 1;
    }
  }

  if (mode.empty() || inputRoot.empty() || treeName.empty() || outputHdf5.empty()) {
    std::cerr << "ERROR: missing required arguments" << std::endl;
    PrintUsage(argv[0]);
    return 1;
  }

  // Load sensor mapping if requested
  const std::vector<int> *sensorIdsPtr = nullptr;
  if (sensorFilter >= 0) {
    if (sensorMappingFile.empty()) {
      std::cerr << "ERROR: --sensor-id requires --sensor-mapping" << std::endl;
      return 1;
    }
    if (!ExtractSensorIds(sensorMappingFile, sensorIds)) {
      std::cerr << "ERROR: failed to load sensor IDs from " << sensorMappingFile << std::endl;
      return 1;
    }
    sensorIdsPtr = &sensorIds;
    std::cout << "Filtering for sensor ID " << sensorFilter << std::endl;
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
      ok = ExportAnalysisFeatures(inputPath, treeName, outputPath, nChannels, sensorFilter, sensorIdsPtr);
    } else {
      std::cerr << "ERROR: unknown mode '" << mode << "'. Use 'raw' or 'analysis'" << std::endl;
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
