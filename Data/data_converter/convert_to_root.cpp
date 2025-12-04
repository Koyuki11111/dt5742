#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <thread>
#include <mutex>

#include "TFile.h"
#include "TTree.h"

#include "wave_converter_config.h"

static const size_t HEADER_WORDS = 8;
static const size_t HEADER_BYTES = HEADER_WORDS * sizeof(uint32_t);

namespace {

struct ChannelHeader {
  uint32_t eventSize = 0;
  uint32_t boardId = 0;
  uint32_t channelId = 0;
  uint32_t eventCounter = 0;
};

struct AsciiEventBlock {
  uint32_t boardId = 0;
  uint32_t channelId = 0;
  uint32_t eventCounter = 0;
  int recordLength = 0;
  std::vector<float> samples;
};

struct BinaryEventData {
  uint32_t boardId = 0;
  uint32_t channelId = 0;
  uint32_t eventCounter = 0;
  std::vector<float> samples;
};

std::string TrimCopy(const std::string &text) {
  const auto begin = std::find_if_not(text.begin(), text.end(),
                                      [](unsigned char ch) { return std::isspace(ch); });
  const auto end = std::find_if_not(text.rbegin(), text.rend(),
                                    [](unsigned char ch) { return std::isspace(ch); })
                       .base();
  if (begin >= end) {
    return "";
  }
  return std::string(begin, end);
}

bool TryParseInt(const std::string &text, int &value) {
  try {
    size_t idx = 0;
    value = std::stoi(text, &idx, 0);
    return idx == text.size();
  } catch (...) {
    return false;
  }
}

bool TryParseUint(const std::string &text, uint32_t &value) {
  try {
    size_t idx = 0;
    value = static_cast<uint32_t>(std::stoul(text, &idx, 0));
    return idx == text.size();
  } catch (...) {
    return false;
  }
}

std::string BuildFileName(const WaveConverterConfig &cfg, int ch) {
  const bool canOverride =
      cfg.enable_special_override && !cfg.special_channel_file.empty() &&
      cfg.special_channel_index >= 0 && cfg.special_channel_index < cfg.n_channels;

  if (canOverride && ch == cfg.special_channel_index) {
    if (!cfg.special_channel_file.empty() && cfg.special_channel_file[0] == '/') {
      return cfg.special_channel_file;
    }
    if (!cfg.input_dir.empty()) {
      std::string dir = cfg.input_dir;
      if (dir.back() != '/') {
        dir += '/';
      }
      return dir + cfg.special_channel_file;
    }
    return cfg.special_channel_file;
  }

  char fname[512];
  std::snprintf(fname, sizeof(fname), cfg.input_pattern.c_str(), ch);

  std::string filename(fname);
  if (!cfg.input_dir.empty() && (filename.empty() || filename[0] != '/')) {
    std::string dir = cfg.input_dir;
    if (dir.back() != '/') {
      dir += '/';
    }
    return dir + filename;
  }
  return filename;
}

bool ReadHeader(std::ifstream &fin, ChannelHeader &out) {
  uint32_t header[HEADER_WORDS] = {0};
  if (!fin.read(reinterpret_cast<char *>(header), HEADER_BYTES)) {
    return false;
  }
  out.eventSize = header[0];
  out.boardId = header[1];
  out.channelId = header[3];
  out.eventCounter = header[4];
  return true;
}

bool LoadAsciiChannelFile(const std::string &path,
                          std::vector<AsciiEventBlock> &events) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::cerr << "ERROR: cannot open ASCII input " << path << std::endl;
    return false;
  }

  AsciiEventBlock current;
  bool inSamples = false;

  auto finalizeBlock = [&]() {
    if (current.samples.empty()) {
      current = AsciiEventBlock{};
      inSamples = false;
      return;
    }
    if (current.recordLength == 0) {
      current.recordLength = static_cast<int>(current.samples.size());
    } else if (current.recordLength !=
               static_cast<int>(current.samples.size())) {
      std::cerr << "WARNING: Record Length mismatch in " << path
                << " (expected " << current.recordLength << " got "
                << current.samples.size() << ")" << std::endl;
    }
    events.push_back(current);
    current = AsciiEventBlock{};
    inSamples = false;
  };

  std::string line;
  while (std::getline(fin, line)) {
    std::string trimmed = TrimCopy(line);
    if (trimmed.empty()) {
      continue;
    }

    const size_t colonPos = trimmed.find(':');
    if (colonPos != std::string::npos && !inSamples) {
      const std::string key = TrimCopy(trimmed.substr(0, colonPos));
      const std::string value = TrimCopy(trimmed.substr(colonPos + 1));
      if (key == "Record Length") {
        TryParseInt(value, current.recordLength);
      } else if (key == "BoardID") {
        TryParseUint(value, current.boardId);
      } else if (key == "Channel") {
        TryParseUint(value, current.channelId);
      } else if (key == "Event Number") {
        TryParseUint(value, current.eventCounter);
      }
      continue;
    }

    if (colonPos != std::string::npos && inSamples) {
      finalizeBlock();
      const std::string key = TrimCopy(trimmed.substr(0, colonPos));
      const std::string value = TrimCopy(trimmed.substr(colonPos + 1));
      if (key == "Record Length") {
        TryParseInt(value, current.recordLength);
      } else if (key == "BoardID") {
        TryParseUint(value, current.boardId);
      } else if (key == "Channel") {
        TryParseUint(value, current.channelId);
      } else if (key == "Event Number") {
        TryParseUint(value, current.eventCounter);
      }
      continue;
    }

    inSamples = true;
    try {
      float val = std::stof(trimmed);
      current.samples.push_back(val);
    } catch (const std::exception &) {
      std::cerr << "WARNING: cannot parse sample \"" << trimmed << "\" in "
                << path << std::endl;
    }
  }

  if (!current.samples.empty()) {
    finalizeBlock();
  }

  if (events.empty()) {
    std::cerr << "ERROR: no events parsed from ASCII input " << path
              << std::endl;
    return false;
  }
  return true;
}

bool ReadChannelChunk(std::ifstream &fin, std::mutex &fileMutex,
                      int chunkSize, std::vector<BinaryEventData> &events,
                      uint8_t &eofReached) {
  std::lock_guard<std::mutex> lock(fileMutex);

  events.clear();
  events.reserve(chunkSize);

  for (int i = 0; i < chunkSize; ++i) {
    ChannelHeader header;
    if (!ReadHeader(fin, header)) {
      eofReached = 1;
      break;
    }

    if (header.eventSize <= HEADER_BYTES) {
      std::cerr << "ERROR: invalid event size " << header.eventSize
                << " in ReadChannelChunk" << std::endl;
      return false;
    }

    const uint32_t payloadBytes = header.eventSize - HEADER_BYTES;
    if (payloadBytes % sizeof(float) != 0) {
      std::cerr << "ERROR: payload not multiple of 4 bytes in ReadChannelChunk"
                << std::endl;
      return false;
    }

    const int nsamples = static_cast<int>(payloadBytes / sizeof(float));
    std::vector<float> buffer(nsamples);
    fin.read(reinterpret_cast<char *>(buffer.data()), payloadBytes);
    if (!fin.good() && !fin.eof()) {
      std::cerr << "ERROR: failed to read payload in ReadChannelChunk" << std::endl;
      return false;
    }

    BinaryEventData evt;
    evt.boardId = header.boardId;
    evt.channelId = header.channelId;
    evt.eventCounter = header.eventCounter;
    evt.samples = std::move(buffer);
    events.push_back(std::move(evt));
  }

  return true;
}

bool ConvertBinaryToRoot(const WaveConverterConfig &cfg) {
  // Build output path with subdirectory structure
  std::string outputPath = BuildOutputPath(cfg.output_dir, "root", cfg.root_file);

  // Extract directory path and create it if needed
  size_t lastSlash = outputPath.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = outputPath.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      return false;
    }
  }

  TFile *file = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outputPath << std::endl;
    return false;
  }

  std::cout << "Creating ROOT file: " << outputPath << std::endl;

  TTree *tree = new TTree(cfg.tree_name.c_str(), "Raw waveforms");

  int eventIdx = 0;
  int nChannelsBranch = cfg.n_channels;
  int nsamplesBranch = 0;
  float samplingNs = static_cast<float>(cfg.tsample_ns);
  float pedTarget = static_cast<float>(cfg.ped_target);
  int pedestalWindow = cfg.pedestal_window;

  std::vector<float> timeAxis;
  std::vector<float> pedestals(cfg.n_channels, 0.0f);
  std::vector<uint32_t> boardIds(cfg.n_channels, 0);
  std::vector<uint32_t> channelIds(cfg.n_channels, 0);
  std::vector<uint32_t> eventCounters(cfg.n_channels, 0);
  std::vector<std::vector<float>> raw(cfg.n_channels);
  std::vector<std::vector<float>> ped(cfg.n_channels);

  tree->Branch("event", &eventIdx, "event/I");
  tree->Branch("n_channels", &nChannelsBranch, "n_channels/I");
  tree->Branch("nsamples", &nsamplesBranch, "nsamples/I");
  tree->Branch("sampling_ns", &samplingNs, "sampling_ns/F");
  tree->Branch("ped_target", &pedTarget, "ped_target/F");
  tree->Branch("pedestal_window", &pedestalWindow, "pedestal_window/I");
  tree->Branch("time_ns", &timeAxis);
  tree->Branch("pedestals", &pedestals);
  tree->Branch("board_ids", &boardIds);
  tree->Branch("channel_ids", &channelIds);
  tree->Branch("event_counters", &eventCounters);

  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    char bnameRaw[32];
    char bnamePed[32];
    std::snprintf(bnameRaw, sizeof(bnameRaw), "ch%02d_raw", ch);
    std::snprintf(bnamePed, sizeof(bnamePed), "ch%02d_ped", ch);
    tree->Branch(bnameRaw, &raw[ch]);
    tree->Branch(bnamePed, &ped[ch]);
  }

  std::vector<std::ifstream> fins(cfg.n_channels);
  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    const std::string fname = BuildFileName(cfg, ch);
    fins[ch].open(fname, std::ios::binary);
    if (!fins[ch].is_open()) {
      std::cerr << "ERROR: cannot open " << fname << std::endl;
      file->Close();
      delete file;
      return false;
    }
    std::cout << "Opened " << fname << std::endl;
  }

  bool running = true;
  int eventCount = 0;
  while (running) {
    ChannelHeader chHeader0;
    if (!ReadHeader(fins[0], chHeader0)) {
      std::cout << "EOF reached at event " << eventCount << " (ch0)" << std::endl;
      break;
    }

    if (chHeader0.eventSize <= HEADER_BYTES) {
      std::cerr << "ERROR: invalid event size " << chHeader0.eventSize
                << " at event " << eventCount << " ch0" << std::endl;
      break;
    }

    const uint32_t payloadBytes0 = chHeader0.eventSize - HEADER_BYTES;
    if (payloadBytes0 % sizeof(float) != 0) {
      std::cerr << "ERROR: payload not multiple of 4 bytes at event " << eventCount
                << " ch0" << std::endl;
      break;
    }

    const int nsamples = static_cast<int>(payloadBytes0 / sizeof(float));
    if (timeAxis.size() != static_cast<size_t>(nsamples)) {
      timeAxis.resize(nsamples);
      for (int i = 0; i < nsamples; ++i) {
        timeAxis[i] = static_cast<float>(i * cfg.tsample_ns);
      }
    }

    boardIds[0] = chHeader0.boardId;
    channelIds[0] = chHeader0.channelId;
    eventCounters[0] = chHeader0.eventCounter;
    nsamplesBranch = nsamples;

    std::vector<ChannelHeader> headers(cfg.n_channels);
    headers[0] = chHeader0;

    for (int ch = 1; ch < cfg.n_channels; ++ch) {
      if (!ReadHeader(fins[ch], headers[ch])) {
        std::cerr << "EOF/read error at event " << eventCount << " channel " << ch
                  << std::endl;
        running = false;
        break;
      }
      if (headers[ch].eventSize != chHeader0.eventSize) {
        std::cerr << "WARNING: event size mismatch event " << eventCount << " ch"
                  << ch << " (" << headers[ch].eventSize
                  << " vs " << chHeader0.eventSize << ")" << std::endl;
      }
      boardIds[ch] = headers[ch].boardId;
      channelIds[ch] = headers[ch].channelId;
      eventCounters[ch] = headers[ch].eventCounter;
    }

    if (!running) {
      break;
    }

    for (int ch = 0; ch < cfg.n_channels; ++ch) {
      const uint32_t payloadBytes = headers[ch].eventSize - HEADER_BYTES;
      const int nsampCh = static_cast<int>(payloadBytes / sizeof(float));
      if (nsampCh != nsamples) {
        std::cerr << "WARNING: nsamples mismatch event " << eventCount << " ch"
                  << ch << " (" << nsampCh << " vs " << nsamples << ")"
                  << std::endl;
      }

      std::vector<float> buffer(nsampCh);
      fins[ch].read(reinterpret_cast<char *>(buffer.data()),
                    nsampCh * sizeof(float));
      if (!fins[ch].good()) {
        std::cerr << "ERROR: failed to read payload at event " << eventCount
                  << " ch" << ch << std::endl;
        running = false;
        break;
      }

      raw[ch].assign(buffer.begin(), buffer.end());

      const int pedWindow = std::max(1, cfg.pedestal_window);
      const int nPed = std::min(nsampCh, pedWindow);
      double pedVal = 0.0;
      for (int i = 0; i < nPed; ++i) {
        pedVal += raw[ch][i];
      }
      pedVal /= static_cast<double>(nPed);
      pedestals[ch] = static_cast<float>(pedVal);

      ped[ch].resize(nsampCh);
      for (int i = 0; i < nsampCh; ++i) {
        ped[ch][i] = raw[ch][i] - pedestals[ch] + pedTarget;
      }
    }

    if (!running) {
      break;
    }

    eventIdx = eventCount;
    tree->Fill();
    ++eventCount;
  }

  for (auto &fin : fins) {
    if (fin.is_open()) {
      fin.close();
    }
  }

  if (eventCount == 0) {
    std::cerr << "ERROR: no events converted from binary input." << std::endl;
    file->Close();
    delete file;
    return false;
  }

  file->cd();
  tree->Write();
  file->Close();
  delete file;

  std::cout << "Stage 1: ROOT file written with " << eventCount << " events." << std::endl;
  return true;
}

bool ConvertBinaryToRootParallel(const WaveConverterConfig &cfg) {
  // Build output path with subdirectory structure
  std::string outputPath = BuildOutputPath(cfg.output_dir, "root", cfg.root_file);

  // Extract directory path and create it if needed
  size_t lastSlash = outputPath.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = outputPath.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      return false;
    }
  }

  TFile *file = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outputPath << std::endl;
    return false;
  }

  std::cout << "Creating ROOT file (parallel mode): " << outputPath << std::endl;
  std::cout << "Chunk size: " << cfg.chunk_size << ", Max cores: " << cfg.max_cores << std::endl;

  TTree *tree = new TTree(cfg.tree_name.c_str(), "Raw waveforms");

  int eventIdx = 0;
  int nChannelsBranch = cfg.n_channels;
  int nsamplesBranch = 0;
  float samplingNs = static_cast<float>(cfg.tsample_ns);
  float pedTarget = static_cast<float>(cfg.ped_target);
  int pedestalWindow = cfg.pedestal_window;

  std::vector<float> timeAxis;
  std::vector<float> pedestals(cfg.n_channels, 0.0f);
  std::vector<uint32_t> boardIds(cfg.n_channels, 0);
  std::vector<uint32_t> channelIds(cfg.n_channels, 0);
  std::vector<uint32_t> eventCounters(cfg.n_channels, 0);
  std::vector<std::vector<float>> raw(cfg.n_channels);
  std::vector<std::vector<float>> ped(cfg.n_channels);

  tree->Branch("event", &eventIdx, "event/I");
  tree->Branch("n_channels", &nChannelsBranch, "n_channels/I");
  tree->Branch("nsamples", &nsamplesBranch, "nsamples/I");
  tree->Branch("sampling_ns", &samplingNs, "sampling_ns/F");
  tree->Branch("ped_target", &pedTarget, "ped_target/F");
  tree->Branch("pedestal_window", &pedestalWindow, "pedestal_window/I");
  tree->Branch("time_ns", &timeAxis);
  tree->Branch("pedestals", &pedestals);
  tree->Branch("board_ids", &boardIds);
  tree->Branch("channel_ids", &channelIds);
  tree->Branch("event_counters", &eventCounters);

  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    char bnameRaw[32];
    char bnamePed[32];
    std::snprintf(bnameRaw, sizeof(bnameRaw), "ch%02d_raw", ch);
    std::snprintf(bnamePed, sizeof(bnamePed), "ch%02d_ped", ch);
    tree->Branch(bnameRaw, &raw[ch]);
    tree->Branch(bnamePed, &ped[ch]);
  }

  // Open all channel files
  std::vector<std::ifstream> fins(cfg.n_channels);
  std::vector<std::mutex> fileMutexes(cfg.n_channels);
  std::vector<uint8_t> channelEof(cfg.n_channels, 0);

  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    const std::string fname = BuildFileName(cfg, ch);
    fins[ch].open(fname, std::ios::binary);
    if (!fins[ch].is_open()) {
      std::cerr << "ERROR: cannot open " << fname << std::endl;
      file->Close();
      delete file;
      return false;
    }
    std::cout << "Opened " << fname << std::endl;
  }

  int totalEventsProcessed = 0;
  int chunkNumber = 0;
  const int pedWindow = std::max(1, cfg.pedestal_window);

  bool allEof = false;
  while (!allEof) {
    // Parallel read: each thread reads chunk_size events from its channel
    std::vector<std::thread> threads;
    std::vector<std::vector<BinaryEventData>> chunkData(cfg.n_channels);

    int maxThreads = std::min(cfg.max_cores, cfg.n_channels);
    for (int ch = 0; ch < cfg.n_channels; ch += maxThreads) {
      threads.clear();
      for (int t = 0; t < maxThreads && (ch + t) < cfg.n_channels; ++t) {
        int chIdx = ch + t;
        threads.emplace_back([&, chIdx]() {
          if (!channelEof[chIdx]) {
            ReadChannelChunk(fins[chIdx], fileMutexes[chIdx], cfg.chunk_size,
                           chunkData[chIdx], channelEof[chIdx]);
          }
        });
      }
      for (auto &thread : threads) {
        thread.join();
      }
    }

    // Check if all channels reached EOF
    allEof = true;
    for (int ch = 0; ch < cfg.n_channels; ++ch) {
      if (!channelEof[ch] || !chunkData[ch].empty()) {
        allEof = false;
        break;
      }
    }

    if (allEof) {
      break;
    }

    // Verify all channels read the same number of events
    int nEventsInChunk = static_cast<int>(chunkData[0].size());
    for (int ch = 1; ch < cfg.n_channels; ++ch) {
      if (static_cast<int>(chunkData[ch].size()) != nEventsInChunk) {
        std::cerr << "ERROR: event count mismatch in chunk " << chunkNumber
                  << " ch0=" << nEventsInChunk << " ch" << ch << "="
                  << chunkData[ch].size() << std::endl;
        for (auto &fin : fins) {
          if (fin.is_open()) {
            fin.close();
          }
        }
        file->Close();
        delete file;
        return false;
      }
    }

    if (nEventsInChunk == 0) {
      break;
    }

    // Fill TTree with this chunk
    for (int evt = 0; evt < nEventsInChunk; ++evt) {
      const int nsamples = static_cast<int>(chunkData[0][evt].samples.size());
      nsamplesBranch = nsamples;

      if (timeAxis.size() != static_cast<size_t>(nsamples)) {
        timeAxis.resize(nsamples);
        for (int i = 0; i < nsamples; ++i) {
          timeAxis[i] = static_cast<float>(i * cfg.tsample_ns);
        }
      }

      for (int ch = 0; ch < cfg.n_channels; ++ch) {
        const auto &evtData = chunkData[ch][evt];

        boardIds[ch] = evtData.boardId;
        channelIds[ch] = evtData.channelId;
        eventCounters[ch] = evtData.eventCounter;
        raw[ch] = evtData.samples;

        // Calculate pedestal
        const int nPed = std::min<int>(evtData.samples.size(), pedWindow);
        double pedVal = 0.0;
        for (int i = 0; i < nPed; ++i) {
          pedVal += raw[ch][i];
        }
        pedVal /= static_cast<double>(std::max(1, nPed));
        pedestals[ch] = static_cast<float>(pedVal);

        // Pedestal-subtracted waveform
        ped[ch].resize(evtData.samples.size());
        for (size_t i = 0; i < evtData.samples.size(); ++i) {
          ped[ch][i] = raw[ch][i] - pedestals[ch] + pedTarget;
        }
      }

      eventIdx = totalEventsProcessed;
      tree->Fill();
      ++totalEventsProcessed;
    }

    ++chunkNumber;
    std::cout << "Processed chunk " << chunkNumber << ": " << nEventsInChunk
              << " events (total: " << totalEventsProcessed << ")" << std::endl;
  }

  for (auto &fin : fins) {
    if (fin.is_open()) {
      fin.close();
    }
  }

  if (totalEventsProcessed == 0) {
    std::cerr << "ERROR: no events converted from binary input." << std::endl;
    file->Close();
    delete file;
    return false;
  }

  file->cd();
  tree->Write();
  file->Close();
  delete file;

  std::cout << "Stage 1: ROOT file written with " << totalEventsProcessed
            << " events (parallel mode)." << std::endl;
  return true;
}

bool ConvertAsciiToRoot(const WaveConverterConfig &cfg) {
  // Build output path with subdirectory structure
  std::string outputPath = BuildOutputPath(cfg.output_dir, "root", cfg.root_file);

  // Extract directory path and create it if needed
  size_t lastSlash = outputPath.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = outputPath.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      return false;
    }
  }

  TFile *file = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outputPath << std::endl;
    return false;
  }

  std::cout << "Creating ROOT file: " << outputPath << std::endl;

  TTree *tree = new TTree(cfg.tree_name.c_str(), "Raw waveforms");

  int eventIdx = 0;
  int nChannelsBranch = cfg.n_channels;
  int nsamplesBranch = 0;
  float samplingNs = static_cast<float>(cfg.tsample_ns);
  float pedTarget = static_cast<float>(cfg.ped_target);
  int pedestalWindow = cfg.pedestal_window;

  std::vector<float> timeAxis;
  std::vector<float> pedestals(cfg.n_channels, 0.0f);
  std::vector<uint32_t> boardIds(cfg.n_channels, 0);
  std::vector<uint32_t> channelIds(cfg.n_channels, 0);
  std::vector<uint32_t> eventCounters(cfg.n_channels, 0);
  std::vector<std::vector<float>> raw(cfg.n_channels);
  std::vector<std::vector<float>> ped(cfg.n_channels);

  tree->Branch("event", &eventIdx, "event/I");
  tree->Branch("n_channels", &nChannelsBranch, "n_channels/I");
  tree->Branch("nsamples", &nsamplesBranch, "nsamples/I");
  tree->Branch("sampling_ns", &samplingNs, "sampling_ns/F");
  tree->Branch("ped_target", &pedTarget, "ped_target/F");
  tree->Branch("pedestal_window", &pedestalWindow, "pedestal_window/I");
  tree->Branch("time_ns", &timeAxis);
  tree->Branch("pedestals", &pedestals);
  tree->Branch("board_ids", &boardIds);
  tree->Branch("channel_ids", &channelIds);
  tree->Branch("event_counters", &eventCounters);

  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    char bnameRaw[32];
    char bnamePed[32];
    std::snprintf(bnameRaw, sizeof(bnameRaw), "ch%02d_raw", ch);
    std::snprintf(bnamePed, sizeof(bnamePed), "ch%02d_ped", ch);
    tree->Branch(bnameRaw, &raw[ch]);
    tree->Branch(bnamePed, &ped[ch]);
  }

  std::vector<std::vector<AsciiEventBlock>> channelEvents(cfg.n_channels);
  size_t expectedEvents = 0;
  for (int ch = 0; ch < cfg.n_channels; ++ch) {
    const std::string fname = BuildFileName(cfg, ch);
    if (!LoadAsciiChannelFile(fname, channelEvents[ch])) {
      file->Close();
      delete file;
      return false;
    }
    std::cout << "Loaded ASCII input " << fname << " with "
              << channelEvents[ch].size() << " event(s)." << std::endl;
    if (ch == 0) {
      expectedEvents = channelEvents[ch].size();
    } else if (channelEvents[ch].size() != expectedEvents) {
      std::cerr << "ERROR: ASCII input " << fname << " has "
                << channelEvents[ch].size()
                << " events but expected " << expectedEvents << std::endl;
      file->Close();
      delete file;
      return false;
    }
  }

  if (expectedEvents == 0) {
    std::cerr << "ERROR: no events found in ASCII inputs." << std::endl;
    file->Close();
    delete file;
    return false;
  }

  const int pedWindow = std::max(1, cfg.pedestal_window);
  for (size_t evt = 0; evt < expectedEvents; ++evt) {
    const size_t nsamples = channelEvents[0][evt].samples.size();
    nsamplesBranch = static_cast<int>(nsamples);
    if (timeAxis.size() != nsamples) {
      timeAxis.resize(nsamples);
      for (size_t i = 0; i < nsamples; ++i) {
        timeAxis[i] = static_cast<float>(i * cfg.tsample_ns);
      }
    }

    for (int ch = 0; ch < cfg.n_channels; ++ch) {
      const auto &block = channelEvents[ch][evt];
      if (block.samples.size() != nsamples) {
        std::cerr << "WARNING: nsamples mismatch event " << evt << " ch" << ch
                  << " (" << block.samples.size() << " vs " << nsamples << ")"
                  << std::endl;
      }
      boardIds[ch] = block.boardId;
      channelIds[ch] = block.channelId;
      eventCounters[ch] = block.eventCounter;
      raw[ch] = block.samples;

      const int nPed = std::min<int>(block.samples.size(), pedWindow);
      double pedVal = 0.0;
      for (int i = 0; i < nPed; ++i) {
        pedVal += raw[ch][i];
      }
      pedVal /= static_cast<double>(std::max(1, nPed));
      pedestals[ch] = static_cast<float>(pedVal);

      ped[ch].resize(block.samples.size());
      for (size_t i = 0; i < block.samples.size(); ++i) {
        ped[ch][i] = raw[ch][i] - pedestals[ch] + pedTarget;
      }
    }

    eventIdx = static_cast<int>(evt);
    tree->Fill();
  }

  file->cd();
  tree->Write();
  file->Close();
  delete file;

  std::cout << "Stage 1: ROOT file written with " << expectedEvents
            << " events (ASCII input)." << std::endl;
  return true;
}

void PrintUsage(const char *prog) {
  std::cout << "Stage 1: Convert binary/ASCII waveform files to ROOT format\n"
            << "Usage: " << prog << " [options]\n"
            << "Options:\n"
            << "  --config PATH       Load settings from JSON file\n"
            << "  --pattern PATTERN   Override input filename pattern\n"
            << "  --channels N        Override number of channels\n"
            << "  --root FILE         Override ROOT output file\n"
            << "  --ascii             Read ASCII waveform text files\n"
            << "  --binary            Force binary waveform decoding (default)\n"
            << "  --parallel          Enable parallel loading (binary mode only)\n"
            << "  --chunk-size N      Set chunk size for parallel loading (default: 1000)\n"
            << "  --max-threads N     Set maximum threads for parallel loading\n"
            << "  -h, --help          Show this help message\n";
}

enum class CliOutcome { kOk, kShowUsage, kError };

CliOutcome ApplyCommandLineArgs(int argc, char **argv, WaveConverterConfig &cfg) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto requireValue = [&](const char *name) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: option " << name << " requires a value" << std::endl;
        return nullptr;
      }
      return argv[++i];
    };

    if (arg == "--help" || arg == "-h") {
      return CliOutcome::kShowUsage;
    } else if (arg == "--config") {
      const char *path = requireValue("--config");
      if (!path) {
        return CliOutcome::kError;
      }
      std::string err;
      if (!LoadConfigFromJson(path, cfg, &err)) {
        std::cerr << "ERROR: " << err << std::endl;
        return CliOutcome::kError;
      }
      std::cout << "Loaded configuration from " << path << std::endl;
    } else if (arg == "--pattern") {
      const char *val = requireValue("--pattern");
      if (!val) {
        return CliOutcome::kError;
      }
      cfg.input_pattern = val;
    } else if (arg == "--channels") {
      const char *val = requireValue("--channels");
      if (!val) {
        return CliOutcome::kError;
      }
      try {
        cfg.n_channels = std::stoi(val);
      } catch (const std::exception &) {
        std::cerr << "ERROR: invalid integer for --channels" << std::endl;
        return CliOutcome::kError;
      }
    } else if (arg == "--root") {
      const char *val = requireValue("--root");
      if (!val) {
        return CliOutcome::kError;
      }
      cfg.root_file = val;
    } else if (arg == "--ascii") {
      cfg.input_is_ascii = true;
    } else if (arg == "--binary") {
      cfg.input_is_ascii = false;
    } else if (arg == "--parallel") {
      // Force parallel mode by ensuring max_cores > 1
      if (cfg.max_cores < 2) {
        cfg.max_cores = 2;
      }
    } else if (arg == "--chunk-size") {
      const char *val = requireValue("--chunk-size");
      if (!val) {
        return CliOutcome::kError;
      }
      try {
        cfg.chunk_size = std::stoi(val);
      } catch (const std::exception &) {
        std::cerr << "ERROR: invalid integer for --chunk-size" << std::endl;
        return CliOutcome::kError;
      }
    } else if (arg == "--max-cores" || arg == "--max-threads") {
      const char *val = requireValue(arg.c_str());
      if (!val) {
        return CliOutcome::kError;
      }
      try {
        cfg.max_cores = std::stoi(val);
      } catch (const std::exception &) {
        std::cerr << "ERROR: invalid integer for " << arg << std::endl;
        return CliOutcome::kError;
      }
    } else {
      std::cerr << "ERROR: unknown option " << arg << std::endl;
      return CliOutcome::kError;
    }
  }
  return CliOutcome::kOk;
}

bool LoadDefaultConfig(WaveConverterConfig &cfg) {
  const std::string defaultPath = "pipeline_config.json";
  std::ifstream fin(defaultPath);
  if (!fin.is_open()) {
    return false;
  }
  fin.close();
  std::string err;
  if (!LoadConfigFromJson(defaultPath, cfg, &err)) {
    std::cerr << "WARNING: failed to load default config: " << err << std::endl;
    return false;
  }
  std::cout << "Loaded default configuration from " << defaultPath << std::endl;
  return true;
}

} // namespace

int main(int argc, char **argv) {
  WaveConverterConfig cfg;
  LoadDefaultConfig(cfg);

  CliOutcome cli = ApplyCommandLineArgs(argc, argv, cfg);
  if (cli == CliOutcome::kShowUsage) {
    PrintUsage(argv[0]);
    return 0;
  }
  if (cli == CliOutcome::kError) {
    PrintUsage(argv[0]);
    return 1;
  }

  try {
    bool ok = false;
    if (cfg.input_is_ascii) {
      ok = ConvertAsciiToRoot(cfg);
    } else if (cfg.max_cores > 1) {
      ok = ConvertBinaryToRootParallel(cfg);
    } else {
      ok = ConvertBinaryToRoot(cfg);
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
