#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <thread>
#include <mutex>

#include "TFile.h"
#include "TTree.h"

#include "config/wave_converter_config.h"
#include "utils/file_io.h"

namespace {

enum class NsamplesPolicy { kStrict, kPad };
enum class EventPolicy { kError, kWarn, kSkip };

bool kSetEventLimit = false;

void CheckEventLimit(const WaveConverterConfig &cfg){
  if ( cfg.max_events() < 0) {
    kSetEventLimit = false;
  } else {
    kSetEventLimit = true;
  }
}
  
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

EventPolicy ResolveEventPolicy(const std::string &policyText) {
  std::string lowered = policyText;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lowered == "warn") {
    return EventPolicy::kWarn;
  }
  if (lowered == "skip") {
    return EventPolicy::kSkip;
  }
  if (lowered != "error") {
    std::cerr << "WARNING: unknown event_policy '" << policyText
              << "', defaulting to 'error'" << std::endl;
  }
  return EventPolicy::kError;
}

std::string BuildFileName(const WaveConverterConfig &cfg, int ch) {
  const bool canOverride =
      cfg.enable_special_override && !cfg.special_channel_file.empty() &&
      cfg.special_channel_index >= 0 && cfg.special_channel_index < cfg.n_channels();

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

bool IsSpecialOverrideChannel(const WaveConverterConfig &cfg, int ch) {
  return cfg.enable_special_override && cfg.special_channel_index == ch &&
         cfg.special_channel_index >= 0 && cfg.special_channel_index < cfg.n_channels();
}

bool EnsureParentDirectory(const std::string &path) {
  size_t lastSlash = path.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = path.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      return false;
    }
  }
  return true;
}

bool ConvertBinaryToRoot(const WaveConverterConfig &cfg) {
  std::string outputPath = BuildOutputPath(cfg.output_dir(), "root", cfg.root_file());
  if (!EnsureParentDirectory(outputPath)) {
    return false;
  }

  TFile *file = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outputPath << std::endl;
    return false;
  }

  std::cout << "Creating ROOT file: " << outputPath << std::endl;

  TTree *tree = new TTree(cfg.tree_name().c_str(), "Raw waveforms");

  int eventIdx = 0;
  int nChannelsBranch = cfg.n_channels();
  int nsamplesBranch = 0;
  float samplingNs = static_cast<float>(cfg.tsample_ns);
  float pedTarget = static_cast<float>(cfg.ped_target);
  int pedestalWindow = cfg.pedestal_window;
  const NsamplesPolicy policy = ResolveNsamplesPolicy(cfg.common.nsamples_policy);
  const EventPolicy eventPolicy = ResolveEventPolicy(cfg.event_policy);

  std::vector<float> timeAxis;
  std::vector<float> pedestals(cfg.n_channels(), 0.0f);
  std::vector<uint32_t> boardIds(cfg.n_channels(), 0);
  std::vector<uint32_t> channelIds(cfg.n_channels(), 0);
  std::vector<uint32_t> eventCounters(cfg.n_channels(), 0);
  std::vector<int> nsamplesPerChannel(cfg.n_channels(), 0);
  std::vector<std::vector<float>> readBuffers(cfg.n_channels());
  std::vector<std::vector<float>> raw(cfg.n_channels());
  std::vector<std::vector<float>> ped(cfg.n_channels());
  bool loggedSpecialChannelIdInfo = false;

  auto defineCommonBranches = [&]() {
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
    tree->Branch("nsamples_per_channel", &nsamplesPerChannel);
  };

  auto defineChannelBranches = [&]() {
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      char bnameRaw[32];
      char bnamePed[32];
      std::snprintf(bnameRaw, sizeof(bnameRaw), "ch%02d_raw", ch);
      std::snprintf(bnamePed, sizeof(bnamePed), "ch%02d_ped", ch);
      tree->Branch(bnameRaw, &raw[ch]);
      tree->Branch(bnamePed, &ped[ch]);
    }
  };

  defineCommonBranches();
  defineChannelBranches();

  std::vector<std::ifstream> fins(cfg.n_channels());
  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
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
  bool encounteredError = false;
  bool loggedNsamplesPadding = false;
  int consistencyWarningCount = 0;
  const int kConsistencyWarnLimit = 20;
  int eventCount = 0;
  std::vector<bool> channelEof(cfg.n_channels(), false);

  while (running) {
    // Try to read headers from all channels
    std::vector<ChannelHeader> headers(cfg.n_channels());
    std::vector<int> samplesThisEvent(cfg.n_channels(), 0);
    bool anyEof = false;

    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (!ReadHeader(fins[ch], headers[ch])) {
        channelEof[ch] = true;
        anyEof = true;
      }
    }

    // If any channel reached EOF, report detailed status and stop
    if (anyEof) {
      std::cout << "INFO: Event count mismatch detected - one or more channels reached EOF" << std::endl;
      std::cout << "      Per-channel status at event " << eventCount << ":" << std::endl;
      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        std::cout << "        ch" << ch << ": "
                  << (channelEof[ch] ? "reached EOF" : "has more data") << std::endl;
      }
      std::cout << "      Processing stopped. Total events processed: " << eventCount << std::endl;
      break;
    }

    // Validate ch0 header
    if (headers[0].eventSize <= HEADER_BYTES) {
      std::cerr << "ERROR: invalid event size " << headers[0].eventSize
                << " at event " << eventCount << " ch0" << std::endl;
      encounteredError = true;
      break;
    }

    const uint32_t payloadBytes0 = headers[0].eventSize - HEADER_BYTES;
    if (payloadBytes0 % sizeof(float) != 0) {
      std::cerr << "ERROR: payload not multiple of 4 bytes at event " << eventCount
                << " ch0" << std::endl;
      encounteredError = true;
      break;
    }

    samplesThisEvent[0] = static_cast<int>(payloadBytes0 / sizeof(float));
    boardIds[0] = headers[0].boardId;
    channelIds[0] = headers[0].channelId;
    eventCounters[0] = headers[0].eventCounter;

    // Validate and process other channels
    for (int ch = 1; ch < cfg.n_channels(); ++ch) {
      const uint32_t payloadBytes = headers[ch].eventSize - HEADER_BYTES;
      if (headers[ch].eventSize != headers[0].eventSize) {
        std::cerr << "WARNING: event size mismatch event " << eventCount << " ch"
                  << ch << " (" << headers[ch].eventSize
                  << " vs " << headers[0].eventSize << ")" << std::endl;
      }
      if (payloadBytes % sizeof(float) != 0) {
        std::cerr << "ERROR: payload not multiple of 4 bytes at event " << eventCount
                  << " ch" << ch << std::endl;
        running = false;
        encounteredError = true;
        break;
      }
      samplesThisEvent[ch] = static_cast<int>(payloadBytes / sizeof(float));
      boardIds[ch] = headers[ch].boardId;
      channelIds[ch] = headers[ch].channelId;
      eventCounters[ch] = headers[ch].eventCounter;
    }

    if (!running) {
      break;
    }

    // Consistency checks across channels
    std::vector<std::string> consistencyIssues;
    for (int ch = 1; ch < cfg.n_channels(); ++ch) {
      if (headers[ch].eventCounter != headers[0].eventCounter) {
        std::ostringstream oss;
        oss << "eventCounter mismatch at event " << eventCount << " ch" << ch
            << " (" << headers[ch].eventCounter << " vs " << headers[0].eventCounter << ")";
        consistencyIssues.push_back(oss.str());
      }
      if (headers[ch].boardId != headers[0].boardId) {
        std::ostringstream oss;
        oss << "boardId mismatch at event " << eventCount << " ch" << ch
            << " (" << headers[ch].boardId << " vs " << headers[0].boardId << ")";
        consistencyIssues.push_back(oss.str());
      }
      if (headers[ch].channelId != static_cast<uint32_t>(ch)) {
        if (IsSpecialOverrideChannel(cfg, ch)) {
          if (!loggedSpecialChannelIdInfo) {
            std::cout << "INFO: special_channel_index " << ch
                      << " allows channelId mismatch (header "
                      << headers[ch].channelId << ")" << std::endl;
            loggedSpecialChannelIdInfo = true;
          }
        } else {
          std::ostringstream oss;
          oss << "channelId mismatch at event " << eventCount << " ch" << ch
              << " (" << headers[ch].channelId << " vs expected " << ch << ")";
          consistencyIssues.push_back(oss.str());
        }
      }
    }

    bool skipEvent = false;
    if (!consistencyIssues.empty()) {
      for (const auto &msg : consistencyIssues) {
        if (eventPolicy == EventPolicy::kWarn &&
            consistencyWarningCount >= kConsistencyWarnLimit) {
          continue;
        }
        std::cerr << ((eventPolicy == EventPolicy::kWarn) ? "WARNING: " : "ERROR: ")
                  << msg << std::endl;
        if (eventPolicy == EventPolicy::kWarn) {
          ++consistencyWarningCount;
          if (consistencyWarningCount == kConsistencyWarnLimit) {
            std::cerr << "WARNING: further consistency warnings suppressed"
                      << std::endl;
          }
        }
      }
      if (eventPolicy == EventPolicy::kError) {
        encounteredError = true;
        break;
      } else if (eventPolicy == EventPolicy::kSkip) {
        skipEvent = true;
      }
    }

    const int maxSamples =
        *std::max_element(samplesThisEvent.begin(), samplesThisEvent.end());
    const bool hasMismatch =
        std::any_of(samplesThisEvent.begin(), samplesThisEvent.end(),
                    [maxSamples](int v) { return v != maxSamples; });

    if (hasMismatch && policy == NsamplesPolicy::kStrict) {
      std::cerr << "ERROR: nsamples mismatch at event " << eventCount
                << ", expected uniform sample counts across channels" << std::endl;
      encounteredError = true;
      break;
    }

    if (hasMismatch && !loggedNsamplesPadding && policy == NsamplesPolicy::kPad) {
      std::cout << "INFO: nsamples mismatch detected at event " << eventCount
                << ", padding shorter channels up to " << maxSamples
                << " samples" << std::endl;
      loggedNsamplesPadding = true;
    }

    if (timeAxis.size() != static_cast<size_t>(maxSamples)) {
      timeAxis.resize(maxSamples);
      for (int i = 0; i < maxSamples; ++i) {
        timeAxis[i] = static_cast<float>(i * cfg.tsample_ns);
      }
    }

    nsamplesBranch = maxSamples;
    nsamplesPerChannel = samplesThisEvent;

    // Read payloads from all channels
    std::vector<bool> readFailed(cfg.n_channels(), false);
    bool anyReadFailed = false;

    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      const int nsampCh = samplesThisEvent[ch];

      auto &buffer = readBuffers[ch];
      buffer.resize(nsampCh);
      fins[ch].read(reinterpret_cast<char *>(buffer.data()),
                    nsampCh * sizeof(float));
      if (!fins[ch].good()) {
        readFailed[ch] = true;
        anyReadFailed = true;
      } else {
        raw[ch] = std::move(buffer);
        buffer.clear();
      }
    }

    // If any channel failed to read payload, report and stop
    if (anyReadFailed) {
      std::cout << "INFO: Payload read failure - one or more channels encountered early EOF" << std::endl;
      std::cout << "      Per-channel status at event " << eventCount << ":" << std::endl;
      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        std::cout << "        ch" << ch << ": "
                  << (readFailed[ch] ? "read failed (EOF)" : "read successful") << std::endl;
      }
      std::cout << "      Processing stopped. Total events processed: " << eventCount << std::endl;
      break;
    }

    // Calculate pedestals and pedestal-subtracted waveforms
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      const int nsampCh = samplesThisEvent[ch];
      const int pedWindow = std::max(1, cfg.pedestal_window);
      const int nPed = std::min(nsampCh, pedWindow);
      double pedVal = 0.0;
      for (int i = 0; i < nPed; ++i) {
        pedVal += raw[ch][i];
      }
      pedVal /= static_cast<double>(nPed);
      pedestals[ch] = static_cast<float>(pedVal);

      raw[ch].resize(maxSamples, pedestals[ch]);
      ped[ch].resize(maxSamples);
      for (int i = 0; i < nsampCh; ++i) {
        ped[ch][i] = raw[ch][i] - pedestals[ch] + pedTarget;
      }
      for (int i = nsampCh; i < maxSamples; ++i) {
        ped[ch][i] = pedTarget;
      }
    }

    eventIdx = eventCount;
    if (!skipEvent) {
      tree->Fill();
    }
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

  if (encounteredError) {
    std::cerr << "ERROR: conversion stopped due to earlier errors." << std::endl;
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
  std::string outputPath = BuildOutputPath(cfg.output_dir(), "root", cfg.root_file());
  if (!EnsureParentDirectory(outputPath)) {
    return false;
  }

  TFile *file = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outputPath << std::endl;
    return false;
  }

  std::cout << "Creating ROOT file (parallel mode): " << outputPath << std::endl;
  std::cout << "Chunk size: " << cfg.chunk_size() << ", Max cores: " << cfg.max_cores() << std::endl;

  CheckEventLimit(cfg);
  if ( kSetEventLimit ) {
    std::cout<<"The number of analyze event size is set = "<<cfg.max_events()<<std::endl;
  } else {
    std::cout<<"All events will be analyzed"<<std::endl;
  }
  
  TTree *tree = new TTree(cfg.tree_name().c_str(), "Raw waveforms");

  int eventIdx = 0;
  int nChannelsBranch = cfg.n_channels();
  int nsamplesBranch = 0;
  float samplingNs = static_cast<float>(cfg.tsample_ns);
  float pedTarget = static_cast<float>(cfg.ped_target);
  int pedestalWindow = cfg.pedestal_window;
  const NsamplesPolicy policy = ResolveNsamplesPolicy(cfg.common.nsamples_policy);
  const EventPolicy eventPolicy = ResolveEventPolicy(cfg.event_policy);

  std::vector<float> timeAxis;
  std::vector<float> pedestals(cfg.n_channels(), 0.0f);
  std::vector<uint32_t> boardIds(cfg.n_channels(), 0);
  std::vector<uint32_t> channelIds(cfg.n_channels(), 0);
  std::vector<uint32_t> eventCounters(cfg.n_channels(), 0);
  std::vector<int> nsamplesPerChannel(cfg.n_channels(), 0);
  std::vector<std::vector<float>> raw(cfg.n_channels());
  std::vector<std::vector<float>> ped(cfg.n_channels());
  bool loggedSpecialChannelIdInfo = false;

  auto defineCommonBranches = [&]() {
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
    tree->Branch("nsamples_per_channel", &nsamplesPerChannel);
  };

  auto defineChannelBranches = [&]() {
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      char bnameRaw[32];
      char bnamePed[32];
      std::snprintf(bnameRaw, sizeof(bnameRaw), "ch%02d_raw", ch);
      std::snprintf(bnamePed, sizeof(bnamePed), "ch%02d_ped", ch);
      tree->Branch(bnameRaw, &raw[ch]);
      tree->Branch(bnamePed, &ped[ch]);
    }
  };

  defineCommonBranches();
  defineChannelBranches();

  // Open all channel files
  std::vector<std::ifstream> fins(cfg.n_channels());
  std::vector<std::mutex> fileMutexes(cfg.n_channels());
  std::vector<uint8_t> channelEof(cfg.n_channels(), 0);

  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
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
  bool loggedNsamplesPadding = false;
  bool loggedEventPolicyInfo = false;
  int consistencyWarningCount = 0;
  const int kConsistencyWarnLimit = 20;

  bool allEof = false;
  while (!allEof) {
    // Parallel read: each thread reads chunk_size events from its channel
    std::vector<std::thread> threads;
    std::vector<std::vector<BinaryEventData>> chunkData(cfg.n_channels());

    int maxThreads = std::min(cfg.max_cores(), cfg.n_channels());
    for (int ch = 0; ch < cfg.n_channels(); ch += maxThreads) {
      threads.clear();
      for (int t = 0; t < maxThreads && (ch + t) < cfg.n_channels(); ++t) {
        int chIdx = ch + t;
        threads.emplace_back([&, chIdx]() {
          if (!channelEof[chIdx]) {
            ReadChannelChunk(fins[chIdx], fileMutexes[chIdx], cfg.chunk_size(),
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
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (!channelEof[ch] || !chunkData[ch].empty()) {
        allEof = false;
        break;
      }
    }

    if (allEof) {
      break;
    }

    // Find minimum number of events across all channels
    int nEventsInChunk = static_cast<int>(chunkData[0].size());
    bool chunkMismatch = false;
    for (int ch = 1; ch < cfg.n_channels(); ++ch) {
      int chEvents = static_cast<int>(chunkData[ch].size());
      if (chEvents != nEventsInChunk) {
        chunkMismatch = true;
        if (chEvents < nEventsInChunk) {
          nEventsInChunk = chEvents;
        }
      }
    }

    if (chunkMismatch) {
      std::cout << "WARNING: Event count mismatch in chunk " << chunkNumber << std::endl;
      std::cout << "         Per-channel event counts:" << std::endl;
      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        std::cout << "           ch" << ch << ": " << chunkData[ch].size()
                  << " events, EOF=" << static_cast<int>(channelEof[ch]) << std::endl;
      }
      std::cout << "         Processing minimum: " << nEventsInChunk << " events" << std::endl;
    }

    if (nEventsInChunk == 0) {
      break;
    }
    
    // Fill TTree with this chunk
    for (int evt = 0; evt < nEventsInChunk; ++evt) {

      if ( kSetEventLimit && evt >= cfg.max_events()) {
	std::cout<<"Reach the events limits (" <<evt<<")"<<std::endl;
	break;
      }
      
      std::vector<int> samplesThisEvent(cfg.n_channels(), 0);
      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        samplesThisEvent[ch] =
            static_cast<int>(chunkData[ch][evt].samples.size());
      }

      std::vector<std::string> consistencyIssues;
      bool skipEvent = false;
      const auto &refEvt = chunkData[0][evt];
      for (int ch = 1; ch < cfg.n_channels(); ++ch) {
        const auto &evtData = chunkData[ch][evt];
        if (evtData.eventCounter != refEvt.eventCounter) {
          std::ostringstream oss;
          oss << "eventCounter mismatch chunk " << chunkNumber << " event " << evt
              << " ch" << ch << " (" << evtData.eventCounter
              << " vs " << refEvt.eventCounter << ")";
          consistencyIssues.push_back(oss.str());
        }
        if (evtData.boardId != refEvt.boardId) {
          std::ostringstream oss;
          oss << "boardId mismatch chunk " << chunkNumber << " event " << evt
              << " ch" << ch << " (" << evtData.boardId
              << " vs " << refEvt.boardId << ")";
          consistencyIssues.push_back(oss.str());
        }
        if (evtData.channelId != static_cast<uint32_t>(ch)) {
          if (IsSpecialOverrideChannel(cfg, ch)) {
            if (!loggedSpecialChannelIdInfo) {
              std::cout << "INFO: special_channel_index " << ch
                        << " allows channelId mismatch (header "
                        << evtData.channelId << ")" << std::endl;
              loggedSpecialChannelIdInfo = true;
            }
          } else {
            std::ostringstream oss;
            oss << "channelId mismatch chunk " << chunkNumber << " event " << evt
                << " ch" << ch << " (" << evtData.channelId
                << " vs expected " << ch << ")";
            consistencyIssues.push_back(oss.str());
          }
        }
      }

      if (!consistencyIssues.empty()) {
        for (const auto &msg : consistencyIssues) {
          if (eventPolicy == EventPolicy::kWarn &&
              consistencyWarningCount >= kConsistencyWarnLimit) {
            continue;
          }
          std::cerr << ((eventPolicy == EventPolicy::kWarn) ? "WARNING: " : "ERROR: ")
                    << msg << std::endl;
          if (eventPolicy == EventPolicy::kWarn) {
            ++consistencyWarningCount;
            if (consistencyWarningCount == kConsistencyWarnLimit) {
              std::cerr << "WARNING: further consistency warnings suppressed"
                        << std::endl;
            }
          }
        }
        if (eventPolicy == EventPolicy::kError) {
          for (auto &fin : fins) {
            if (fin.is_open()) {
              fin.close();
            }
          }
          file->Close();
          delete file;
          return false;
        } else if (eventPolicy == EventPolicy::kSkip) {
          skipEvent = true;
          if (!loggedEventPolicyInfo) {
            std::cout << "INFO: skipping inconsistent events per event_policy=skip" << std::endl;
            loggedEventPolicyInfo = true;
          }
        }
      }

      const int maxSamples =
          *std::max_element(samplesThisEvent.begin(), samplesThisEvent.end());
      const bool hasMismatch =
          std::any_of(samplesThisEvent.begin(), samplesThisEvent.end(),
                      [maxSamples](int v) { return v != maxSamples; });

      if (hasMismatch && policy == NsamplesPolicy::kStrict) {
        std::cerr << "ERROR: nsamples mismatch in chunk " << chunkNumber
                  << " event " << evt << ", expected uniform sample counts"
                  << std::endl;
        for (auto &fin : fins) {
          if (fin.is_open()) {
            fin.close();
          }
        }
        file->Close();
        delete file;
        return false;
      }

      if (hasMismatch && !loggedNsamplesPadding &&
          policy == NsamplesPolicy::kPad) {
        std::cout << "INFO: nsamples mismatch detected in chunk " << chunkNumber
                  << ", padding shorter channels up to " << maxSamples
                  << " samples" << std::endl;
        loggedNsamplesPadding = true;
      }

      nsamplesBranch = maxSamples;
      nsamplesPerChannel = samplesThisEvent;

      if (timeAxis.size() != static_cast<size_t>(maxSamples)) {
        timeAxis.resize(maxSamples);
        for (int i = 0; i < maxSamples; ++i) {
          timeAxis[i] = static_cast<float>(i * cfg.tsample_ns);
        }
      }

      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        auto &evtData = chunkData[ch][evt];
        const int nsampCh = samplesThisEvent[ch];

        boardIds[ch] = evtData.boardId;
        channelIds[ch] = evtData.channelId;
        eventCounters[ch] = evtData.eventCounter;
	raw[ch] = std::move(evtData.samples);
	evtData.samples.clear();

        // Calculate pedestal
        const int nPed = std::min<int>(raw[ch].size(), pedWindow);
        double pedVal = 0.0;
        for (int i = 0; i < nPed; ++i) {
	  pedVal += raw[ch][i];
        }
        pedVal /= static_cast<double>(std::max(1, nPed));
        pedestals[ch] = static_cast<float>(pedVal);

        // Pedestal-subtracted waveform
        raw[ch].resize(maxSamples, pedestals[ch]);
        ped[ch].resize(maxSamples);
        for (int i = 0; i < nsampCh; ++i) {
          ped[ch][i] = raw[ch][i] - pedestals[ch] + pedTarget;
        }
        for (int i = nsampCh; i < maxSamples; ++i) {
          ped[ch][i] = pedTarget;
        }
      }

      if (!skipEvent) {
        eventIdx = totalEventsProcessed;
        tree->Fill();
        ++totalEventsProcessed;
      }
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
  std::string outputPath = BuildOutputPath(cfg.output_dir(), "root", cfg.root_file());
  if (!EnsureParentDirectory(outputPath)) {
    return false;
  }

  TFile *file = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outputPath << std::endl;
    return false;
  }

  CheckEventLimit(cfg);
  if ( kSetEventLimit ) {
    std::cout<<"The number of analyze event size is set = "<<cfg.max_events()<<std::endl;
  } else {
    std::cout<<"All events will be analyzed"<<std::endl;
  }
  
  std::cout << "Creating ROOT file: " << outputPath << std::endl;
  
  TTree *tree = new TTree(cfg.tree_name().c_str(), "Raw waveforms");

  int eventIdx = 0;
  int nChannelsBranch = cfg.n_channels();
  int nsamplesBranch = 0;
  float samplingNs = static_cast<float>(cfg.tsample_ns);
  float pedTarget = static_cast<float>(cfg.ped_target);
  int pedestalWindow = cfg.pedestal_window;
  const NsamplesPolicy policy = ResolveNsamplesPolicy(cfg.common.nsamples_policy);
  const EventPolicy eventPolicy = ResolveEventPolicy(cfg.event_policy);

  std::vector<float> timeAxis;
  std::vector<float> pedestals(cfg.n_channels(), 0.0f);
  std::vector<uint32_t> boardIds(cfg.n_channels(), 0);
  std::vector<uint32_t> channelIds(cfg.n_channels(), 0);
  std::vector<uint32_t> eventCounters(cfg.n_channels(), 0);
  std::vector<int> nsamplesPerChannel(cfg.n_channels(), 0);
  std::vector<std::vector<float>> raw(cfg.n_channels());
  std::vector<std::vector<float>> ped(cfg.n_channels());
  bool loggedSpecialChannelIdInfo = false;

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
  tree->Branch("nsamples_per_channel", &nsamplesPerChannel);

  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
    char bnameRaw[32];
    char bnamePed[32];
    std::snprintf(bnameRaw, sizeof(bnameRaw), "ch%02d_raw", ch);
    std::snprintf(bnamePed, sizeof(bnamePed), "ch%02d_ped", ch);
    tree->Branch(bnameRaw, &raw[ch]);
    tree->Branch(bnamePed, &ped[ch]);
  }

  std::vector<std::vector<AsciiEventBlock>> channelEvents(cfg.n_channels());
  size_t minEvents = 0;
  size_t maxEvents = 0;
  bool eventCountMismatch = false;

  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
    const std::string fname = BuildFileName(cfg, ch);
    if (!LoadAsciiChannelFile(fname, channelEvents[ch])) {
      file->Close();
      delete file;
      return false;
    }
    std::cout << "Loaded ASCII input " << fname << " with "
              << channelEvents[ch].size() << " event(s)." << std::endl;

    size_t evtCount = channelEvents[ch].size();
    if (ch == 0) {
      minEvents = evtCount;
      maxEvents = evtCount;
    } else {
      if (evtCount < minEvents) minEvents = evtCount;
      if (evtCount > maxEvents) maxEvents = evtCount;
      if (evtCount != channelEvents[0].size()) {
        eventCountMismatch = true;
      }
    }
  }

  if (eventCountMismatch) {
    std::cout << "WARNING: Event count mismatch detected across channels." << std::endl;
    std::cout << "         Will process only the minimum number of events: "
              << minEvents << std::endl;
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (channelEvents[ch].size() != minEvents) {
        std::cout << "         Channel " << ch << " has "
                  << channelEvents[ch].size() << " events (will use first "
                  << minEvents << ")" << std::endl;
      }
    }
  }

  if (minEvents == 0) {
    std::cerr << "ERROR: no events found in ASCII inputs." << std::endl;
    file->Close();
    delete file;
    return false;
  }

  size_t expectedEvents = minEvents;

  const int pedWindow = std::max(1, cfg.pedestal_window);
  bool loggedNsamplesPadding = false;
  bool loggedEventPolicyInfo = false;
  for (size_t evt = 0; evt < expectedEvents; ++evt) {

    if ( kSetEventLimit && evt >= cfg.max_events()) {
      std::cout<<"Reach the events limits (" <<evt<<")"<<std::endl;
      break;
    }

    std::vector<int> samplesThisEvent(cfg.n_channels(), 0);
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      samplesThisEvent[ch] = static_cast<int>(channelEvents[ch][evt].samples.size());
    }

    const int maxSamples =
        *std::max_element(samplesThisEvent.begin(), samplesThisEvent.end());
    const bool hasMismatch =
        std::any_of(samplesThisEvent.begin(), samplesThisEvent.end(),
                    [maxSamples](int v) { return v != maxSamples; });

    if (hasMismatch && policy == NsamplesPolicy::kStrict) {
      std::cerr << "ERROR: nsamples mismatch at ASCII event " << evt
                << ", expected uniform sample counts across channels" << std::endl;
      file->Close();
      delete file;
      return false;
    }

    if (hasMismatch && !loggedNsamplesPadding &&
        policy == NsamplesPolicy::kPad) {
      std::cout << "INFO: nsamples mismatch detected at ASCII event " << evt
                << ", padding shorter channels up to " << maxSamples
                << " samples" << std::endl;
      loggedNsamplesPadding = true;
    }

    nsamplesBranch = maxSamples;
    nsamplesPerChannel = samplesThisEvent;
    if (timeAxis.size() != static_cast<size_t>(maxSamples)) {
      timeAxis.resize(maxSamples);
      for (size_t i = 0; i < static_cast<size_t>(maxSamples); ++i) {
        timeAxis[i] = static_cast<float>(i * cfg.tsample_ns);
      }
    }

    std::vector<std::string> consistencyIssues;
    const auto &refBlock = channelEvents[0][evt];
    for (int ch = 1; ch < cfg.n_channels(); ++ch) {
      const auto &block = channelEvents[ch][evt];
      if (block.eventCounter != refBlock.eventCounter) {
        std::ostringstream oss;
        oss << "eventCounter mismatch ASCII event " << evt << " ch" << ch
            << " (" << block.eventCounter << " vs " << refBlock.eventCounter << ")";
        consistencyIssues.push_back(oss.str());
      }
      if (block.boardId != refBlock.boardId) {
        std::ostringstream oss;
        oss << "boardId mismatch ASCII event " << evt << " ch" << ch
            << " (" << block.boardId << " vs " << refBlock.boardId << ")";
        consistencyIssues.push_back(oss.str());
      }
      if (block.channelId != static_cast<uint32_t>(ch)) {
        if (IsSpecialOverrideChannel(cfg, ch)) {
          if (!loggedSpecialChannelIdInfo) {
            std::cout << "INFO: special_channel_index " << ch
                      << " allows channelId mismatch (header "
                      << block.channelId << ")" << std::endl;
            loggedSpecialChannelIdInfo = true;
          }
        } else {
          std::ostringstream oss;
          oss << "channelId mismatch ASCII event " << evt << " ch" << ch
              << " (" << block.channelId << " vs expected " << ch << ")";
          consistencyIssues.push_back(oss.str());
        }
      }
    }

    bool skipEvent = false;
    if (!consistencyIssues.empty()) {
      for (const auto &msg : consistencyIssues) {
        std::cerr << ((eventPolicy == EventPolicy::kWarn) ? "WARNING: " : "ERROR: ")
                  << msg << std::endl;
      }
      if (eventPolicy == EventPolicy::kError) {
        file->Close();
        delete file;
        return false;
      } else if (eventPolicy == EventPolicy::kSkip) {
        skipEvent = true;
        if (!loggedEventPolicyInfo) {
          std::cout << "INFO: skipping inconsistent events per event_policy=skip" << std::endl;
          loggedEventPolicyInfo = true;
        }
      }
    }

    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      const auto &block = channelEvents[ch][evt];
      const int nsampCh = samplesThisEvent[ch];
      boardIds[ch] = block.boardId;
      channelIds[ch] = block.channelId;
      eventCounters[ch] = block.eventCounter;
      raw[ch] = std::move(block.samples);

      const int nPed = std::min<int>(raw[ch].size(), pedWindow);
      double pedVal = 0.0;
      for (int i = 0; i < nPed; ++i) {
        pedVal += raw[ch][i];
      }
      pedVal /= static_cast<double>(std::max(1, nPed));
      pedestals[ch] = static_cast<float>(pedVal);

      raw[ch].resize(maxSamples, pedestals[ch]);
      ped[ch].resize(maxSamples);
      for (int i = 0; i < nsampCh; ++i) {
        ped[ch][i] = raw[ch][i] - pedestals[ch] + pedTarget;
      }
      for (int i = nsampCh; i < maxSamples; ++i) {
        ped[ch][i] = pedTarget;
      }
    }

    if (!skipEvent) {
      eventIdx = static_cast<int>(evt);
      tree->Fill();
    }
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
            << "  --nsamples-policy POLICY  nsamples handling: 'strict' or 'pad'\n"
            << "  --event-policy POLICY     event consistency: 'error', 'warn', or 'skip'\n"
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
        cfg.set_n_channels(std::stoi(val));
      } catch (const std::exception &) {
        std::cerr << "ERROR: invalid integer for --channels" << std::endl;
        return CliOutcome::kError;
      }
    } else if (arg == "--root") {
      const char *val = requireValue("--root");
      if (!val) {
        return CliOutcome::kError;
      }
      cfg.set_root_file(val);
    } else if (arg == "--nsamples-policy") {
      const char *val = requireValue("--nsamples-policy");
      if (!val) {
        return CliOutcome::kError;
      }
      cfg.common.nsamples_policy = val;
    } else if (arg == "--event-policy") {
      const char *val = requireValue("--event-policy");
      if (!val) {
        return CliOutcome::kError;
      }
      cfg.event_policy = val;
    } else if (arg == "--ascii") {
      cfg.input_is_ascii = true;
    } else if (arg == "--binary") {
      cfg.input_is_ascii = false;
    } else if (arg == "--parallel") {
      // Force parallel mode by ensuring max_cores > 1
      if (cfg.max_cores() < 2) {
        cfg.set_max_cores(2);
      }
    } else if (arg == "--chunk-size") {
      const char *val = requireValue("--chunk-size");
      if (!val) {
        return CliOutcome::kError;
      }
      try {
        cfg.set_chunk_size(std::stoi(val));
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
        cfg.set_max_cores(std::stoi(val));
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
  const std::string defaultPath = "converter_config.json";
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
    } else if (cfg.max_cores() > 1) {
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
