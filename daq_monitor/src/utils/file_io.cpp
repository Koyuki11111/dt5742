#include "utils/file_io.h"

#include <iostream>

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
