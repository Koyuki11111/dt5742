#pragma once

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <string>
#include <vector>

static const size_t HEADER_WORDS = 8;
static const size_t HEADER_BYTES = HEADER_WORDS * sizeof(uint32_t);

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

std::string TrimCopy(const std::string &text);
bool TryParseInt(const std::string &text, int &value);
bool TryParseUint(const std::string &text, uint32_t &value);
bool ReadHeader(std::ifstream &fin, ChannelHeader &out);
bool LoadAsciiChannelFile(const std::string &path, std::vector<AsciiEventBlock> &events);
bool ReadChannelChunk(std::ifstream &fin, std::mutex &fileMutex, int chunkSize,
                      std::vector<BinaryEventData> &events, uint8_t &eofReached);
