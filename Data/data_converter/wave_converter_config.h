#pragma once

#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

struct WaveConverterConfig {
  // Common fields (from pipeline_config.json "common" section)
  std::string output_dir = "output";
  int n_channels = 16;
  int max_cores = 8;
  int chunk_size = 100;
  std::string temp_dir = "./temp";
  std::string root_file = "waveforms.root";      // waveforms_root in config
  std::string tree_name = "Waveforms";            // waveforms_tree in config

  // Stage1 specific fields
  std::string input_pattern = "wave_%d.dat";
  std::string input_dir = ".";
  bool input_is_ascii = false;
  std::string special_channel_file = "TR_0_0.dat";
  bool enable_special_override = true;
  int special_channel_index = 3;
  double tsample_ns = 0.2;
  int pedestal_window = 100;
  double ped_target = 3500.0;
};

inline bool ExtractStringValue(const std::string &text, const std::string &key,
                               std::string &value) {
  const std::string pattern =
      "\"" + key + "\"\\s*:\\s*\"((?:\\\\\"|[^\"])*)\"";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    value = match[1].str();
    // Unescape quotation marks
    std::string restored;
    restored.reserve(value.size());
    for (size_t i = 0; i < value.size(); ++i) {
      if (value[i] == '\\' && i + 1 < value.size()) {
        char next = value[i + 1];
        if (next == '\\' || next == '"') {
          restored.push_back(next);
          ++i;
          continue;
        }
      }
      restored.push_back(value[i]);
    }
    value = restored;
    return true;
  }
  return false;
}

inline bool ExtractNumberValue(const std::string &text, const std::string &key,
                               double &value) {
  const std::string pattern =
      "\"" + key + "\"\\s*:\\s*(-?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?)";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    value = std::stod(match[1].str());
    return true;
  }
  return false;
}

inline bool ExtractBoolValue(const std::string &text, const std::string &key,
                             bool &value) {
  const std::string pattern = "\"" + key + "\"\\s*:\\s*(true|false)";
  std::regex re(pattern);
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    value = (match[1] == "true");
    return true;
  }
  return false;
}

// Extract a section from JSON (e.g., "common", "stage1", "stage2")
inline bool ExtractSection(const std::string &text, const std::string &section,
                          std::string &sectionContent) {
  // Pattern: "section": { ... }
  // We need to handle nested braces properly
  std::string pattern = "\"" + section + "\"\\s*:\\s*\\{";
  std::regex re(pattern);
  std::smatch match;

  if (!std::regex_search(text, match, re)) {
    return false;
  }

  // Find the starting position of the section content
  size_t start = match.position() + match.length();

  // Count braces to find the end of this section
  int braceCount = 1;
  size_t end = start;
  bool inString = false;
  bool escaped = false;

  while (end < text.size() && braceCount > 0) {
    char c = text[end];

    if (escaped) {
      escaped = false;
      end++;
      continue;
    }

    if (c == '\\') {
      escaped = true;
      end++;
      continue;
    }

    if (c == '"') {
      inString = !inString;
    } else if (!inString) {
      if (c == '{') {
        braceCount++;
      } else if (c == '}') {
        braceCount--;
      }
    }

    end++;
  }

  if (braceCount != 0) {
    return false;  // Mismatched braces
  }

  // Extract the content (excluding the final '}')
  sectionContent = text.substr(start, end - start - 1);
  return true;
}

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

// Helper function to build output path with subdirectory
inline std::string BuildOutputPath(const std::string &output_dir,
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

inline bool LoadConfigFromJson(const std::string &path,
                               WaveConverterConfig &cfg,
                               std::string *errorMessage = nullptr) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    if (errorMessage) {
      *errorMessage = "cannot open config file: " + path;
    }
    return false;
  }

  std::string content((std::istreambuf_iterator<char>(fin)),
                      std::istreambuf_iterator<char>());

  if (content.empty()) {
    if (errorMessage) {
      *errorMessage = "config file is empty: " + path;
    }
    return false;
  }

  // Extract sections
  std::string commonSection, stage1Section;
  bool hasCommon = ExtractSection(content, "common", commonSection);
  bool hasStage1 = ExtractSection(content, "stage1", stage1Section);

  if (!hasCommon && !hasStage1) {
    if (errorMessage) {
      *errorMessage = "config file missing 'common' and 'stage1' sections";
    }
    return false;
  }

  // Parse common section
  if (hasCommon) {
    std::string strValue;
    double numValue = 0.0;

    if (ExtractStringValue(commonSection, "output_dir", strValue)) {
      cfg.output_dir = strValue;
    }
    if (ExtractNumberValue(commonSection, "n_channels", numValue)) {
      cfg.n_channels = static_cast<int>(numValue);
    }
    if (ExtractNumberValue(commonSection, "max_cores", numValue)) {
      cfg.max_cores = static_cast<int>(numValue);
    }
    if (ExtractNumberValue(commonSection, "chunk_size", numValue)) {
      cfg.chunk_size = static_cast<int>(numValue);
    }
    if (ExtractStringValue(commonSection, "temp_dir", strValue)) {
      cfg.temp_dir = strValue;
    }
    if (ExtractStringValue(commonSection, "waveforms_root", strValue)) {
      cfg.root_file = strValue;
    }
    if (ExtractStringValue(commonSection, "waveforms_tree", strValue)) {
      cfg.tree_name = strValue;
    }
  }

  // Parse stage1 section
  if (hasStage1) {
    std::string strValue;
    double numValue = 0.0;
    bool boolValue = false;

    if (ExtractStringValue(stage1Section, "input_pattern", strValue)) {
      cfg.input_pattern = strValue;
    }
    if (ExtractStringValue(stage1Section, "input_dir", strValue)) {
      cfg.input_dir = strValue;
    }
    if (ExtractBoolValue(stage1Section, "input_is_ascii", boolValue)) {
      cfg.input_is_ascii = boolValue;
    }
    if (ExtractStringValue(stage1Section, "special_channel_file", strValue)) {
      cfg.special_channel_file = strValue;
    }
    if (ExtractBoolValue(stage1Section, "enable_special_override", boolValue)) {
      cfg.enable_special_override = boolValue;
    }
    if (ExtractNumberValue(stage1Section, "special_channel_index", numValue)) {
      cfg.special_channel_index = static_cast<int>(numValue);
    }
    if (ExtractNumberValue(stage1Section, "tsample_ns", numValue)) {
      cfg.tsample_ns = numValue;
    }
    if (ExtractNumberValue(stage1Section, "pedestal_window", numValue)) {
      cfg.pedestal_window = static_cast<int>(numValue);
    }
    if (ExtractNumberValue(stage1Section, "ped_target", numValue)) {
      cfg.ped_target = numValue;
    }
  }

  return true;
}


