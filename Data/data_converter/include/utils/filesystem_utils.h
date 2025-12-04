#pragma once

#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <string>

// Recursively create directory path if it does not exist.
inline bool CreateDirectoryIfNeeded(const std::string &path) {
  if (path.empty() || path == ".") {
    return true;
  }

  struct stat info;
  if (stat(path.c_str(), &info) == 0) {
    if (info.st_mode & S_IFDIR) {
      return true;
    }
    std::cerr << "ERROR: Path exists but is not a directory: " << path << std::endl;
    return false;
  }

  size_t pos = path.find_last_of('/');
  if (pos != std::string::npos && pos > 0) {
    std::string parent = path.substr(0, pos);
    if (!CreateDirectoryIfNeeded(parent)) {
      return false;
    }
  }

  if (mkdir(path.c_str(), 0755) != 0) {
    std::cerr << "ERROR: Failed to create directory: " << path << std::endl;
    return false;
  }

  return true;
}

// Build path: output_dir/subdir/filename (handles absolute filename and ".").
inline std::string BuildOutputPath(const std::string &output_dir,
                                   const std::string &subdir,
                                   const std::string &filename) {
  if (!filename.empty() && filename[0] == '/') {
    return filename;
  }
  if (output_dir.empty() || output_dir == ".") {
    return filename;
  }

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

// Backward-compatible alias for BuildOutputPath used in export stage.
inline std::string BuildPath(const std::string &output_dir,
                             const std::string &subdir,
                             const std::string &filename) {
  return BuildOutputPath(output_dir, subdir, filename);
}
