#pragma once

#include <string>
#include <vector>

#include "simdjson.h"

inline bool ParseJsonFile(const std::string &path,
                          simdjson::dom::parser &parser,
                          simdjson::dom::element &root,
                          std::string *errorMessage = nullptr) {
  auto padded = simdjson::padded_string::load(path);
  if (padded.error()) {
    if (errorMessage) {
      *errorMessage = "cannot open config file: " + path;
    }
    return false;
  }

  auto doc = parser.parse(padded.value());
  if (doc.error()) {
    if (errorMessage) {
      *errorMessage = "JSON parse error in " + path + ": " +
                      std::string(simdjson::error_message(doc.error()));
    }
    return false;
  }

  root = doc.value();
  return true;
}

inline bool GetObject(const simdjson::dom::element &parent,
                      const std::string &key,
                      simdjson::dom::element &out) {
  auto obj = parent[key];
  if (obj.error()) {
    return false;
  }
  out = obj.value();
  return out.is_object();
}

inline bool GetString(const simdjson::dom::element &parent,
                      const std::string &key,
                      std::string &out) {
  auto val = parent[key];
  if (val.error()) {
    return false;
  }
  auto str = val.get_string();
  if (str.error()) {
    return false;
  }
  out = std::string(str.value());
  return true;
}

inline bool GetNumber(const simdjson::dom::element &parent,
                      const std::string &key,
                      double &out) {
  auto val = parent[key];
  if (val.error()) {
    return false;
  }
  if (val.is<int64_t>()) {
    out = static_cast<double>(int64_t(val));
    return true;
  }
  if (val.is<uint64_t>()) {
    out = static_cast<double>(uint64_t(val));
    return true;
  }
  if (val.is<double>()) {
    out = double(val);
    return true;
  }
  return false;
}

inline bool GetBool(const simdjson::dom::element &parent,
                    const std::string &key,
                    bool &out) {
  auto val = parent[key];
  if (val.error()) {
    return false;
  }
  auto b = val.get<bool>();
  if (b.error()) {
    return false;
  }
  out = b.value();
  return true;
}

inline bool GetFloatArray(const simdjson::dom::element &parent,
                          const std::string &key,
                          std::vector<float> &out) {
  auto val = parent[key];
  if (val.error()) {
    return false;
  }
  auto arr = val.get_array();
  if (arr.error()) {
    return false;
  }
  out.clear();
  for (auto v : arr.value()) {
    double num = 0.0;
    if (v.is<int64_t>()) {
      num = static_cast<double>(int64_t(v));
    } else if (v.is<uint64_t>()) {
      num = static_cast<double>(uint64_t(v));
    } else if (v.is<double>()) {
      num = double(v);
    } else {
      continue;
    }
    out.push_back(static_cast<float>(num));
  }
  return true;
}

inline bool GetIntArray(const simdjson::dom::element &parent,
                        const std::string &key,
                        std::vector<int> &out) {
  auto val = parent[key];
  if (val.error()) {
    return false;
  }
  auto arr = val.get_array();
  if (arr.error()) {
    return false;
  }
  out.clear();
  for (auto v : arr.value()) {
    if (v.is<int64_t>()) {
      out.push_back(static_cast<int>(int64_t(v)));
    } else if (v.is<uint64_t>()) {
      out.push_back(static_cast<int>(uint64_t(v)));
    }
  }
  return true;
}
