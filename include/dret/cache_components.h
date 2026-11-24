//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/4/26.
//

#pragma once

#include <filesystem>
#include <fstream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

namespace dret {

// A cached object is stored as its components: contiguous ranges of the object's
// own serialization, each in a file of its own, so every object that contains a
// component reads the same file. The object's format does not change -- loading
// concatenates the files in order and hands the stream to the object's load().
//
// A component's key names every parameter its construction reads (built with
// PrefixedKey from the common keys) and its type hash names its encoding. Two
// objects share a component exactly when their keys and types agree; equal bytes
// never decide it.
struct Component {
  std::string key;         // cache key, parameters included
  std::string type;        // type hash of the component's own type
  std::size_t bytes = 0;   // length of its range in the object's serialization
};

template <typename T>
std::string TypeHash() {
  return sdsl::util::class_to_hash(T{});
}

// Bytes t serializes to. Called on a base subobject it gives the length of that
// base's range, since every class here serializes its bases first, in order.
template <typename T>
std::size_t SerializedSize(const T& t) {
  return sdsl::size_in_bytes(t);
}

// The file a component lives in. For a single-component object this is the file
// sdsl::store_to_cache(object, key, config, true) writes.
inline std::string ComponentFile(const Component& t_component, const sdsl::cache_config& t_config) {
  return sdsl::cache_file_name(t_component.key + "_" + t_component.type, t_config);
}

inline bool ComponentsExist(const std::vector<Component>& t_components, const sdsl::cache_config& t_config) {
  for (const auto& c : t_components)
    if (!std::filesystem::exists(ComponentFile(c, t_config))) return false;
  return true;
}

namespace internal {

inline std::string ReadFile(const std::string& t_path) {
  std::ifstream in(t_path, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// Write through a temporary and rename, so a concurrent reader (or a second
// builder of the same shared component) never sees a partial file.
inline void WriteFileAtomic(const std::string& t_path, std::string_view t_bytes) {
  const auto tmp = t_path + ".tmp." + std::to_string(::getpid());
  {
    std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
    out.write(t_bytes.data(), static_cast<std::streamsize>(t_bytes.size()));
    if (!out) throw std::runtime_error("cannot write cache file " + tmp);
  }
  std::filesystem::rename(tmp, t_path);
}

}  // namespace internal

// Store t as t_components (in serialization order, with their byte lengths set).
// A component file that already exists is kept: its key carries every parameter
// its construction reads, so it holds the same structure. Its length must match,
// or the key is missing a parameter -- a bug in the layout, which throws rather
// than overwrite. The bytes themselves are not compared: SDSL vectors leave the
// unused bits of their last word uninitialized, so two builds of one structure
// can differ there (on page, GCDA's and GCDA-differential's sampled trees did).
template <typename T>
void StoreComponents(const T& t, const std::vector<Component>& t_components, const sdsl::cache_config& t_config) {
  std::ostringstream buffer;
  sdsl::serialize(t, buffer);
  const auto bytes = std::move(buffer).str();

  const auto total = std::accumulate(t_components.begin(), t_components.end(), std::size_t{0},
                                     [](std::size_t s, const Component& c) { return s + c.bytes; });
  if (total != bytes.size())
    throw std::logic_error("component layout covers " + std::to_string(total) + " of " +
                           std::to_string(bytes.size()) + " serialized bytes (" + t_components.front().key + ")");

  std::size_t offset = 0;
  for (const auto& c : t_components) {
    const std::string_view part(bytes.data() + offset, c.bytes);
    offset += c.bytes;
    const auto path = ComponentFile(c, t_config);
    if (std::filesystem::exists(path)) {
      if (std::filesystem::file_size(path) != part.size())
        throw std::logic_error("cache component " + path + " exists with another length (" +
                               std::to_string(std::filesystem::file_size(path)) + " vs " +
                               std::to_string(part.size()) + " bytes): its key misses a parameter of its construction");
      continue;
    }
    internal::WriteFileAtomic(path, part);
  }
}

// Load t from t_components. False if any is missing. Throws if the components
// do not add up to exactly one object (a layout that does not match the format).
template <typename T>
bool LoadComponents(T& t, const std::vector<Component>& t_components, const sdsl::cache_config& t_config) {
  if (!ComponentsExist(t_components, t_config)) return false;
  std::string bytes;
  for (const auto& c : t_components) bytes += internal::ReadFile(ComponentFile(c, t_config));
  std::istringstream in(bytes);
  t.load(in);
  if (!in || static_cast<std::size_t>(in.tellg()) != bytes.size())
    throw std::logic_error("cache components of " + t_components.front().key +
                           " do not form one object (read " + std::to_string(in ? std::size_t(in.tellg()) : 0) +
                           " of " + std::to_string(bytes.size()) + " bytes)");
  return true;
}

}  // namespace dret
