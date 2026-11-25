//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The one place an index's structures are written to and read from a stream. Every component
// serializes and loads through these functions (FieldsComponent for the field lists in dret/fields.h,
// IndexBaseWithExternalStorage for the loads, and the hand-written SAPhi and RLCSA), so the stream
// format is chosen here: today it is SDSL's own, and another one (e.g. cereal archives, which SDSL
// supports when cereal is installed) would start by changing these bodies.
//
// What does not go through here:
//   - sizes: sdsl::size_in_bytes and the structure tree, which measure a structure, not an archive;
//   - the cache files construction writes and the Config loads read (sdsl::store_to_cache /
//     load_from_cache), and the cache components (dret/cache_components.h), whose layouts are byte
//     ranges of SDSL's serialization;
//   - the inner structures' own serialize/load (grammars, PDL tree cores, sr-index, rlcsa), which
//     these functions call.
//

#pragma once

#include <cstddef>
#include <istream>
#include <ostream>
#include <string>

#include <sdsl/io.hpp>
#include <sdsl/structure_tree.hpp>

namespace dret::archive {

// A structure (anything with serialize/load).
template <typename T>
std::size_t Write(const T& t_item, std::ostream& t_out, sdsl::structure_tree_node* t_v, const std::string& t_name) {
  return sdsl::serialize(t_item, t_out, t_v, t_name);
}

// A structure a component points to; a null one is written as a default-constructed T, which Read
// reads back.
template <typename T>
std::size_t WriteOrEmpty(const T* t_item,
                         std::ostream& t_out,
                         sdsl::structure_tree_node* t_v,
                         const std::string& t_name) {
  return t_item ? Write(*t_item, t_out, t_v, t_name) : sdsl::serialize_empty_object<T>(t_out, t_v, t_name);
}

template <typename T>
void Read(T& t_item, std::istream& t_in) {
  sdsl::load(t_item, t_in);
}

// A plain value (an integer, a string).
template <typename T>
std::size_t WriteValue(const T& t_value,
                       std::ostream& t_out,
                       sdsl::structure_tree_node* t_v,
                       const std::string& t_name) {
  return sdsl::write_member(t_value, t_out, t_v, t_name);
}

template <typename T>
void ReadValue(T& t_value, std::istream& t_in) {
  sdsl::read_member(t_value, t_in);
}

}  // namespace dret::archive
