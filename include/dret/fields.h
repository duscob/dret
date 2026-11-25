//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Components that list their fields once. A component's load, serialization and size report each walk
// its structures, and kept by hand they drift apart: a structure saved under one key and loaded under
// another, one left out of the size, one read back in a different order than it was written. A
// FieldsComponent states the list once, as a static Fields(self), and derives the three from it, in
// the list's order:
//
//   Item(name, ptr, key)            a structure in the shared storage under key(keys), with the type
//                                   hash (loadItemPtr); sized under name
//   Rank<TBv>(name, ptr, key)       the rank / select support of the bitvector at key(keys)
//   Select<TBv>(name, ptr, key)       (loadBVRank / loadBVSelect)
//   Components(name, ptr, comps)    a structure cached as components, comps(keys) (loadComponentsPtr)
//   NDoc(n)                         the number of documents (loadNDoc); serialized, not sized
//   Part(name, part)                a sub-component, loaded, serialized and sized by its own methods
//   SriPart(name, part)             the same for an sr-index structure, which loads from a stream with
//                                   its own default keys, not the component's
//
// Each takes how it is sized (SizeAs, a template argument): kItem (one entry, name, sdsl::size_in_bytes), kCollect (the
// structure's breakdown, collectSizes, under a prefix), kReport (its own GetSizeReport, under a
// prefix) or kNone. A null structure is serialized as an empty object and not sized.
//
// The reads and writes themselves go through dret/archive.h, which fixes the stream format.
//
// Fields(self) is a static template over the component, const for serialize and the size report,
// mutable for the load, so the same list reads and writes the members. A component that derives
// something from what it loaded does so in OnLoad().
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "dret/archive.h"
#include "dret/config.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret {

enum class SizeAs { kItem, kCollect, kReport, kNone };

namespace fields {

enum class Kind { kItem, kRank, kSelect, kComponents };

// A structure the component points to; TPtr is the member (const T*), const or not.
template <Kind kKind, SizeAs kSizeAs, typename TBv, typename TPtr, typename TKey>
struct Ref {
  const char* name;
  TPtr& ptr;
  TKey key;
  const char* prefix;
};

template <typename TN>
struct NDocField {
  TN& n;
};

template <SizeAs kSizeAs, bool kWithKeys, typename TPart>
struct PartField {
  const char* name;
  TPart& part;
  const char* prefix;
};

}  // namespace fields

template <SizeAs kSizeAs = SizeAs::kItem, typename TPtr, typename TKey>
auto Item(const char* t_name, TPtr& t_ptr, TKey t_key, const char* t_prefix = "") {
  return fields::Ref<fields::Kind::kItem, kSizeAs, void, TPtr, TKey>{t_name, t_ptr, std::move(t_key), t_prefix};
}

template <typename TBv, typename TPtr, typename TKey>
auto Rank(const char* t_name, TPtr& t_ptr, TKey t_key) {
  return fields::Ref<fields::Kind::kRank, SizeAs::kItem, TBv, TPtr, TKey>{t_name, t_ptr, std::move(t_key), ""};
}

template <typename TBv, typename TPtr, typename TKey>
auto Select(const char* t_name, TPtr& t_ptr, TKey t_key) {
  return fields::Ref<fields::Kind::kSelect, SizeAs::kItem, TBv, TPtr, TKey>{t_name, t_ptr, std::move(t_key), ""};
}

template <SizeAs kSizeAs = SizeAs::kItem, typename TPtr, typename TComponents>
auto Components(const char* t_name, TPtr& t_ptr, TComponents t_components, const char* t_prefix = "") {
  return fields::Ref<fields::Kind::kComponents, kSizeAs, void, TPtr, TComponents>{
      t_name, t_ptr, std::move(t_components), t_prefix};
}

template <typename TN>
auto NDoc(TN& t_n) {
  return fields::NDocField<TN>{t_n};
}

template <SizeAs kSizeAs = SizeAs::kReport, typename TPart>
auto Part(const char* t_name, TPart& t_part, const char* t_prefix = "") {
  return fields::PartField<kSizeAs, true, TPart>{t_name, t_part, t_prefix};
}

// An sr-index structure (the r-index, the sr-index): from a stream it loads with its own default keys.
template <SizeAs kSizeAs = SizeAs::kItem, typename TPart>
auto SriPart(const char* t_name, TPart& t_part) {
  return fields::PartField<kSizeAs, false, TPart>{t_name, t_part, ""};
}

// The base of a component with a field list. TDerived declares
//   template <typename TSelf> static auto Fields(TSelf& t_self);   (a std::tuple of the fields above)
// and, if it derives anything from what it loaded, void OnLoad().
template <typename TDerived, typename TStorage, uint8_t t_width>
class FieldsComponent : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  using Base::Base;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(self()));
    size_type written = 0;
    std::apply([&](const auto&... f) { ((written += serializeField(f, out, child)), ...); }, TDerived::Fields(self()));
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    std::apply([&](const auto&... f) { (sizeField(f, r), ...); }, TDerived::Fields(self()));
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    auto& derived = static_cast<TDerived&>(*this);
    std::apply([&](auto&&... f) { (loadField(f, t_source, t_keys), ...); }, TDerived::Fields(derived));
    if constexpr (requires { derived.OnLoad(); }) derived.OnLoad();
  }

 private:
  const TDerived& self() const { return static_cast<const TDerived&>(*this); }

  template <typename TPtr>
  using Pointee = std::remove_cvref_t<decltype(*std::declval<TPtr>())>;

  //~~~~~~~ load

  template <fields::Kind kKind, SizeAs kSizeAs, typename TBv, typename TPtr, typename TKey>
  void loadField(fields::Ref<kKind, kSizeAs, TBv, TPtr, TKey>& t_f, TSource& t_source, const JSON& t_keys) {
    using T = Pointee<TPtr>;
    if constexpr (kKind == fields::Kind::kItem) {
      t_f.ptr = this->template loadItemPtr<T>(t_f.key(t_keys), t_source, true);
    } else if constexpr (kKind == fields::Kind::kRank) {
      t_f.ptr = &this->template loadBVRank<TBv>(t_f.key(t_keys), t_source, true).get();
    } else if constexpr (kKind == fields::Kind::kSelect) {
      t_f.ptr = &this->template loadBVSelect<TBv>(t_f.key(t_keys), t_source, true).get();
    } else {
      t_f.ptr = this->template loadComponentsPtr<T>(t_f.key(t_keys), t_source);
    }
  }

  template <typename TN>
  void loadField(fields::NDocField<TN>& t_f, TSource& t_source, const JSON& t_keys) {
    t_f.n = this->loadNDoc(t_keys, t_source);
  }

  // From a stream, a part takes the component's keys, an sr-index one its own.
  template <SizeAs kSizeAs, bool kWithKeys, typename TPart>
  void loadField(fields::PartField<kSizeAs, kWithKeys, TPart>& t_f, TSource& t_source, const JSON& t_keys) {
    std::visit(
        [&](auto&& tt_source) {
          using TS = std::decay_t<decltype(tt_source.get())>;
          if constexpr (std::is_base_of_v<std::istream, TS> && kWithKeys)
            t_f.part.load(tt_source.get(), t_keys);
          else
            t_f.part.load(tt_source.get());
        },
        t_source);
  }

  //~~~~~~~ serialize

  template <fields::Kind kKind, SizeAs kSizeAs, typename TBv, typename TPtr, typename TKey>
  static size_type serializeField(const fields::Ref<kKind, kSizeAs, TBv, TPtr, TKey>& t_f,
                                  std::ostream& out,
                                  sdsl::structure_tree_node* v) {
    return archive::WriteOrEmpty<Pointee<TPtr>>(t_f.ptr, out, v, t_f.name);
  }

  template <typename TN>
  static size_type serializeField(const fields::NDocField<TN>& t_f, std::ostream& out, sdsl::structure_tree_node* v) {
    sdsl::int_vector<64> n_doc(1, t_f.n);
    return archive::Write(n_doc, out, v, "n_doc");
  }

  template <SizeAs kSizeAs, bool kWithKeys, typename TPart>
  static size_type serializeField(const fields::PartField<kSizeAs, kWithKeys, TPart>& t_f,
                                  std::ostream& out,
                                  sdsl::structure_tree_node* v) {
    return t_f.part.serialize(out, v, t_f.name);
  }

  //~~~~~~~ size report

  template <fields::Kind kKind, SizeAs kSizeAs, typename TBv, typename TPtr, typename TKey>
  static void sizeField(const fields::Ref<kKind, kSizeAs, TBv, TPtr, TKey>& t_f, SizeReport& t_report) {
    if (!t_f.ptr) return;
    if constexpr (kSizeAs == SizeAs::kItem) {
      append(t_report, t_f.name, sdsl::size_in_bytes(*t_f.ptr));
    } else if constexpr (kSizeAs == SizeAs::kCollect) {
      collectSizes(t_report, *t_f.ptr, t_f.prefix);
    } else if constexpr (kSizeAs == SizeAs::kReport) {
      for (const auto& [name, bytes] : t_f.ptr->GetSizeReport()) append(t_report, t_f.prefix + name, bytes);
    }
  }

  template <typename TN>
  static void sizeField(const fields::NDocField<TN>&, SizeReport&) {}

  template <SizeAs kSizeAs, bool kWithKeys, typename TPart>
  static void sizeField(const fields::PartField<kSizeAs, kWithKeys, TPart>& t_f, SizeReport& t_report) {
    if constexpr (kSizeAs == SizeAs::kItem) {
      append(t_report, t_f.name, sdsl::size_in_bytes(t_f.part));
    } else if constexpr (kSizeAs == SizeAs::kReport) {
      for (const auto& [name, bytes] : t_f.part.GetSizeReport()) append(t_report, t_f.prefix + name, bytes);
    }
  }
};

}  // namespace dret
