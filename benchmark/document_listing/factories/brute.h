//
// Brute-baseline factory header — type aliases for the dret::DocListIdxBrute
// family wrapping an r-index or sr-index.
//
// Brute baselines bypass the construct-then-load auto-build pattern: their
// underlying sri::RIndex / sri::SrIndexValidArea caches are managed
// externally (in the query factory: via the lazy rIndex() / srIndex()
// accessors; in bm_build_items: via DocListIdxBrute::construct() which
// builds the locate-index from the raw data file). The factory.h facade
// special-cases brute and does not call construct() on the index itself —
// only load(). bm_build_items uses DocListIdxBrute<> directly with
// BM_ConstructBruteIdx<> which exercises construct(idx, config).
//
// Brute force over a raw-range backend (dret::DocListIdxBruteBackend: the
// plain DA or the RLCSA) does follow construct-then-load, as the RMQ cores do.
//

#pragma once

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_brute_backend.h"
#include "dret/rmq/rmq_get_doc_rlcsa.h"

namespace bench::factories::brute {

// r-index-backed brute baseline.
template <typename TStorage>
using IdxR = dret::DocListIdxBrute<
    TStorage,
    dret::Alphabet<>,
    sri::RIndex<TStorage>>;

// sr-index-backed brute baseline. The sampling parameter is threaded
// through the index constructor as a runtime argument.
template <typename TStorage>
using IdxSr = dret::DocListIdxBrute<
    TStorage,
    dret::Alphabet<>,
    sri::SrIndexValidArea<TStorage>>;

// The same shape, named for construct mode. DocListIdxBrute's own default
// locate index is sri::SrIdxGeneric<SrIndexValidArea<...>, 16>, which binds
// the subsample rate at compile time; construct mode has to cover whatever
// sampling_size values the spec lists, so it names SrIndexValidArea directly
// and passes the rate to the constructor. That is all SrIdxGeneric does
// internally -- it exists only to expose a one-argument (storage) ctor.
template <typename TStorage>
using IdxSrConstruct = IdxSr<TStorage>;

// Brute force over the plain DA and over the RLCSA.
template <typename TStorage>
using IdxDA = dret::DocListIdxBruteBackend<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    dret::rmq::GetDocDA<TStorage, dret::Alphabet<>::int_width>>;

template <typename TStorage>
using IdxRLCSA = dret::DocListIdxBruteBackend<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    dret::rmq::GetDocRLCSA<TStorage, dret::Alphabet<>::int_width>>;

// A backend that reads nothing: brute force over it is the r-index count alone, the step every index
// that counts with sri::RIndexCount starts with (all but Brute-RI and Brute-SRI), timed in isolation.
// The width of each range goes to a sink, so the count cannot be optimized away; nothing is reported.
template <typename TStorage, uint8_t t_width>
class GetDocNone : public dret::IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = dret::IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit GetDocNone(const TStorage& t_storage) : Base(t_storage) {}

  GetDocNone() = default;

  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& /*unused*/) const {
    sink_ += e - b;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return 0;
  }

  dret::SizeReport GetSizeReport() const { return {}; }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& /*unused*/, const dret::JSON& /*unused*/) override {}

  mutable std::size_t sink_ = 0;
};

template <typename TStorage, uint8_t t_width>
void construct(GetDocNone<TStorage, t_width>& /*unused*/, dret::Config& /*unused*/) {}

template <typename TStorage>
using IdxCount = dret::DocListIdxBruteBackend<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    GetDocNone<TStorage, dret::Alphabet<>::int_width>>;

template <typename TIndex, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeBackend(TStorage t_storage, dret::Config& t_config) {
  auto idx = std::make_shared<TIndex>(t_storage);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

}  // namespace bench::factories::brute
