//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/23/20.
//

#ifndef DRET_BENCHMARK_FACTORY_H_
#define DRET_BENCHMARK_FACTORY_H_

#include <cstddef>
#include <memory>

#include <sdsl/config.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>

#include <rindex/r_index.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include <doc_freq_index.h>
#include <doc_freq_index_brute.h>
#include <doc_freq_index_rmq.h>

#include "../tool/definitions.h"

#include "csa_wrapper.h"
#include "sa_wrapper.h"
#include "doc_isas_wrapper.h"

const char *KEY_GRM_SPAN_SUMS = "grm_sums";
const char *KEY_GRM_SAMPLE_VALUES = "grm_sample_val";
const char *KEY_GRM_SAMPLE_POSITIONS = "grm_sample_pos";
const char *KEY_GRM_SAMPLE_ROOTS_POSITIONS = "grm_sample_root_pos";

template<typename T>
auto Load(T &_v, const std::string &_key, const sdsl::cache_config &_config, const std::string &_msg) {
  if (!sdsl::load_from_cache(_v, _key, _config)) {
    std::cerr << "ERROR: Failed to load " + _msg + " (" + cache_file_name(_key, _config) + ")" << std::endl;
    exit(1);
  }
}

template<typename BitVectors>
class IsMarkedWrapper {
 public:
  IsMarkedWrapper(const BitVectors *_marked = nullptr) : bit_vectors_{_marked} {}

  template<typename SuffixDoc>
  auto operator()(std::size_t /*_i*/, const SuffixDoc &_suffix_doc, dret::OccurrenceSide _side) const {
    return bit_vectors_[_side == dret::OccurrenceSide::LEFTMOST ? 0 : 1][_suffix_doc.second];
  }

 private:
  const BitVectors *bit_vectors_;
};

template<typename BitVectors>
class MarkWrapper {
 public:
  MarkWrapper(BitVectors *_marked = nullptr) : bit_vectors_{_marked} {}

  void operator()(std::size_t _doc, dret::OccurrenceSide _side) {
    bit_vectors_[_side == dret::OccurrenceSide::LEFTMOST ? 0 : 1][_doc] = 1;
  }

 private:
  BitVectors *bit_vectors_;
};

template<typename BitVectors>
class UnmarkWrapper {
 public:
  UnmarkWrapper(BitVectors *_marked = nullptr) : bit_vectors_{_marked} {}

  void operator()(std::size_t _doc) {
    bit_vectors_[0][_doc] = 0;
    bit_vectors_[1][_doc] = 0;
  }

 private:
  BitVectors *bit_vectors_;
};

class Factory {
 public:
  using BitVectorCompact = sdsl::sd_vector<>;
  using BitVectorCompactRank = BitVectorCompact::rank_1_type;
  using BitVectorCompactSelect = BitVectorCompact::select_1_type;

  using RangeMinQuery = sdsl::rmq_succinct_sct<true>;
  using RangeMaxQuery = sdsl::rmq_succinct_sct<false>;

  enum class IndexEnum { Brute = 0, SADA, ILCP, CILCP };

  Factory(const sdsl::cache_config &_config) : config_{_config} {
    // Loading CSAs
    Load(r_idx_, KEY_R_INDEX, config_, "R-Index");
    csa_wrappers_.emplace_back(new RIndexWrapper<decltype(r_idx_)>(&r_idx_));

    // Loading document ending marks
    Load(doc_endings_, KEY_DOC_END, config_, "Document Endings");
    doc_endings_rank_ = BitVectorCompactRank(&doc_endings_);
    doc_endings_select_ = BitVectorCompactSelect(&doc_endings_);

    n_doc_ = doc_endings_rank_(doc_endings_.size());
    doc_marked_[0].resize(n_doc_, false);
    doc_marked_[1].resize(n_doc_, false);

    is_marked_ = decltype(is_marked_)(doc_marked_);
    mark_ = decltype(mark_)(doc_marked_);
    unmark_ = decltype(unmark_)(doc_marked_);

    // Loading differential suffix array
    LoadDifferentialSLP(KEY_DSA_RAW,
                        dsa_.slp_,
                        dsa_.roots_,
                        dsa_.seq_size_,
                        dsa_.seq_diff_base_,
                        dsa_.span_sums_,
                        dsa_.sums_diff_base_,
                        dsa_.samples_,
                        dsa_.sample_roots_pos_,
                        dsa_.samples_pos_,
                        dsa_.samples_pos_rank_,
                        dsa_.samples_pos_select_,
                        dsa_.dslp_);

    sa_wrappers_.emplace_back(
        new DSAWrapper<std::remove_pointer<decltype(dsa_.dslp_)>::type, BitVectorCompactRank>(
            dsa_.dslp_,
            &doc_endings_rank_));

    // Document inverse suffix arrays
    // Loading differential document inverse suffix arrays
    LoadDifferentialSLP(KEY_DOC_DISAS_RAW,
                        doc_disas_.slp_,
                        doc_disas_.roots_,
                        doc_disas_.seq_size_,
                        doc_disas_.seq_diff_base_,
                        doc_disas_.span_sums_,
                        doc_disas_.sums_diff_base_,
                        doc_disas_.samples_,
                        doc_disas_.sample_roots_pos_,
                        doc_disas_.samples_pos_,
                        doc_disas_.samples_pos_rank_,
                        doc_disas_.samples_pos_select_,
                        doc_disas_.dslp_);

    doc_isas_wrappers_.emplace_back(new DocDISAsWrapper<std::remove_pointer<decltype(doc_disas_.dslp_)>::type>(
        doc_disas_.dslp_));

    compute_doc_freq_suff_wrappers_.emplace_back(
        new dret::ComputeDocFrequencyWithSuffixesDocISAs<
            DocISAsWrapper, decltype(doc_endings_rank_), decltype(unmark_)>(
            doc_isas_wrappers_[0].get(), &doc_endings_rank_, &unmark_)
    );


    // Loading Sada components
    Load(core_sada_.left_rmq, KEY_SADA_RMINQ, config_, "Sada RMinQ");
    Load(core_sada_.right_rmq, KEY_SADA_RMAXQ, config_, "Sada RMaxQ");


    // Loading ILCP components
    Load(core_ilcp_.left_rmq, KEY_ILCP_BACKWARD_RMQ, config_, "ILCP Left RMQ");
    Load(core_ilcp_.right_rmq, KEY_ILCP_FORWARD_RMQ, config_, "ILCP Right RMQ");
    {
      sdsl::bit_vector tmp_bv;
      std::string keys[2] = {KEY_ILCP_BACKWARD_RUN_HEADS, KEY_ILCP_FORWARD_RUN_HEADS};
      for (std::size_t i = 0; i < 2; ++i) {
        Load(tmp_bv, keys[i], config_, "ILCP Left/Right Run Heads");
        core_ilcp_.run_heads[i] =
            std::remove_reference<std::remove_pointer<decltype(core_ilcp_.run_heads[i])>::type>::type(tmp_bv);
        core_ilcp_.run_heads_rank[i].set_vector(&core_ilcp_.run_heads[i]);
        core_ilcp_.run_heads_select[i].set_vector(&core_ilcp_.run_heads[i]);
      }
    }


    // Loading CILCP components
    Load(core_cilcp_.left_rmq, KEY_CILCP_BACKWARD_RMQ, config_, "CILCP Left RMQ");
    Load(core_cilcp_.right_rmq, KEY_CILCP_FORWARD_RMQ, config_, "CILCP Right RMQ");
    {
      sdsl::bit_vector tmp_bv;
      std::string keys[2] = {KEY_CILCP_BACKWARD_RUN_HEADS, KEY_CILCP_FORWARD_RUN_HEADS};
      for (std::size_t i = 0; i < 2; ++i) {
        Load(tmp_bv, keys[i], config_, "CILCP Left/Right Run Heads");
        core_cilcp_.run_heads[i] =
            std::remove_reference<std::remove_pointer<decltype(core_cilcp_.run_heads[i])>::type>::type(tmp_bv);
        core_cilcp_.run_heads_rank[i].set_vector(&core_cilcp_.run_heads[i]);
        core_cilcp_.run_heads_select[i].set_vector(&core_cilcp_.run_heads[i]);
      }
    }

    // RMQ get values functors
    get_values_rmq_functors_.emplace_back(new dret::GetValuesRMQSadaFunctor<SAWrapper>(sa_wrappers_[0]));
    get_values_rmq_functors_.emplace_back(
        new dret::GetValuesRMQILCPFunctor<SAWrapper, decltype(core_ilcp_)>(sa_wrappers_[0], &core_ilcp_));
    get_values_rmq_functors_.emplace_back(
        new dret::GetValuesRMQCILCPFunctor<SAWrapper, decltype(core_cilcp_)>(sa_wrappers_[0], &core_cilcp_));

    // RMQ reporters
    rmq_reporters.emplace_back(new dret::RMQSadaReporter<decltype(mark_)>(&mark_));
    rmq_reporters.emplace_back(
        new dret::RMQILCPReporter<decltype(mark_), SAWrapper, decltype(core_ilcp_)>(
            &mark_, sa_wrappers_[0], &core_ilcp_));
    rmq_reporters.emplace_back(
        new dret::RMQCILCPReporter<decltype(mark_), SAWrapper, decltype(core_cilcp_)>(
            &mark_, sa_wrappers_[0], &core_cilcp_));
  }

  struct Config {
    IndexEnum index;
  };

  dret::DocFreqIndex *Build(const Config &_config) const {
    switch (_config.index) {
      case IndexEnum::Brute: {
        return dret::MakeNewDocFreqIndexBrute(*csa_wrappers_[0], doc_endings_rank_);
      }
      case IndexEnum::SADA: {
        return dret::MakePtrDocFreqIndexBasicScheme(
            csa_wrappers_[0],
            std::shared_ptr<dret::ComputeSuffixesByDocFunctor>(
                dret::MakeNewComputeSuffixesByDocRMQFunctor(core_sada_,
                                                            *get_values_rmq_functors_[0],
                                                            *rmq_reporters[0],
                                                            is_marked_)),
            compute_doc_freq_suff_wrappers_[0]
        );
      }
      case IndexEnum::ILCP: {
        return dret::MakePtrDocFreqIndexBasicScheme(
            csa_wrappers_[0],
            std::shared_ptr<dret::ComputeSuffixesByDocFunctor>(
                dret::MakeNewComputeSuffixesByDocRMQFunctor(core_ilcp_,
                                                            *get_values_rmq_functors_[1],
                                                            *rmq_reporters[1],
                                                            is_marked_)),
            compute_doc_freq_suff_wrappers_[0]
        );
      }
      case IndexEnum::CILCP: {
        return dret::MakePtrDocFreqIndexBasicScheme(
            csa_wrappers_[0],
            std::shared_ptr<dret::ComputeSuffixesByDocFunctor>(
                dret::MakeNewComputeSuffixesByDocRMQFunctor(core_cilcp_,
                                                            *get_values_rmq_functors_[2],
                                                            *rmq_reporters[2],
                                                            is_marked_)),
            compute_doc_freq_suff_wrappers_[0]
        );
      }
    }
    exit(4);
  }

 private:
  sdsl::cache_config config_;

  // CSA
  std::vector<std::shared_ptr<CSAWrapper>> csa_wrappers_;
  ri::r_index<> r_idx_;

  // Document ending marks
  BitVectorCompact doc_endings_;
  BitVectorCompactRank doc_endings_rank_;
  BitVectorCompactSelect doc_endings_select_;

  // Documents
  std::size_t n_doc_;

  using MarkedDocs = std::vector<bool>;
  MarkedDocs doc_marked_[2];
  IsMarkedWrapper<MarkedDocs> is_marked_;
  MarkWrapper<MarkedDocs> mark_;
  UnmarkWrapper<MarkedDocs> unmark_;

  struct DSLPWrapper {
    grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>> slp_;
    sdsl::int_vector<> roots_;
    std::size_t seq_size_;
    uint64_t seq_diff_base_;

    sdsl::int_vector<> span_sums_;
    uint64_t sums_diff_base_;

    sdsl::int_vector<> samples_;
    sdsl::int_vector<> sample_roots_pos_;
    BitVectorCompact samples_pos_;
    BitVectorCompactRank samples_pos_rank_;
    BitVectorCompactSelect samples_pos_select_;

    grammar::DifferentialSLPWrapper<decltype(slp_),
                                    decltype(roots_),
                                    decltype(span_sums_),
                                    decltype(samples_),
                                    decltype(sample_roots_pos_),
                                    decltype(samples_pos_)> *dslp_ = nullptr;
  };

  // Suffix array
  std::vector<SAWrapper *> sa_wrappers_;
  // Differential suffix array
  DSLPWrapper dsa_;

  // Wrapper to compute document frequency with suffixes
  std::vector<std::shared_ptr<dret::ComputeDocFrequencyWithSuffixes>> compute_doc_freq_suff_wrappers_;
  // Document inverse suffix arrays
  std::vector<std::unique_ptr<DocISAsWrapper>> doc_isas_wrappers_;
  // Differential Document Inverse Suffix Arrays
  DSLPWrapper doc_disas_;

  // Sada components
  dret::RMQAlgoCoreSada<RangeMinQuery, RangeMaxQuery> core_sada_;

  // ILCP components
  //TODO Experiment with other bitvectors
  dret::RMQAlgoCoreILCP<RangeMinQuery, RangeMinQuery, sdsl::rrr_vector<>> core_ilcp_;

  // CILCP
  dret::RMQAlgoCoreCILCP<RangeMinQuery, RangeMinQuery, sdsl::sd_vector<>> core_cilcp_;

  // RMQ get values
  std::vector<std::unique_ptr<dret::GetValuesRMQFunctor>> get_values_rmq_functors_;

  // RMQ reporters
  std::vector<std::unique_ptr<dret::RMQReporter>> rmq_reporters;

  template<typename SLP, typename Roots, typename SpanSums, typename Samples, typename SampleRootsPos, typename SamplePos, typename SamplePosRank, typename SamplePosSelect, typename DSLP>
  void LoadDifferentialSLP(const std::string &_key,
                           SLP &_slp,
                           Roots &_roots,
                           std::size_t &_seq_size,
                           uint64_t &_seq_diff_base,
                           SpanSums &_span_sums,
                           uint64_t &_sums_diff_base,
                           Samples &_samples,
                           SampleRootsPos &_sample_roots_pos,
                           SamplePos &_sample_pos,
                           SamplePosRank &_sample_pos_rank,
                           SamplePosSelect &_sample_pos_select,
                           DSLP &_dslp) const {
    std::string _basefile = cache_file_name(_key, config_);
    auto bit_compress = [](sdsl::int_vector<> &_v) { sdsl::util::bit_compress(_v); };

    {
      grammar::SLP<> tmp_slp;
      std::vector<std::size_t> tmp_roots;
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(tmp_slp);

      auto report_compact_seq = [&tmp_roots](const auto &_var) {
        tmp_roots.emplace_back(_var);
      };

      re_pair_reader.Read(_basefile, slp_wrapper, report_compact_seq);

      _slp = SLP(tmp_slp, bit_compress, bit_compress);

      grammar::Construct(_roots, tmp_roots);
      bit_compress(_roots);
    }

    {
      auto filename = _basefile + ".info";
      std::ifstream in(filename);
      if (!in) {
        std::cerr << "ERROR: Failed to load " + filename << std::endl;
        exit(1);
      }
      in >> _seq_size;

      int64_t minimal;
      in >> minimal;
      _seq_diff_base = minimal < 0 ? std::abs(minimal) : 0;
    }

    Load(_span_sums, KEY_GRM_SPAN_SUMS + ("_" + _key), config_, "DSA Span Sums");

    {
      auto filename = cache_file_name(KEY_GRM_SPAN_SUMS + ("_" + _key), config_) + ".info";
      std::ifstream in(filename);
      if (!in) {
        std::cerr << "ERROR: Failed to load " + filename << std::endl;
        exit(1);
      }
      int64_t minimal;
      in >> minimal;
      _sums_diff_base = minimal < 0 ? std::abs(minimal) : 0;
    }

    Load(_samples, KEY_GRM_SAMPLE_VALUES + ("_" + _key), config_, "DSA Samples");

    Load(_sample_roots_pos, KEY_GRM_SAMPLE_ROOTS_POSITIONS + ("_" + _key), config_, "DSA Samples Roots pos");

    Load(_sample_pos, KEY_GRM_SAMPLE_POSITIONS + ("_" + _key), config_, "DSA Samples Pos");
    _sample_pos_rank = SamplePosRank(&_sample_pos);
    _sample_pos_select = SamplePosSelect(&_sample_pos);

    _dslp = new typename std::remove_pointer<DSLP>::type(_seq_size,
                                                         _slp,
                                                         _roots,
                                                         _seq_diff_base,
                                                         _span_sums,
                                                         _sums_diff_base,
                                                         _samples,
                                                         _sample_roots_pos,
                                                         _sample_pos,
                                                         _sample_pos_rank,
                                                         _sample_pos_select);
  }
};

#endif //DRET_BENCHMARK_FACTORY_H_
