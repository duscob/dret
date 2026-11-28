//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/23/20.
//

#ifndef DRET_BENCHMARK_FACTORY_H_
#define DRET_BENCHMARK_FACTORY_H_

#include <cstddef>
#include <memory>

#include <sdsl/config.hpp>
#include <sdsl/bit_vectors.hpp>

#include <r_index.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>
#include <grammar/sampled_slp.h>

#include <dret/doc_freq/doc_freq_base.h>
#include <dret/doc_freq/doc_freq_gcda.h>

#include "../../tool/definitions.h"

#include "../csa_wrapper.h"

#include "doc_isas_wrapper.h"
#include "sa_wrapper.h"

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

  enum class IndexEnum {
    GCDA_ISAs = 0,
    GCDA_WT_DA
  };

  Factory(const sdsl::cache_config &_config) : config_{_config} {
    std::ofstream breakdown_csv("breakdown.csv");
    breakdown_csv << "R-Index-Basic,R-Index,Doc-Endings,DSA,Doc-DISAs,WT-DA,GCDA"
                  << std::endl;

    // Loading CSAs
    Load(r_idx_, KEY_R_INDEX, config_, "R-Index");
    size_r_idx_ = sdsl::size_in_bytes(r_idx_);
    size_r_idx_basic_ = r_idx_.size_in_bytes_basic();
    breakdown_csv << size_r_idx_basic_ << "," << size_r_idx_ << ",";

    csa_wrappers_.emplace_back(new RIndexWrapper<decltype(r_idx_)>(&r_idx_));

    seq_size_ = r_idx_.text_size();


    // Loading document ending marks
    Load(doc_endings_, KEY_DOC_END, config_, "Document Endings");
    doc_endings_rank_ = BitVectorCompactRank(&doc_endings_);
    doc_endings_select_ = BitVectorCompactSelect(&doc_endings_);
    size_doc_endings_ = sdsl::size_in_bytes(doc_endings_) + sdsl::size_in_bytes(doc_endings_rank_)
        + sdsl::size_in_bytes(doc_endings_select_);

    breakdown_csv << size_doc_endings_ << ",";

    n_doc_ = doc_endings_rank_(doc_endings_.size());
    doc_marked_[0].resize(n_doc_, false);
    doc_marked_[1].resize(n_doc_, false);

    unmark_ = decltype(unmark_)(doc_marked_);


    // Loading differential suffix array
    dsa_.size_in_bytes_ = LoadDifferentialSLP(
        KEY_DSA_RAW,
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

    breakdown_csv << dsa_.size_in_bytes_ << ",";

    sa_wrappers_.emplace_back(
        new DSAWrapper<std::remove_pointer<decltype(dsa_.dslp_)>::type, BitVectorCompactRank>(
            dsa_.dslp_,
            &doc_endings_rank_));


    // Document inverse suffix arrays
    // Loading differential document inverse suffix arrays
    doc_disas_.size_in_bytes_ = LoadDifferentialSLP(
        KEY_DOC_DISAS_RAW,
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

    breakdown_csv << doc_disas_.size_in_bytes_ << ",";

    doc_isas_wrappers_.emplace_back(new DocDISAsWrapper<std::remove_pointer<decltype(doc_disas_.dslp_)>::type>(
        doc_disas_.dslp_));

    compute_doc_freq_suff_wrappers_.emplace_back(
        new dret::ComputeDocFrequencyWithSuffixesDocISAs<
            DocISAsWrapper, decltype(doc_endings_rank_), decltype(unmark_)>(
            doc_isas_wrappers_[0].get(), &doc_endings_rank_, &unmark_)
    );


    // Loading wavelet tree on document array
//    Load(wt_da_huff_, KEY_DA_WT + std::string("_huff"), config_, "WT DA");
//    size_wt_da_huff_ = sdsl::size_in_bytes(wt_da_huff_);
//    sa_wrappers_.emplace_back(new WTOnDAWrapper<sdsl::wt_huff_int<>>(&wt_da_huff_));

    Load(wt_da_ap_, KEY_DA_WT + std::string("_ap"), config_, "WT DA");
    size_wt_da_ap_ = sdsl::size_in_bytes(wt_da_ap_);
    sa_wrappers_.emplace_back(new WTOnDAWrapper<sdsl::wt_ap<>>(&wt_da_ap_));

    breakdown_csv << size_wt_da_ap_ << ",";

    compute_doc_freq_suff_wrappers_.emplace_back(
        new dret::ComputeDocFrequencyWithSuffixesOccs<decltype(unmark_)>(&unmark_));



    // Loading GCDA components
    auto bit_compress = [](sdsl::int_vector<> &_v) { sdsl::util::bit_compress(_v); };
    {
      grammar::SLP<> tmp_slp;
      std::vector<std::size_t> tmp_roots;
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(tmp_slp);

      auto report_compact_seq = [&tmp_roots](const auto &_var) {
        tmp_roots.emplace_back(_var);
      };

      re_pair_reader.Read(cache_file_name(KEY_DA_RAW, config_), slp_wrapper, report_compact_seq);

      gcda_slp_ = decltype(gcda_slp_)(tmp_slp, bit_compress, bit_compress);

      grammar::Construct(gcda_roots_, tmp_roots);
      bit_compress(gcda_roots_);
    }
    size_gcda_ = sdsl::size_in_bytes(gcda_slp_) + sdsl::size_in_bytes(gcda_roots_);

    {
      sdsl::bit_vector tmp_bv;
      Load(tmp_bv, KEY_GCDA_CSEQ_HEADS, config_, "GCDA Root Heads");

      gcda_root_head_bv_ = decltype(gcda_root_head_bv_)(tmp_bv);
      gcda_root_head_bv_rank_.set_vector(&gcda_root_head_bv_);
      gcda_root_head_bv_select_.set_vector(&gcda_root_head_bv_);
    }
    size_gcda_ += sdsl::size_in_bytes(gcda_root_head_bv_)
        + sdsl::size_in_bytes(gcda_root_head_bv_rank_)
        + sdsl::size_in_bytes(gcda_root_head_bv_select_);

    Load(gcda_occs_bvs_.first, KEY_GCDA_FIRST_OCCS, config_, "GCDA First Occs BVs");
    Load(gcda_occs_bvs_.second, KEY_GCDA_LAST_OCCS, config_, "GCDA Last Occs BVs");
    size_gcda_ += sdsl::size_in_bytes(gcda_occs_bvs_.first) + sdsl::size_in_bytes(gcda_occs_bvs_.second);

    breakdown_csv << size_gcda_ << std::endl;
  }

  auto SequenceSize() const {
    return seq_size_;
  }

  struct Config {
    IndexEnum index;
  };

  std::pair<dret::DocFreqIndex *, std::size_t> Build(const Config &_config) const {
    switch (_config.index) {
      case IndexEnum::GCDA_ISAs: {
        return {
            dret::MakePtrDocFreqIndexBasicScheme(
                csa_wrappers_[0],
                std::shared_ptr<dret::ComputeSuffixesByDocFunctor>(
                    dret::MakePtrComputeSuffixesByDocGCDAFunctor(gcda_slp_,
                                                                 gcda_roots_,
                                                                 gcda_root_head_bv_,
                                                                 gcda_root_head_bv_rank_,
                                                                 gcda_root_head_bv_select_,
                                                                 gcda_occs_bvs_.first,
                                                                 gcda_occs_bvs_.second,
                                                                 *sa_wrappers_[0])),
                compute_doc_freq_suff_wrappers_[0]),
            size_r_idx_basic_ + size_doc_endings_ + dsa_.size_in_bytes_ + doc_disas_.size_in_bytes_ + size_gcda_
        };
      }

      case IndexEnum::GCDA_WT_DA: {
        return {
            dret::MakePtrDocFreqIndexBasicScheme(
                csa_wrappers_[0],
                std::shared_ptr<dret::ComputeSuffixesByDocFunctor>(
                    dret::MakePtrComputeSuffixesByDocGCDAFunctor(gcda_slp_,
                                                                 gcda_roots_,
                                                                 gcda_root_head_bv_,
                                                                 gcda_root_head_bv_rank_,
                                                                 gcda_root_head_bv_select_,
                                                                 gcda_occs_bvs_.first,
                                                                 gcda_occs_bvs_.second,
                                                                 *sa_wrappers_[1])),
                compute_doc_freq_suff_wrappers_[1]),
            size_r_idx_basic_ + size_doc_endings_ + size_wt_da_ap_ + size_gcda_
        };
      }
    }
    exit(4);
  }

 private:
  sdsl::cache_config config_;

  std::size_t seq_size_;

  // CSA
  std::vector<std::shared_ptr<CSAWrapper>> csa_wrappers_;
  ri::r_index<> r_idx_;
  std::size_t size_r_idx_;
  std::size_t size_r_idx_basic_;

  // Document ending marks
  BitVectorCompact doc_endings_;
  BitVectorCompactRank doc_endings_rank_;
  BitVectorCompactSelect doc_endings_select_;
  std::size_t size_doc_endings_;

  // Documents
  std::size_t n_doc_;

  using MarkedDocs = std::vector<bool>;
  MarkedDocs doc_marked_[2];
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

    std::size_t size_in_bytes_;
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

  // Wavetree on document array
//  sdsl::wt_huff_int<> wt_da_huff_;
//  std::size_t size_wt_da_huff_;

  sdsl::wt_ap<> wt_da_ap_;
  std::size_t size_wt_da_ap_;


  // GCDA components
  grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>> gcda_slp_;
  sdsl::int_vector<> gcda_roots_;

  BitVectorCompact gcda_root_head_bv_;
  BitVectorCompactRank gcda_root_head_bv_rank_;
  BitVectorCompactSelect gcda_root_head_bv_select_;

  using OccsBV = std::vector<sdsl::bit_vector>;
  std::pair<OccsBV, OccsBV> gcda_occs_bvs_;
  std::size_t size_gcda_;


  template<typename SLP, typename Roots, typename SpanSums, typename Samples, typename SampleRootsPos, typename SamplePos, typename SamplePosRank, typename SamplePosSelect, typename DSLP>
  std::size_t LoadDifferentialSLP(const std::string &_key,
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

    std::size_t size_in_bytes = 0;
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
    size_in_bytes += sdsl::size_in_bytes(_slp) + sdsl::size_in_bytes(_roots);

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
    size_in_bytes += sizeof(_seq_size) + sizeof(_seq_diff_base);

    Load(_span_sums, KEY_GRM_SPAN_SUMS + ("_" + _key), config_, "DSA Span Sums");
    size_in_bytes += sdsl::size_in_bytes(_span_sums);

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
    size_in_bytes += sizeof(_sums_diff_base);

    Load(_samples, KEY_GRM_SAMPLE_VALUES + ("_" + _key), config_, "DSA Samples");
    size_in_bytes += sdsl::size_in_bytes(_samples);

    Load(_sample_roots_pos, KEY_GRM_SAMPLE_ROOTS_POSITIONS + ("_" + _key), config_, "DSA Samples Roots pos");
    size_in_bytes += sdsl::size_in_bytes(_sample_roots_pos);

    Load(_sample_pos, KEY_GRM_SAMPLE_POSITIONS + ("_" + _key), config_, "DSA Samples Pos");
    _sample_pos_rank = SamplePosRank(&_sample_pos);
    _sample_pos_select = SamplePosSelect(&_sample_pos);
    size_in_bytes += sdsl::size_in_bytes(_sample_pos) + sdsl::size_in_bytes(_sample_pos_rank)
        + sdsl::size_in_bytes(_sample_pos_select);

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

    return size_in_bytes;
  }
};

#endif //DRET_BENCHMARK_FACTORY_H_
