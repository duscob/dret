//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/23/20.
//

#ifndef DRET_DOC_FREQ_INDEX_H_
#define DRET_DOC_FREQ_INDEX_H_

#include <cstddef>
#include <unordered_map>
#include <memory>

namespace dret {

class ComputeSuffixesByDocFunctor {
 public:
  virtual std::vector<std::size_t> operator()(std::size_t _sp, std::size_t _ep) const = 0;
};

class ComputeDocFrequencyWithSuffixes {
 public:
  virtual void operator()(const std::vector<std::size_t> &_suffixes,
                          const std::function<void(std::size_t, std::size_t)> &_report_doc_freq) const = 0;

  virtual std::unordered_map<std::size_t, std::size_t> operator()(const std::vector<std::size_t> &_suffixes) const {
    std::unordered_map<std::size_t, std::size_t> doc_freq;
    auto report = [&doc_freq](std::size_t doc, std::size_t freq) {
      doc_freq.emplace(doc, freq);
    };

    (*this)(_suffixes, report);

    return doc_freq;
  }
};

template<typename GetSuffixPosInDoc, typename GetDoc, typename UnmarkDoc>
class ComputeDocFrequencyWithSuffixesDocISAs : public ComputeDocFrequencyWithSuffixes {
 public:
  ComputeDocFrequencyWithSuffixesDocISAs(const GetSuffixPosInDoc *_get_suffix_pos_in_doc,
                                         const GetDoc *_get_doc,
                                         UnmarkDoc *_unmark_doc)
      : get_suffix_pos_in_doc_{_get_suffix_pos_in_doc}, get_doc_{_get_doc}, unmark_doc_{_unmark_doc} {
  }

  void operator()(const std::vector<std::size_t> &_suffixes,
                  const std::function<void(std::size_t, std::size_t)> &_report_doc_freq) const override {
    for (int i = 0; i < _suffixes.size(); i += 2) {
      auto suffix_1 = _suffixes[i];
      auto suffix_2 = _suffixes[i + 1];
      auto doc = (*get_doc_)(suffix_1);

      if (suffix_1 == suffix_2) { // Pattern occurs exactly once
        _report_doc_freq(doc, 1);
      } else { // Pattern occurs more than once
        auto values = (*get_suffix_pos_in_doc_)(doc, suffix_1, suffix_2);
        auto &doc_sp = values.first;
        auto &doc_ep = values.second;

        if (doc_sp > doc_ep) {
          std::swap(doc_sp, doc_ep);
        }

        _report_doc_freq(doc, doc_ep - doc_sp + 1);
      }

      (*unmark_doc_)(doc);
    }
  }

 private:
  const GetSuffixPosInDoc *get_suffix_pos_in_doc_;
  const GetDoc *get_doc_;
  UnmarkDoc *unmark_doc_;
};

//***********
// Document Frequency Index
//***********
class DocFreqIndex {
 public:
  virtual std::unordered_map<std::size_t, std::size_t> Search(const std::string &_pattern) const = 0;
};

template<typename CSA, typename GetSuffixesByDoc, typename GetDocFreq>
class DocFreqIndexBasicScheme : public DocFreqIndex {
 public:
  DocFreqIndexBasicScheme(const std::shared_ptr<CSA> &_csa,
                          const std::shared_ptr<GetSuffixesByDoc> &_get_suffixes_by_doc,
                          const std::shared_ptr<GetDocFreq> &_get_doc_freq) :
      csa_{_csa}, get_suffixes_by_doc_{_get_suffixes_by_doc}, get_doc_freq_{_get_doc_freq} {
  }

  std::unordered_map<std::size_t, std::size_t> Search(const std::string &_pattern) const override {
    auto values = csa_->Search(_pattern);
    auto &sp = values.first;
    auto &ep = values.second;

    auto suffixes = (*get_suffixes_by_doc_)(sp, ep);

    return (*get_doc_freq_)(suffixes);
  }

 private:
  std::shared_ptr<CSA> csa_;
  std::shared_ptr<GetSuffixesByDoc> get_suffixes_by_doc_;
  std::shared_ptr<GetDocFreq> get_doc_freq_;
};

template<typename CSA, typename GetSuffixesByDoc, typename GetDocFreq>
auto MakeDocFreqIndexBasicScheme(const std::shared_ptr<CSA> &_csa,
                                 const std::shared_ptr<GetSuffixesByDoc> &_get_suffixes_by_doc,
                                 const std::shared_ptr<GetDocFreq> &_get_doc_freq) {
  return DocFreqIndexBasicScheme<CSA, GetSuffixesByDoc, GetDocFreq>(_csa, _get_suffixes_by_doc, _get_doc_freq);
}

template<typename CSA, typename GetSuffixesByDoc, typename GetDocFreq>
auto MakePtrDocFreqIndexBasicScheme(const std::shared_ptr<CSA> &_csa,
                                 const std::shared_ptr<GetSuffixesByDoc> &_get_suffixes_by_doc,
                                 const std::shared_ptr<GetDocFreq> &_get_doc_freq) {
  return new DocFreqIndexBasicScheme<CSA, GetSuffixesByDoc, GetDocFreq>(_csa, _get_suffixes_by_doc, _get_doc_freq);
}

}

#endif //DRET_DOC_FREQ_INDEX_H_
