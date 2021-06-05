//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/2/21.
//

#ifndef DRET_DOC_LIST_INDEX_H_
#define DRET_DOC_LIST_INDEX_H_

#include <string>
#include <functional>
#include <memory>

namespace dret {

class DocListIndex {
 public:
  using TPattern = std::string;
  using TDocId = std::size_t;

  virtual void Search(const TPattern &t_pattern, const std::function<void(TDocId)> &t_report) = 0;
};

template<typename TLocate, typename TGetDoc>
class DocListIndexBrute : public DocListIndex {
 public:
  DocListIndexBrute(const TLocate &t_locate, const TGetDoc &t_get_doc) : locate_{t_locate}, get_doc_{t_get_doc} {}

  void Search(const TPattern &t_pattern, const std::function<void(TDocId)> &t_report) override {
    auto occurrences = locate_(t_pattern);

    for (const auto &item : occurrences) {
      t_report(get_doc_(item));
    }
  }

 private:
  TLocate locate_;
  TGetDoc get_doc_;
};

template<typename TComputeSARange, typename TComputeDocs>
class DocListIndexBasicScheme : public DocListIndex {
 public:
  DocListIndexBasicScheme(const TComputeSARange &t_csa,
                          const TComputeDocs &t_compute_docs)
      : compute_sa_range_{t_csa}, compute_docs_{t_compute_docs} {
  }

  void Search(const TPattern &t_pattern, const std::function<void(TDocId)> &t_report) override {
    auto[sp, ep] = compute_sa_range_(t_pattern);

    compute_docs_(sp, ep, t_report);
  }

 private:
  TComputeSARange compute_sa_range_;
  TComputeDocs compute_docs_;
};

}

#endif //DRET_DOC_LIST_INDEX_H_
