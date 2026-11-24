//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Differential test of the CILCP listing traversal of Cobas, Makinen and Rossi
// (SPIRE 2020, "CMR20"), ported line by line from its implementation (dret_copy,
// 2020-05-26: tool/build_items_ilcp.cpp for the partition, and
// include/dret/doc_freq_index_rmq.h for the traversal), against brute force.
//
// The question it answers (journal, Section "Document-Aware Compression of the
// Interleaved LCP"): CMR20's traversal stops by MARKING on the merged
// (document-aware) partition, and it isolates the run that straddles the range
// end as a subrange of its own, but still gates that run by the marks. The
// journal's Proposition prop:marking-interior says marking is unsound on the
// merged partition only at that run, and repairs it by expanding the run
// unconditionally. Is CMR20's isolation enough?
//
// Variants, each run in both directions CMR20 uses (LEFTMOST over the backward
// ILCP, RIGHTMOST over the forward one); every pass alone must list every
// document:
//   cmr20        CMR20 exactly: merged partition, marking, boundary run isolated.
//   cmr20-nosplit  the same without the isolation (control: should FAIL, or the
//                test has no teeth).
//   ilcp         CMR20's plain ILCP core (unmerged runs, marking): should pass.
//
// Collections follow cilcp_random_diff_test.cpp: 3-7 documents of 1-10 letters
// over {A, B, C} (one-character documents are what make runs straddle a query
// boundary), every pattern of length 1-4. Seed and count via DRET_DIFF_SEED /
// DRET_DIFF_COLLECTIONS.
//
// Standalone (its own main, no gtest), registered with ctest; it exits non-zero
// only if cmr20 or ilcp gives a wrong answer (cmr20-nosplit is the control and
// is expected to fail). sdsl here is header-only, so it also builds alone:
//   g++ -std=c++17 -O2 -I cmake-build-release/_deps/sdsl-src/include
//       test/cmr20_traversal_diff_test.cpp -o /tmp/cmr20_diff && /tmp/cmr20_diff
//

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/rmq_succinct_sct.hpp>

namespace {

using RangeMinQuery = sdsl::rmq_succinct_sct<true>;  // CMR20's type (leftmost minimum)
enum Side { LEFTMOST = 0, RIGHTMOST = 1 };

std::size_t EnvOr(const char* name, std::size_t fallback) {
  if (const char* v = std::getenv(name)) {
    const auto n = std::strtoull(v, nullptr, 10);
    if (n > 0) return static_cast<std::size_t>(n);
  }
  return fallback;
}

std::vector<std::string> RandomCollection(std::mt19937& gen) {
  std::uniform_int_distribution<std::size_t> n_docs(3, 7);
  std::uniform_int_distribution<std::size_t> doc_len(1, 10);
  std::uniform_int_distribution<int> letter(0, 2);
  std::vector<std::string> docs(n_docs(gen));
  for (auto& doc : docs) {
    doc.resize(doc_len(gen));
    for (auto& c : doc) c = static_cast<char>('A' + letter(gen));
  }
  return docs;
}

std::vector<std::string> AllPatterns(const std::vector<std::string>& docs) {
  std::set<std::string> pats;
  for (const auto& doc : docs)
    for (std::size_t i = 0; i < doc.size(); ++i)
      for (std::size_t m = 1; m <= 4 && i + m <= doc.size(); ++m) pats.insert(doc.substr(i, m));
  return {pats.begin(), pats.end()};
}

std::vector<std::size_t> SuffixArray(const std::string& t) {
  std::vector<std::size_t> sa(t.size());
  for (std::size_t i = 0; i < sa.size(); ++i) sa[i] = i;
  std::sort(sa.begin(), sa.end(), [&t](std::size_t a, std::size_t b) {
    return t.compare(a, std::string::npos, t, b, std::string::npos) < 0;
  });
  return sa;
}

// sdsl's construct_lcp_kasai convention: lcp[0] = 0, lcp[k] = lcp(SA[k-1], SA[k]).
std::vector<std::size_t> Lcp(const std::string& t, const std::vector<std::size_t>& sa) {
  std::vector<std::size_t> lcp(sa.size(), 0);
  for (std::size_t k = 1; k < sa.size(); ++k) {
    std::size_t l = 0;
    while (sa[k - 1] + l < t.size() && sa[k] + l < t.size() && t[sa[k - 1] + l] == t[sa[k] + l]) ++l;
    lcp[k] = l;
  }
  return lcp;
}

// One direction of a run partition: run heads, rank/select over them, RMQ over the run values.
struct Runs {
  std::vector<std::size_t> heads;  // positions of the run heads, increasing
  std::vector<std::size_t> values;
  RangeMinQuery rmq;
  std::size_t n = 0;
  std::size_t rank(std::size_t i) const {  // run heads in [0, i)
    return std::lower_bound(heads.begin(), heads.end(), i) - heads.begin();
  }
  std::size_t select(std::size_t k) const {  // k-th head, 1-based; one past the last is n
    return k <= heads.size() ? heads[k - 1] : n;
  }
};

struct Index {
  std::size_t n_docs = 0;
  std::vector<std::size_t> sa, da;
  std::vector<std::size_t> ilcp[2];
  Runs ilcp_runs[2], cilcp_runs[2];
  std::string text;
};

// build_items_ilcp.cpp, ILCP runs: a new run where the value changes.
Runs IlcpRuns(const std::vector<std::size_t>& ilcp) {
  Runs r;
  r.n = ilcp.size();
  r.values.push_back(ilcp[0]);
  r.heads.push_back(0);
  for (std::size_t i = 1; i < ilcp.size(); ++i)
    if (ilcp[i - 1] != ilcp[i]) { r.values.push_back(ilcp[i]); r.heads.push_back(i); }
  r.rmq = RangeMinQuery(&r.values);
  return r;
}

// build_items_ilcp.cpp, CILCP runs, verbatim (including how the outer ++i interacts with the inner loops).
Runs CilcpRuns(const std::vector<std::size_t>& ilcp, const std::vector<std::size_t>& da) {
  Runs r;
  r.n = ilcp.size();
  r.values.push_back(ilcp[0]);
  r.heads.push_back(0);
  for (std::size_t i = 1; i < ilcp.size(); ++i) {
    std::size_t l = ilcp[i - 1];
    std::size_t d = da[i - 1];
    if (l == ilcp[i]) {
      while (++i < ilcp.size() && l == ilcp[i]) {}
    } else if (d == da[i]) {
      do {
        l = std::min(l, ilcp[i]);
      } while (++i < ilcp.size() && d == da[i]);
    }
    if (i < ilcp.size()) {
      r.values[r.values.size() - 1] = l;
      r.values.push_back(ilcp[i]);
      r.heads.push_back(i);
    }
  }
  r.rmq = RangeMinQuery(&r.values);
  return r;
}

Index Build(const std::vector<std::string>& docs) {
  Index x;
  x.n_docs = docs.size();
  for (const auto& d : docs) { x.text += d; x.text += '\1'; }
  x.text += '\0';
  x.sa = SuffixArray(x.text);
  const std::size_t n = x.text.size();

  std::vector<std::size_t> doc_of(n);
  for (std::size_t p = 0, d = 0; p < n; ++p) {
    doc_of[p] = d < docs.size() ? d : docs.size();  // the final '\0' gets a dummy document
    if (x.text[p] == '\1') ++d;
  }
  x.da.resize(n);
  for (std::size_t i = 0; i < n; ++i) x.da[i] = doc_of[x.sa[i]];

  // Per-document LCP over doc + '\0', as build_items_ilcp.cpp builds it.
  std::vector<std::vector<std::size_t>> lcp_docs, sa_docs;
  for (const auto& d : docs) {
    std::string t = d + '\0';
    sa_docs.push_back(SuffixArray(t));
    lcp_docs.push_back(Lcp(t, sa_docs.back()));
  }
  x.ilcp[0].assign(n, 0);
  x.ilcp[1].assign(n, 0);
  std::vector<std::size_t> idx(docs.size(), 0);
  for (std::size_t i = 0; i < n; ++i) {
    const auto d = x.da[i];
    if (d == docs.size()) continue;  // the dummy '\0' suffix: value 0
    auto& k = idx[d];
    const auto& l = lcp_docs[d];
    x.ilcp[0][i] = l[k];
    x.ilcp[1][i] = l[++k % l.size()];
  }
  for (std::size_t d = 0; d < docs.size(); ++d)  // the global order of a document's suffixes is its own
    if (idx[d] != lcp_docs[d].size()) { std::cerr << "suffix count mismatch\n"; std::exit(2); }

  for (int dir = 0; dir < 2; ++dir) {
    x.ilcp_runs[dir] = IlcpRuns(x.ilcp[dir]);
    x.cilcp_runs[dir] = CilcpRuns(x.ilcp[dir], x.da);
  }
  return x;
}

enum class Variant { kCmr20, kCmr20NoSplit, kIlcp };
const char* Name(Variant v) {
  return v == Variant::kCmr20 ? "cmr20" : v == Variant::kCmr20NoSplit ? "cmr20-nosplit" : "ilcp";
}

// GetExtremeOccurrencesRMQ + RMQAlgoCore{ILCP,CILCP} + GetValuesRMQ{ILCP,CILCP}Functor +
// RMQ{ILCP,CILCP}Reporter, one side. Returns the distinct documents reported.
std::vector<std::size_t> Traverse(const Index& x, Variant v, Side side, std::size_t sp, std::size_t ep) {
  const int dir = side;
  const bool merged = v != Variant::kIlcp;
  const Runs& R = merged ? x.cilcp_runs[dir] : x.ilcp_runs[dir];
  std::vector<char> marked(x.n_docs + 1, 0);
  std::set<std::size_t> out;

  using Range = std::pair<long long, long long>;
  std::stack<Range> stack;
  // set_first_ranges
  const long long a = static_cast<long long>(R.rank(sp + 1)) - 1;
  const long long b = static_cast<long long>(R.rank(ep + 1)) - 1;
  if (v == Variant::kCmr20 && a < b) {
    if (side == LEFTMOST) { stack.emplace(b, b); stack.emplace(a, b - 1); }
    else { stack.emplace(a, a); stack.emplace(a + 1, b); }
  } else {
    stack.emplace(a, b);
  }

  // GetValues: the position a run stands for.
  auto value_pos = [&](std::size_t i) -> std::size_t {
    if (!merged || dir == 0) return std::max(sp, R.select(i + 1));  // the ILCP functor ignores the side
    return std::min(ep, R.select(i + 2) - 1);
  };

  while (!stack.empty()) {
    const auto [rsp, rep] = stack.top();
    stack.pop();
    if (rsp > rep) continue;
    const std::size_t i = R.rmq(static_cast<std::size_t>(rsp), static_cast<std::size_t>(rep));
    const std::size_t pos = value_pos(i);
    const std::size_t doc = x.da[pos];
    if (marked[doc]) continue;  // is_reported: the marks

    // report: the selected position, then the rest of the run inside the range.
    out.insert(doc);
    marked[doc] = 1;
    if (!merged) {
      const std::size_t bb = std::max(sp, R.select(i + 1)) + 1;
      const std::size_t ee = std::min(ep, R.select(i + 2) - 1);
      for (std::size_t k = bb; k <= ee; ++k) { out.insert(x.da[k]); marked[x.da[k]] = 1; }
    } else {
      const long long bb = static_cast<long long>(std::max(sp, R.select(i + 1))) + 1 - dir;
      const long long ee = static_cast<long long>(std::min(ep, R.select(i + 2) - 1)) - dir;
      for (long long k = bb; k <= ee; ++k) {
        const auto d2 = x.da[static_cast<std::size_t>(k)];
        if (doc == d2) break;
        out.insert(d2);
        marked[d2] = 1;
      }
    }

    const long long idx = static_cast<long long>(i);
    if (side == LEFTMOST) stack.emplace(idx + 1, rep);
    stack.emplace(rsp, idx - 1);
    if (side == RIGHTMOST) stack.emplace(idx + 1, rep);
  }
  out.erase(x.n_docs);  // never reached by a pattern over {A, B, C}; kept out of the comparison anyway
  return {out.begin(), out.end()};
}

std::vector<std::size_t> BruteForce(const std::vector<std::string>& docs, const std::string& p) {
  std::set<std::size_t> s;
  for (std::size_t d = 0; d < docs.size(); ++d)
    if (docs[d].find(p) != std::string::npos) s.insert(d);
  return {s.begin(), s.end()};
}

std::string Show(const std::vector<std::size_t>& v) {
  std::string s = "{";
  for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + std::to_string(v[i]);
  return s + "}";
}

}  // namespace

int main() {
  const auto seed = EnvOr("DRET_DIFF_SEED", 20260907u);
  const auto n_collections = EnvOr("DRET_DIFF_COLLECTIONS", 700u);
  const Variant variants[] = {Variant::kCmr20, Variant::kCmr20NoSplit, Variant::kIlcp};

  int status = 0;
  for (const auto v : variants) {
    for (const Side side : {LEFTMOST, RIGHTMOST}) {
      std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));
      std::size_t n_queries = 0, n_wrong = 0, n_bad_colls = 0;
      for (std::size_t c = 0; c < n_collections; ++c) {
        const auto docs = RandomCollection(gen);
        const Index x = Build(docs);
        bool bad = false;
        for (const auto& p : AllPatterns(docs)) {
          // [sp, ep]: the suffixes that start with p.
          std::size_t sp = x.sa.size(), ep = 0;
          for (std::size_t i = 0; i < x.sa.size(); ++i)
            if (x.text.compare(x.sa[i], p.size(), p) == 0) { sp = std::min(sp, i); ep = i; }
          const auto got = Traverse(x, v, side, sp, ep);
          const auto want = BruteForce(docs, p);
          ++n_queries;
          if (got != want) {
            bad = true;
            if (++n_wrong <= 3) {
              std::cout << "  " << Name(v) << (side == LEFTMOST ? " left" : " right") << ": collection " << c
                        << ", pattern \"" << p << "\", docs [";
              for (std::size_t d = 0; d < docs.size(); ++d) std::cout << (d ? ", " : "") << '"' << docs[d] << '"';
              std::cout << "]  got " << Show(got) << " want " << Show(want) << '\n';
            }
          }
        }
        n_bad_colls += bad;
      }
      std::cout << Name(v) << (side == LEFTMOST ? " LEFTMOST " : " RIGHTMOST") << ": " << n_wrong
                << " wrong of " << n_queries << " queries, in " << n_bad_colls << " of " << n_collections
                << " collections (seed " << seed << ")\n";
      if (v != Variant::kCmr20NoSplit && n_wrong) status = 1;
    }
  }
  return status;
}
