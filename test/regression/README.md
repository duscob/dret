# Regression harness

The safety net for the redesign. Every index `bm_doc_list` can build (`sweep.json`, one cell of
each family, representation, backend, codec and policy) is built from scratch on two small
collections and queried. A run is reduced to a manifest per collection:

- the md5 of each index's answers (one sorted line of document ids per pattern),
- its resident and serialized sizes,
- its query time (median of the repetitions).

`baseline/` holds the manifests recorded at the start of the redesign, on `redesign/base`.

## Collections

Both are cut, in document order, from the journal collections (`experiments/data/` is local and
not tracked; the manifests store the sha256 of each subset, so a check against different data
fails at once):

```bash
R=test/regression/regression.py
python3 $R make-data --source experiments/data/revision/data \
    --out experiments/data/regression/revision-small --max-bytes 8000000
python3 $R make-data --source experiments/data/page/data \
    --out experiments/data/regression/page-small --max-bytes 8000000 --max-doc-bytes 1000000
```

`revision-small` has 847 short documents, `page-small` 35 long ones (pages up to 1 MB). Both are
queried with `experiments/data/revision/patterns`.

## Running a check

Release build of `bm_doc_list`, then, per collection (about 7 minutes each on the laptop):

```bash
python3 test/regression/regression.py run --bin cmake-build-release/benchmark/bm_doc_list \
    --data experiments/data/regression/revision-small/data \
    --patterns experiments/data/revision/patterns \
    --check test/regression/baseline/revision-small.json
```

It fails when
1. the indexes disagree with each other (the count-only index is left out: it reports nothing),
2. a benchmark of the baseline is missing, or its answers changed,
3. a size changed (skip with `--allow-size-change` when a change is meant to move sizes),
4. the query pass created, removed or wrote a file in the index directory (query mode only loads),
5. with `--time-tolerance F`, a query got slower by more than the fraction F.

Query times are otherwise only reported, as the median and extreme ratios against the baseline.
Laptop times are noisy (frequency scaling); gate on time only on a quiet machine.

A change that renames benchmarks passes `--names NAMES.json`, a map from the baseline's names to
the new ones.

## Sizes that are not reproducible

Two fresh builds of the same code, in any directories, agree on every answer and every size.
`UNSTABLE_SIZES` in `regression.py` lists the sizes to leave out of the comparison when that does
not hold, and is empty. It had two entries, both fixed in the redesign:

- **PDL-BC** (`codec=bc`), both sizes: the vendored `vnmextract` seeded its min-hashes with the
  time, so each build mined other bicliques. Builds are seeded now (`BCCodec::MiningParams::seed`).
- **RLCSA backend**, serialized size: it serialized the path of its files, so the size depended on
  the directory the index lived in. It serializes the files' contents now, so the serialized size
  is the RLCSA's own.

## Gate for merging into `redesign/base`

- `ctest` passes (including `cilcp_random_diff_test` and `cmr20_traversal_diff_test`);
- this check passes on both collections;
- for a change to a query path: the time ratios, and on shannon a timing comparison on Page and
  Revision.
