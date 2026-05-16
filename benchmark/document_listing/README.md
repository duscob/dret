# Document-listing benchmarks

Two complementary binaries cover the doc-list index family.

| Binary | Measures | When to use |
|---|---|---|
| `bm_query_doc_list` | Query throughput per index variant; reports size breakdown via `GetSizeReport()`. Auto-builds any missing cache artefacts on first run. | Day-to-day "how fast is this index?" runs. |
| `bm_build_items`    | Construction time per index variant; emits `construction-*.html` / `construction-*.json` sdsl memory logs and `sizes-*.json` size breakdowns. | When you specifically want construction-time numbers or memory traces. |

Both binaries share the same axis enums (`bench::axes::*`) and the same parse/name helpers (`bench::axes::EnumTraits<E>` and `bench::axes::ParseCSV<E>`) — so every flag value listed below behaves the same across both.

## Quick start

```bash
# One-shot query benchmark (builds caches inline if missing):
./bm_query_doc_list \
  --data_dir=/path/to/data --data_name=wiki \
  --patterns=/path/to/queries.txt \
  --gcda_slp_variants=default,compact_bp \
  --min_block_size=256 --max_block_size=1024 \
  --min_storing_factor=2 --max_storing_factor=8 \
  --benchmark_out_format=csv --benchmark_out=results.csv
```

```bash
# If you want construction-time + memory logs as well:
./bm_build_items \
  --data=/path/to/data/wiki \
  --gcda_slp_variants=default,compact_bp \
  --min_block_size=256 --max_block_size=1024 \
  --min_storing_factor=2 --max_storing_factor=8 \
  --benchmark_out_format=json --benchmark_out=build.json
```

## CLI axes (shared by both binaries)

| Flag | Values | Meaning |
|---|---|---|
| `--rmq_get_doc_variants` | `da`, `slp`, `slp_ns`, `dslp` | RMQ family's raw-range doc lookup. |
| `--gcda_slp_variants`    | `default`, `compact_bp`, `compact_louds`, `cslp` | GCDA's TSLP choice (also reused by RMQ-SLP for cache sharing). |
| `--bare_slp_variants`    | `default`, `raw`, `dv`, `vv` | Non-sampled `grammar::SLP<>` container variants (SLP-NS family). |
| `--pdl_variants`         | `plain`, `rp`, `bc` (or empty) | PDL stored-set codec. Empty disables PDL. |
| `--pdl_get_doc_variants` | `da`, `slp`, `dslp` | PDL raw-range get-doc backing. (`slp_ns` is rejected.) |
| `--pdl_storage_policy`   | `occurrence_weighted`, `all_internal`, `leaves_only` | PDL tree-construction storage policy. |
| `--min_block_size` / `--max_block_size` | powers of 2 | Block-size sweep for the GCDA / DGCDA / PDL families. |
| `--min_storing_factor` / `--max_storing_factor` | powers of 2 | Storing-factor sweep, same families. |

## Output formats (Google Benchmark native)

`benchmark::Initialize` picks these up automatically:

- **Console** (default): the usual pretty-printed table to stdout.
- **CSV**: `--benchmark_out_format=csv --benchmark_out=path.csv`
- **JSON**: `--benchmark_out_format=json --benchmark_out=path.json`

Each output row carries the per-component size breakdown emitted by `GetSizeReport()` as named counters (e.g. `slp_basic_slp`, `slp_sets_objs`, `count_idx`), plus the totals `Size(bytes)`, `Bits_x_Symbol`, `Time_x_Pattern`, etc.

`bm_build_items` additionally writes a `sizes-<bench-name>.json` and a memory-monitor `construction-<bench-name>.{html,json}` next to the working directory for each benchmark cell.

## Where the code lives

```
benchmark/document_listing/
├── axes.h                 — six axis enums at namespace bench::axes scope
├── enum_traits.h          — EnumTraits<E>::ParseCSV / ::Name, plus ParseCSV<E>(csv) helper
├── factory.h              — Factory<>: builds + caches index instances, auto-builds missing cache artefacts
├── bm_query_doc_list.cpp  — query benchmark binary
└── bm_build_items.cpp     — construction benchmark binary
```

To add a new axis value (e.g. a new TSLP variant), add it to `axes.h`, add the parse/name handling to `enum_traits.h`'s relevant `EnumTraits<>` specialisation, then wire it into the `Factory::MakeIndex` switch in `factory.h` and into the sweep loop in whichever binary cares about it. A future phase will pull the type-alias library out of `factory.h` into per-family headers under `factories/` so adding a variant is a one-file edit.
