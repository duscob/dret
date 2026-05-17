# Document-listing benchmarks

Three binaries cover the doc-list index family.

| Binary | Measures | Input | When to use |
|---|---|---|---|
| **`bm_doc_list`** *(recommended)* | Query throughput **or** construction time (selectable per spec); reports size breakdown via `GetSizeReport()`. Auto-builds missing cache artefacts on first run. | JSON sweep spec (`--spec=path.json`) | Day-to-day runs. One declarative file describes the whole sweep; the same binary handles build + query workflows. |
| `bm_query_doc_list` | Query throughput, gflags-driven sweep. | gflags (`--gcda_slp_variants=...`, etc.) | Legacy. Still works; kept alongside while migrating to `bm_doc_list`. |
| `bm_build_items`    | Construction time, gflags-driven sweep. Emits `construction-*.{html,json}` sdsl memory logs and `sizes-*.json` size breakdowns. | gflags (same as above) | Legacy. Use when you specifically want construction-time **memory traces** (which `bm_doc_list` doesn't emit yet). |

All three binaries share the same axis enums (`bench::axes::*`) and parse/name helpers (`bench::axes::EnumTraits<E>` and `bench::axes::ParseCSV<E>`), so axis names behave the same across the JSON spec and the legacy gflags.

## `bm_doc_list` — JSON-driven sweep

See `example_sweep.json` for a complete reference. The spec is a single JSON document with four top-level keys:

```json
{
  "mode":     "query",           // "query" | "construct"
  "dataset":  { "dir": "...", "name": "...", "patterns": "...", "data_width": 8 },
  "workload": { "reps": 10, "min_time": 0.0 },
  "sweep":    [ ...family blocks... ]
}
```

Each entry in `sweep` is one family block. The driver expands each block's Cartesian product against the family's axes and registers one Google Benchmark per cell:

| Family | Axes |
|---|---|
| `gcda`   | `tslp` (`light` / `compact_bp` / `compact_louds` / `combined`), `block_size`, `storing_factor` |
| `dgcda`  | `tslp` (`default` / `otf` / `crl` / `ev` / `dv` / `vv`), `block_size`, `storing_factor` |
| `slp_ns` | `tslp` (`iv` / `raw` / `dv` / `vv`) |
| `rmq`    | `core` (`sada` / `ilcp` / `cilcp`), `get_doc` (`da` / `slp` / `slp_ns` / `dslp`), `gcda_slp`, `bare_slp`, `block_size`, `storing_factor` |
| `pdl`    | `codec` (`plain` / `rp` / `bc`), `get_doc` (`da` / `slp` / `dslp`), `policy`, `block_size`, `storing_factor` |
| `brute`  | `kind` (`r-index` / `sr-index`), `sampling_size` |

Missing axes default to the family's canonical value (`light`, `default`, `iv`, etc.). Unknown axis values throw at parse time with a helpful message.

```bash
./bm_doc_list --spec=example_sweep.json --benchmark_out_format=csv --benchmark_out=results.csv
```

Switch the spec's `mode` to `"construct"` to time `construct()` instead of `Search()` — same sweep, same expansion logic.

### Construct-mode `rebuild` (clean timing)

`dret::*::construct()` is idempotent: on the second run with a warm cache it short-circuits in microseconds, which makes the reported time meaningless. Two safeguards in construct mode:

- **`workload.rebuild = true`** in the spec — for each cell, the binary deletes that cell's variant-specific cache files (via a per-cell prefix glob like `<basename>_<bs>-<sf>_gcda_*`) before each timed `construct()` iteration. Shared artefacts (`Text` / `SA` / `DocEnds` / `DA` / `LCP`) are kept, so the (expensive) shared infrastructure isn't rebuilt every iteration. Deletion happens under `state.PauseTiming()` so it doesn't count toward the measurement.

  ```json
  "workload": { "rebuild": true, "min_time": 0.0 }
  ```

  Supported for the **GCDA**, **DGCDA**, **SLP-NS**, **PDL**, and **RMQ** families. Brute baselines ignore `rebuild` — their r-index caches are managed externally.

  **RMQ scope:** each RMQ cell wipes only the per-core RMQ artefacts it builds (`<basename>_sada_rmq_*` / `<basename>_ilcp_rmq_*` + `<basename>_ilcp_run_heads_*` / `<basename>_cilcp_rmq_*` + `<basename>_cilcp_run_heads_*`, plus `<basename>_rmq_n_doc_*`). The SLP / DSLP / DA caches are shared with the corresponding GCDA / DGCDA / brute builds via SDSL type-hashing and stay warm. The reported time is therefore "RMQ-core construction overhead on top of an already-built SLP / DA", which is the right thing to compare across RMQ variants — the SLP build is amortised across the wider sweep.

- **`first_construct_ns` counter** — every construct-mode cell records the wall-clock time of its *first* `construct()` iteration as a Google Benchmark counter, so it shows up in the CSV/JSON output regardless of how many iterations GBenchmark's auto-tuner picked.

- **Stderr warmth warning** — when `rebuild` is off and the first iteration returns in under 1 ms, the binary prints a one-line warning to stderr identifying the cell. Quick visual signal that the measurement is bogus.

### Construct-mode `memory_trace` (sdsl memory log + size sidecar)

Set `"memory_trace": true` in the workload block to make `bm_doc_list` emit, for every construct-mode benchmark cell:

- `construction-<bench>.html` and `construction-<bench>.json` — the sdsl memory-monitor traces (an event-by-event log of allocations during `construct()`).
- `sizes-<bench>.json` — the per-component size breakdown produced by `GetSizeReport()` (also available as inline GBenchmark counters, but the JSON sidecar is convenient for post-processing).

Files land in the working directory keyed by the GBenchmark cell name (e.g. `DocListGCDA-CompactBP-bs1024-sf4`). This is the same artefact set `bm_build_items` emits unconditionally, gated here on the spec flag so a large sweep doesn't litter the working directory. Default is `false`.

## Legacy: `bm_query_doc_list` / `bm_build_items`

These older binaries take gflags directly (one comma-separated list per axis). They predate the JSON spec; kept alongside for one cycle.

```bash
# One-shot query benchmark (builds caches inline if missing):
./bm_query_doc_list \
  --data_dir=/path/to/data --data_name=wiki \
  --patterns=/path/to/queries.txt \
  --gcda_slp_variants=light,compact_bp \
  --min_block_size=256 --max_block_size=1024 \
  --min_storing_factor=2 --max_storing_factor=8 \
  --benchmark_out_format=csv --benchmark_out=results.csv
```

```bash
# If you want construction-time + memory logs as well:
./bm_build_items \
  --data=/path/to/data/wiki \
  --gcda_slp_variants=light,compact_bp \
  --min_block_size=256 --max_block_size=1024 \
  --min_storing_factor=2 --max_storing_factor=8 \
  --benchmark_out_format=json --benchmark_out=build.json
```

### CLI axes (legacy binaries; same values appear in the `bm_doc_list` JSON spec)

| Flag | Values | Meaning |
|---|---|---|
| `--rmq_get_doc_variants` | `da`, `slp`, `slp_ns`, `dslp` | RMQ family's raw-range doc lookup. |
| `--gcda_slp_variants`    | `light`, `compact_bp`, `compact_louds`, `combined` | GCDA's TSLP choice (also reused by RMQ-SLP for cache sharing). `light` = `grammar::LightSLP<...>` (the family default). |
| `--bare_slp_variants`    | `iv`, `raw`, `dv`, `vv` | Non-sampled `grammar::SLP<>` container variants. `iv` = `int_vector` (the family default), `raw` = `std::vector<uint32_t>`, `dv` = `dac_vector`, `vv` = `vlc_vector`. |
| `--dgcda_slp_variants`   | `default`, `otf`, `crl`, `ev`, `dv`, `vv` | DGCDA's TSLP choice (DifferentialLightSLP span-length and inner-container variants). |
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
├── spec.h                 — JSON spec types + parser for bm_doc_list
├── factory.h              — Factory<>: builds + caches index instances; auto-builds missing cache artefacts
├── factories/             — one header per family (brute, gcda, dgcda, slp_ns, rmq, pdl)
├── bm_doc_list.cpp            — JSON-spec-driven benchmark binary (recommended entry point)
├── bm_query_doc_list.cpp  — legacy gflags-driven query binary
├── bm_build_items.cpp     — legacy gflags-driven construction binary (kept for memory-monitor traces)
└── example_sweep.json     — annotated example spec for bm_doc_list
```

To add a new axis value (e.g. a new TSLP variant): one line in `axes.h`, one arm in `enum_traits.h`'s `EnumTraits<>` specialisation, one alias + one Make() arm in the relevant `factories/<family>.h`, and one register-block in each binary that participates in the sweep.
