# Document-listing benchmarks

One binary, `bm_doc_list`, benchmarks every document-listing index: query throughput or construction
time (selected per spec), with the size breakdown each index reports through `GetSizeReport()`. A
JSON spec (`-spec=path.json`) declares the whole sweep; missing cache artefacts are built on the first
run. The axis enums (`bench::axes::*`) and their parsers (`bench::axes::EnumTraits<E>`) define the
names the spec accepts.

Compiled with `-DDRET_DOC_LIST_PROFILE` (the `cmake-build-profile` tree), query mode also reports the
phase profile of the sampled-tree indexes (GCDA, GCDA over a backend, PDL) as counters per query:
`prof_count_ns_x_q`, `prof_cover_ns_x_q`, `prof_expand_ns_x_q`, `prof_combine_ns_x_q`,
`prof_nodes_x_q`, `prof_rawpos_x_q`, `prof_docs_x_q` (see `dret/doc_list/search_profile.h`).

## `bm_doc_list` — JSON-driven sweep

See `example_sweep.json` for a complete reference. The spec is a single JSON document with four top-level keys:

```json
{
  "mode":     "query",           // "query" | "construct"
  "dataset":  { "dir": "...", "name": "...", "patterns": "...", "data_width": 8 },
  "workload": { "reps": 10, "min_time": 0.0, "task": "listing" },   // task: "listing" | "frequency"
  "sweep":    [ ...family blocks... ]
}
```

Each entry in `sweep` is one family block. The driver expands each block's Cartesian product against the family's axes and registers one Google Benchmark per cell:

| Family | Axes |
|---|---|
| `gcda`   | `tslp` (`light` / `compact_bp` / `compact_louds` / `combined`), `sets` (`rp` / `plain`), `counts` (`none` / `plain`), `block_size`, `storing_factor` |
| `dgcda`  | `tslp` (`default` / `otf` / `crl` / `ev` / `dv` / `vv`), `block_size`, `storing_factor` |
| `slp_ns` | `tslp` (`iv` / `raw` / `dv` / `vv`) |
| `rmq`    | `core` (`sada` / `ilcp` / `cilcp`), `get_doc` (`da` / `slp` / `slp_ns` / `dslp`), `gcda_slp`, `bare_slp`, `block_size`, `storing_factor` |
| `pdl`    | `codec` (`plain` / `rp` / `bc`), `get_doc` (`da` / `slp` / `dslp`), `policy`, `block_size`, `storing_factor` |
| `brute`  | `kind` (`r-index` / `sr-index`), `sampling_size` |
| `rmq_freq`  | `core` (`sada` / `ilcp` / `cilcp`), `oracle` (`rank` / `doc_isa`), `block_size` |
| `gcda_occs` | `oracle` (`rank` / `doc_isa`), `block_size` |

With `"task": "frequency"` a query returns each document with its frequency, and the answers dump
(`-print_results_dir`) writes `d:f` per document. Only the indexes that compute frequencies can run
it (the others are skipped with an error): `brute` (all but `count`), `slp_ns`, `gcda` with
`counts: plain`, `rmq_freq` and `gcda_occs`. The last two need a frequency oracle: `rank`, over a
wavelet tree on the document array, which is then also their backend; or `doc_isa`, the
per-document inverse suffix arrays over a differential suffix array (their backend), whose grammars
sample every `block_size` positions.

Missing axes default to the family's canonical value (`light`, `default`, `iv`, etc.). Unknown axis values throw at parse time with a helpful message.

```bash
./bm_doc_list --spec=example_sweep.json --benchmark_out_format=csv --benchmark_out=results.csv
```

Switch the spec's `mode` to `"construct"` to time `construct()` instead of `Search()`. Both modes register the same cells under the same names (`FAMILY[k=v|...]`): a family turns its sweep block into cells once (`factories/cells.h`), and construct mode builds every index the query mode would load.

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

Files land in the working directory keyed by the benchmark name, sanitized for a file name (every character but `[A-Za-z0-9.-]` becomes `_`). Gated on the spec flag so a large sweep doesn't litter the working directory. Default is `false`.

## Output formats (Google Benchmark native)

`benchmark::Initialize` picks these up automatically:

- **Console** (default): the usual pretty-printed table to stdout.
- **CSV**: `--benchmark_out_format=csv --benchmark_out=path.csv`
- **JSON**: `--benchmark_out_format=json --benchmark_out=path.json`

Each output row carries the per-component size breakdown emitted by `GetSizeReport()` as named counters (e.g. `slp_basic_slp`, `slp_sets_objs`, `count_idx`), plus the totals `Size(bytes)`, `Bits_x_Symbol`, `Time_x_Pattern`, etc.


## Where the code lives

```
benchmark/document_retrieval/
├── axes.h                 — the axis enums at namespace bench::axes scope
├── enum_traits.h          — EnumTraits<E>::ParseCSV / ::Name, plus ParseCSV<E>(csv) helper
├── spec.h                 — JSON spec types + parser for bm_doc_list
├── factory.h              — Factory<>: query mode's shared storage; builds (if missing) and loads each cell's index
├── factories/
│   ├── cells.h            — Cell (name, how to make the index, rebuild prefixes) and Cells(sweep block) per family
│   ├── names.h            — the FAMILY[k=v|...] name helpers
│   ├── <family>.h         — the family's types (brute, gcda, dgcda, slp_ns, rmq, pdl)
│   └── <family>.cpp       — its cells: enumeration, names, types (RMQ: rmq.cpp + rmq_<core>.cpp; PDL: pdl.cpp + pdl_<codec>.cpp)
├── bm_doc_list.cpp        — the driver: spec, one query kernel, one construct kernel, registration of the cells
└── example_sweep.json     — annotated example spec for bm_doc_list
```

To add a new axis value (e.g. a new TSLP variant): one line in `axes.h`, one arm in `enum_traits.h`'s `EnumTraits<>` specialisation, an alias in the family's `factories/<family>.h`, and its arm and name keys in the family's `.cpp`; both modes pick it up. A new family adds its sweep block to `spec.h`, a `Cells()` overload in `factories/cells.h`, its source file, and that file in `bm_doc_list_family_sources` (`benchmark/CMakeLists.txt`).

Each family is compiled in its own source file so that a change to one family cannot move another's query times through GCC's inlining choices (in one large translation unit it did, by a few percent; see `factories/cells.h`).

