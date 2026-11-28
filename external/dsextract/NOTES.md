# dsextract (vendored)

Cecilia Hernandez's dense-subgraph (biclique) extractor described in
**"Compression of Web and Social graphs via dense subgraphs"** (Hernandez
& Navarro, SPIRE 2012). Used by drl's PDL-BC and consumed by dret's
`pdl::BCCodec` to find shared document subsets that become named rules
in the BC layout.

## Source

Vendored from drl's pre-built copy at
`drl/external/build/dsextract/` (which itself comes from the upstream
tarball at `drl/external/doclist-env/software/dsextract.tgz`). No public
upstream URL is known.

## Patches relative to the upstream tarball

`Shingles.cpp` was modified to drop boost: `boost/tokenizer.hpp` +
`boost/functional/hash.hpp` are unused in the actual code path, and
`boost::hash<string>` was replaced with `std::hash<string>`. Same patch
applied in drl's `external/build/dsextract/`.

The Makefile still has `-I/usr/include/boost` in CFLAGS but no source
file references boost after the patch, so the missing path is harmless.

`vnmextract` takes an optional 8th argument, the seed of the min-hash
functions (`Shingles::seed()`). Upstream calls `srand(time(NULL))` at every
pass (in the `Shingles` constructor and in `cleanAll`), so a run was not
reproducible. With the seed, pass `i` calls `srand(seed + i)` instead: a run
is a function of its input, and each pass still gets other hash functions, as
upstream's passes did, since each takes seconds. Without the seed the
behaviour is upstream's. dret's `BCCodec` always passes one
(`MiningParams::seed`).
