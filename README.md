DRet - Document Retrieval Library
=================================

DRet is a C++ library containing Document RETrieval indexes for repetitive scenarios. Currently, it contains solutions for the following problem: 

>**Document listing with frequencies**: Given a document collection and a pattern, return the set of documents where the pattern occurs and their frequencies.

The implemented indexes are described in:

>Cobas D, Mäkinen V, Rossi M. Tailoring r-index for document listing towards metagenomics applications. To appear in Proceedings of the SPIRE 2020: 27th International Symposium on String Processing and Information Retrieval.

Getting Started
-----

### Downloading
To download our project, execute:
```shell
$ git clone https://github.com/duscob/dret.git
```

### Compiling
DRet is a C++ template library using standard C++14.
As a build system, the project uses [CMake](https://cmake.org/).

>##### Dependencies
>
>The DRet library requires the following external libraries:
>  * [Succinct Data Structure Library (SDSL)](https://github.com/simongog/sdsl-lite "SDSL's GitHub repository")
>  * [Grammar Library](https://github.com/duscob/grammar "Grammar's GitHub repository")
>
>`SDSL` must be installed on your system. In the case of `grammar` library, our project downloads it from GitHub.


##### Building

To build our solution, you can create a build folder on the source directory and move to it.
```shell
$ cd dret
$ mkdir build
$ cd build
```

Then build and compile the project using the commands `cmake` and `make`, respectively.
```shell
$ cmake ..
$ make
```


### Benchmarks

`bm_doc_list` builds and queries every index from a JSON spec: the collection, the mode (`construct`
or `query`) and one block per index family. In query mode, `workload.task` selects what a query
returns: the documents of the pattern (`listing`, the default) or each document with its frequency
(`frequency`). Every index computes listing; frequencies come from the brute-force indexes, GCDA with
`counts`, and the `rmq_freq` (SADA / ILCP / CILCP) and `gcda_occs` families, the last two with a
frequency oracle (`rank` over a wavelet tree on the document array, or `doc_isa`, the per-document
inverse suffix arrays over a differential suffix array). See `benchmark/document_listing/README.md`
for the spec, and `test/regression/` for a sweep of every index in each task.

>##### Dependencies
>
>The grammar-compressed indexes run the balanced Re-Pair for integers by Gonzalo Navarro
>([here](https://users.dcc.uchile.cl/~gnavarro/software/repair.tgz)), built with the project.
>The benchmarks require the [GFlags](https://gflags.github.io/gflags/) library.

The 2020 tools and scripts of the paper above (`build_items*`, `build_bm.sh`, `run_bm.sh` and
`bm_query_doc_freq_index`) were replaced by these indexes; they are in the history up to commit
`52b4525`.
