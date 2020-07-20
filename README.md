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


### Tools and Benchmarks

We include a set of tools to preprocess the datasets and build the document retrieval indexes.
To facilitate the building process, we also provide some shell scripts to orchestrate the flow execution.

>##### Dependencies
>
>Notice that some of the implemented solutions are grammar-compression based. In our experimental setup, we used the balanced re-pair implementation for integers provided by Gonzalo Navarro [here](https://users.dcc.uchile.cl/~gnavarro/software/repair.tgz).
>
>Our tools and benchmarks also require that the [GFlags](https://gflags.github.io/gflags/) library is installed on your system.


##### Building the Indexes

To build the indexes, you can execute the following script

```shell
$ <source_dir>/benchmark/build_bm.sh <datasets_dir> <source_dir>/build <irepair_exe> <grammar_source_dir>/build 
```

where
 * *source_dir* is the source directory of our project.
 * *datasets_dir* is the parent directory of the datasets. Each subdirectory is a different dataset containing two files: `data` and `patterns`. The file `data` is a binary file that represents the collection, where each document is separated by the symbol **0** (zero). The file must not contain symbols **1** and **2**.
 * *irepair_exe* is the absolute path of the balanced re-pair executable for integer sequences.

All paths must be absolute.

This script creates a folder for each dataset in the current directory and stores on it the data structures required by each index.


##### Benchmarking the Indexes
To run the benchmarks of the indexes, you can execute the following script.

```shell
$ <source_dir>/benchmark/run_bm.sh <indexes_dir> <datasets_dir> 
```

where
 * *indexes_dir* is the directory that stores the indexes built for the collection, i.e., the current directory where you executed the previous command.
 * *datasets_dir* is the parent directory of the datasets. Each subdirectory is a different dataset containing two files: `data` and `patterns`. The file `patterns` is a text file that contains the query patterns, one per line.

All paths must be absolute.

This script creates a folder for each dataset in the current directory, and stores on it the query results by each index and the required space and time.

