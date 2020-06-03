#!/usr/bin/env bash

dir_idxs=$1
dir_colls=$2
benchmark=$3

for coll in "$dir_idxs"/*; do
  coll_name=$(basename "$coll")
  echo "Collection $coll_name"

  mkdir -p "$coll_name"
  cd "$coll_name" || exit

  $benchmark --data_dir $dir_idxs/$coll_name --patterns $dir_colls/$coll_name/patterns --benchmark_counters_tabular=true --benchmark_out_format=csv --benchmark_out=$coll_name.csv --print_result

  cd ..
done
