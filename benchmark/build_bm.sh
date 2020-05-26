#!/usr/bin/env bash

dir_colls=$1
dir_tools=$2
repair=$3
dir_grm_tools=$4

for coll in "$dir_colls"/*; do
  coll_name=$(basename "$coll")
  echo "Collection $coll_name"

  mkdir -p "$coll_name"
  cd "$coll_name" || exit

  ls -la "$coll"/data
  # Build common items
  "$dir_tools"/build_items --data "$coll"/data

  # Build Sada items
  "$dir_tools"/build_items_sada --data "$coll"/data

  # Build ILCP items
  "$dir_tools"/build_items_ilcp --data "$coll"/data

  for raw_file in *raw*.sdsl; do
    if [[ ! -f "$raw_file".R || ! -f "$raw_file".C ]]; then
      echo "Re-pair $raw_file"
      "$repair" "$raw_file"
    fi
  done

  for dslp_file in dsa_raw_data.sdsl doc_disas_raw_data.sdsl; do
    "$dir_grm_tools"/build_dslp_span_sums --data "$dslp_file"

    "$dir_grm_tools"/build_dslp_samples --data "$dslp_file" --max_size 1024
  done

  #  basename=(`cat $dir_colls/$coll/sample_rates.txt`)
  #  echo $basename

  #  $benchmark --data $dir_colls/$coll/$basename --patterns $dir_colls/$coll/patterns --benchmark_counters_tabular=true --benchmark_out_format=csv --benchmark_out=$coll.csv

  cd ..
done
