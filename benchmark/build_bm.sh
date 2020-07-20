#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

dir_colls=$1
repair=$2
dir_grm_tools=${3:-"$SCRIPT_DIR/../build/grammar-build"}
dir_tools=${4:-"$SCRIPT_DIR/../build"}

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

  for raw_file in dsa_raw_data.sdsl doc_disas_raw_data.sdsl da_raw_data.sdsl da_rle_raw_data.sdsl; do
    if [[ ! -f "$raw_file".R || ! -f "$raw_file".C ]]; then
      echo "Re-pair $raw_file"
      "$repair" "$raw_file"
    fi
  done

  # Build GCDA items
  "$dir_tools"/build_items_gcda --data "$coll"/data

  for dslp_file in dsa_raw_data.sdsl doc_disas_raw_data.sdsl; do
    "$dir_grm_tools"/build_dslp_span_sums --data "$dslp_file"

    "$dir_grm_tools"/build_dslp_samples --data "$dslp_file" --max_size 1024
  done

  cd ..
done
