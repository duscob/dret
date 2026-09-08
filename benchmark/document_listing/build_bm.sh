#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

dir_colls=$1
repair=$2
dir_grm_tools=${3:-"$SCRIPT_DIR/../build/grammar-build"}
dir_tools=${4:-"$SCRIPT_DIR/../build"}
rmq_get_doc_variants=${5:-da}

for coll in "$dir_colls"/*; do
  coll_name=$(basename "$coll")
  echo "Collection $coll_name"

  mkdir -p "$coll_name"
  cd "$coll_name" || exit

  ls -la "$coll"/data
  # Build common items (avoiding following components)
  touch "da_wt_huff_data.sdsl" "dsa_raw_data.sdsl" "doc_isas_data.sdsl" "doc_disas_raw_data.sdsl" "csa_data.sdsl" "bwt_data.sdsl"
  "$dir_tools"/build_items --data "$coll"/data

#  # Build Sada items
#  "$dir_tools"/build_items_sada --data "$coll"/data


#  # Build ILCP items
#  "$dir_tools"/build_items_ilcp --data "$coll"/data


  # Build GCDA items (avoiding following components)
  ln -s da_raw_data.sdsl da
  ln -s da_rle_raw_data.sdsl da_rle
  for raw_file in da da_rle; do
    if [[ ! -f "$raw_file".R || ! -f "$raw_file".C ]]; then
      echo "Re-pair $raw_file"
      "$repair" "$raw_file"
    fi

    "$dir_grm_tools"/bm_build_slp_partition --data="$raw_file" --benchmark_out_format=csv --benchmark_out=build_slp_partition.csv --benchmark_counters_tabular=true

#    "$dir_grm_tools"/bm_build_compact_slp --data="$raw_file" --benchmark_out_format=csv --benchmark_out=build_compact_slp.csv --benchmark_counters_tabular=true
  done

#  touch "gcda_f_occs_data.sdsl" "gcda_l_occs_data.sdsl"
#  "$dir_tools"/build_items_gcda --data "$coll"/data
#
#  for dslp_file in dsa_raw_data.sdsl doc_disas_raw_data.sdsl; do
#    "$dir_grm_tools"/build_dslp_span_sums --data "$dslp_file"
#
#    "$dir_grm_tools"/build_dslp_samples --data "$dslp_file" --max_size 1024
#  done

  if [[ "$rmq_get_doc_variants" != "da" ]]; then
    "$dir_tools"/bm_build_items \
      --data "$coll"/data \
      --rmq_get_doc_variants="$rmq_get_doc_variants" \
      --benchmark_counters_tabular=true \
      --benchmark_out_format=csv \
      --benchmark_out=build_rmq_get_doc.csv
  fi

  cd ..
done
