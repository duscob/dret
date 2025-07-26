//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include "doc_list_index.h"

namespace dret {

template <typename TStorage,
          typename TComputeSARange,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
class DLSampledTreeScheme : public DocListIndex {};

//~~~~~~~


template <typename TStorage,
          typename TComputeSARange,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
void constructItems(
    DLSampledTreeScheme<TStorage, TComputeSARange, TComputeCover, TGetDocs, TGetDocSet, TMergeSets>& t_index,
    const std::string& t_data_path,
    Config& t_config) {}

//~~~~~~~

}  // namespace dret
