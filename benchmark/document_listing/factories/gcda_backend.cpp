//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode factory of GCDA over a document-array backend, compiled on its own (see query.h).
//

#include "gcda.h"
#include "query.h"

namespace bench::factories::query {

Built GCDABackend(Storage t_storage,
                  dret::Config& t_config,
                  uint32_t t_block_size,
                  float t_storing_factor,
                  axes::GetDocEnum t_get_doc,
                  axes::GCDASetsCodec t_sets) {
  return gcda::MakeBackend(t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_sets);
}

}  // namespace bench::factories::query
