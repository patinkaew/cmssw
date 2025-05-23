#ifndef UtilAlgos_StatusSelector_h
#define UtilAlgos_StatusSelector_h
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/Utils/interface/StatusSelector.h"

namespace reco {
  namespace modules {

    template <>
    struct ParameterAdapter<StatusSelector> {
      static StatusSelector make(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC) {
        return StatusSelector(cfg.getParameter<std::vector<int> >("status"));
      }

      static void fillPSetDescription(edm::ParameterSetDescription& desc) { desc.add<std::vector<int> >("status", {}); }
    };

  }  // namespace modules
}  // namespace reco

#endif
