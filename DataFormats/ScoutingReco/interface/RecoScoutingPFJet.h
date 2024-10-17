#ifndef DataFormats_RecoScoutingPFJet_h
#define DataFormats_RecoScoutingPFJet_h

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"

namespace reco {
  class ScoutingPFJet : public reco::PFJet {
  public:
    ScoutingPFJet(const Run3ScoutingPFJet &scoutingPFJet);

  private:
    makeSpecific(const Run3ScoutingPFJet &scoutingPFJet);

  };
} // namespace reco

#endif
