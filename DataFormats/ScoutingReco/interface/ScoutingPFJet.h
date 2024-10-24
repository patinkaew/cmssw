#ifndef DataFormats_ScoutingReco_ScoutingPFJet_h
#define DataFormats_ScoutingReco_ScoutingPFJet_h

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"

namespace reco {
  class ScoutingPFJet : public reco::PFJet {
  public:
    //default constructor
    ScoutingPFJet() {}
    // constructor from Run3ScoutingPFJet
    ScoutingPFJet(const Run3ScoutingPFJet &scoutingPFJet);

    size_t numberOfDaughters() const override;

  private:
      reco::PFJet::Specific makeSpecific(const Run3ScoutingPFJet &scoutingPFJet);

  };

  typedef std::vector<ScoutingPFJet> ScoutingPFJetCollection;
} // namespace reco
#endif
