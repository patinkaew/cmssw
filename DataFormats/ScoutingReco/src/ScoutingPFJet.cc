#include "DataFormats/ScoutingReco/interface/ScoutingPFJet.h"

#include "DataFormats/Candidate/interface/Candidate.h"

// constructor
reco::ScoutingPFJet::ScoutingPFJet(const Run3ScoutingPFJet &scoutingPFJet)
  : reco::PFJet(
      reco::Candidate::LorentzVector(reco::Candidate::PolarLorentzVector(scoutingPFJet.pt(), scoutingPFJet.eta(), scoutingPFJet.phi(), scoutingPFJet.m())),
      reco::Candidate::Point(0, 0, 0),
      makeSpecific(scoutingPFJet)
      ) {

  // set jetArea
  setJetArea(scoutingPFJet.jetArea());
}

// makeSpecific
reco::PFJet::Specific reco::ScoutingPFJet::makeSpecific(const Run3ScoutingPFJet &scoutingPFJet) {
    reco::PFJet::Specific specific;
    specific.mChargedHadronEnergy = scoutingPFJet.chargedHadronEnergy();
    specific.mNeutralHadronEnergy = scoutingPFJet.neutralHadronEnergy();
    specific.mPhotonEnergy = scoutingPFJet.photonEnergy();
    specific.mElectronEnergy = scoutingPFJet.electronEnergy();
    specific.mMuonEnergy = scoutingPFJet.muonEnergy();
    specific.mHFHadronEnergy = scoutingPFJet.HFHadronEnergy();
    specific.mHFEMEnergy = scoutingPFJet.HFEMEnergy();

    specific.mChargedHadronMultiplicity = scoutingPFJet.chargedHadronMultiplicity();
    specific.mNeutralHadronMultiplicity = scoutingPFJet.neutralHadronMultiplicity();
    specific.mPhotonMultiplicity = scoutingPFJet.photonMultiplicity();
    specific.mElectronMultiplicity = scoutingPFJet.electronMultiplicity();
    specific.mMuonMultiplicity = scoutingPFJet.muonMultiplicity();
    specific.mHFHadronMultiplicity = scoutingPFJet.HFHadronMultiplicity();
    specific.mHFEMMultiplicity = scoutingPFJet.HFEMMultiplicity();

    specific.mChargedEmEnergy = scoutingPFJet.electronEnergy();
    specific.mChargedMuEnergy = scoutingPFJet.muonEnergy();
    specific.mNeutralEmEnergy = scoutingPFJet.photonEnergy();
    specific.mChargedMultiplicity = scoutingPFJet.chargedHadronMultiplicity() + scoutingPFJet.electronMultiplicity() + scoutingPFJet.muonMultiplicity();
    specific.mNeutralMultiplicity = scoutingPFJet.neutralHadronMultiplicity() + scoutingPFJet.photonMultiplicity() + scoutingPFJet.HFHadronMultiplicity() + scoutingPFJet.HFEMMultiplicity();

    specific.mHOEnergy = scoutingPFJet.HOEnergy();

    return specific;
}

// override numberOfDaugthers
size_t reco::ScoutingPFJet::numberOfDaughters() const {
  return chargedMultiplicity() + neutralMultiplicity();
}
