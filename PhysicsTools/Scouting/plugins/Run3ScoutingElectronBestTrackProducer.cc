// -*- C++ -*-
//
// Package:    PhysicsTools/Scouting
// Class:      ScoutingElectronBestTrackProducer
//
/**
 Description: Choose the most suitable track for a given scouting electron
 Implementation:
     Allows for ID selections on the tracks before associating them to the electrons
*/
//
// Original Author:  Abanti Ranadhir Sahasransu
// Adapted by:       Patin Inkaew
//
//

// system include files
#include <limits>
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"

//
// class declaration
//

class ScoutingElectronBestTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit ScoutingElectronBestTrackProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  template <typename T>
  void putValueMap(edm::Event&, edm::Handle<Run3ScoutingElectronCollection>&, const std::vector<T>&, const std::string&);

  const edm::EDGetTokenT<Run3ScoutingElectronCollection> run3ScoutingElectronToken_;
  std::vector<double> trackPtMin_;
  std::vector<double> trackChi2OverNdofMax_;
  std::vector<double> relativeEnergyDifferenceMax_;
  std::vector<double> deltaPhiMax_;
  bool produceBestTrackIndex_;
  bool produceBestTrackVariables_;
};

//
// constructors and destructor
//
ScoutingElectronBestTrackProducer::ScoutingElectronBestTrackProducer(const edm::ParameterSet& iConfig)
    : run3ScoutingElectronToken_(consumes(iConfig.getParameter<edm::InputTag>("Run3ScoutingElectron"))),
      produceBestTrackIndex_(iConfig.getParameter<bool>("produceBestTrackIndex")),
      produceBestTrackVariables_(iConfig.getParameter<bool>("produceBestTrackVariables")) {
 
  trackPtMin_ = iConfig.getParameter<std::vector<double>>("TrackPtMin");
  trackChi2OverNdofMax_ = iConfig.getParameter<std::vector<double>>("TrackChi2OverNdofMax");
  relativeEnergyDifferenceMax_ = iConfig.getParameter<std::vector<double>>("RelativeEnergyDifferenceMax");
  deltaPhiMax_ = iConfig.getParameter<std::vector<double>>("DeltaPhiMax");

  if (trackPtMin_.size() != 2) {
    throw cms::Exception("ScoutingElectronBestTrack")
        << "TrackPtMin must have exactly 2 elements for EB and EE respectively!";
  }
  if (trackChi2OverNdofMax_.size() != 2) {
    throw cms::Exception("ScoutingElectronBestTrack")
        << "TrackChi2OverNdofMax must have exactly 2 elements for EB and EE respectively!";
  }
  if (relativeEnergyDifferenceMax_.size() != 2) {
    throw cms::Exception("ScoutingElectronBestTrack")
        << "RelativeEnergyDifferenceMax must have exactly 2 elements for EB and EE respectively!";
  }
  if (deltaPhiMax_.size() != 2) {
    throw cms::Exception("ScoutingElectronBestTrack")
        << "DeltaPhiMax must have exactly 2 elements for EB and EE respectively!";
  }

  if (produceBestTrackIndex_) {
    produces<edm::ValueMap<int>>("bestTrack-index");
  }

  if (produceBestTrackVariables_) {
    produces<edm::ValueMap<float>>("bestTrack-d0");
    produces<edm::ValueMap<float>>("bestTrack-dz");
    produces<edm::ValueMap<float>>("bestTrack-pt");
    produces<edm::ValueMap<float>>("bestTrack-eta");
    produces<edm::ValueMap<float>>("bestTrack-phi");
    produces<edm::ValueMap<float>>("bestTrack-pMode");
    produces<edm::ValueMap<float>>("bestTrack-etaMode");
    produces<edm::ValueMap<float>>("bestTrack-phiMode");
    produces<edm::ValueMap<float>>("bestTrack-qoverpModeError");
    produces<edm::ValueMap<float>>("bestTrack-chi2overndf");
    produces<edm::ValueMap<int>>("bestTrack-charge");
  }
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void ScoutingElectronBestTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<Run3ScoutingElectronCollection> run3ScoutingElectronHandle;
  iEvent.getByToken(run3ScoutingElectronToken_, run3ScoutingElectronHandle);

  if (!run3ScoutingElectronHandle.isValid()) {
    // Handle the absence as a warning
    edm::LogWarning("ScoutingElectronBestTrack") << "No Run3ScoutingElectron collection found in the event!";
    return;
  }
  
  size_t numElectrons = run3ScoutingElectronHandle->size();
  std::vector<int> bestTrack_indices(numElectrons, -1);

  for (size_t iElectron = 0; iElectron < numElectrons; iElectron++) {
    const auto& electron = run3ScoutingElectronHandle->at(iElectron);
    const math::PtEtaPhiMLorentzVector cluster(electron.pt(), electron.eta(), electron.phi(), 0.0005);

    int bestTrack_index = -1;
    double bestTrack_ediff = std::numeric_limits<double>::max();
    
    // find bestTrack index
    for (unsigned int i = 0; i < electron.trkpt().size(); ++i) {
      const unsigned int eta_idx = (abs(electron.trketa().at(i)) < 1.479) ? 0 : 1; // EB=0, EE=1
      if (electron.trkpt().at(i) < trackPtMin_.at(eta_idx))
        continue;
      if (electron.trkchi2overndf().at(i) > trackChi2OverNdofMax_.at(eta_idx))
        continue;

      const math::PtEtaPhiMLorentzVector gsftrack(
          electron.trkpt().at(i), electron.trketa().at(i), electron.trkphi().at(i), 0.0005);

      if (deltaPhi(cluster.phi(), gsftrack.phi()) > deltaPhiMax_.at(eta_idx))
        continue;

      const double track_ediff = abs((cluster.energy() - gsftrack.energy()) / cluster.energy());
      if (track_ediff > relativeEnergyDifferenceMax_.at(eta_idx))
        continue;

      if (track_ediff < bestTrack_ediff) {
        bestTrack_ediff = track_ediff;
        bestTrack_index = i;
      }
    }
    bestTrack_indices[iElectron] = bestTrack_index;
  }

  if (produceBestTrackIndex_) {
    putValueMap<int>(iEvent, run3ScoutingElectronHandle, bestTrack_indices, "bestTrack-index");
  }

  if (produceBestTrackVariables_) {
    // translate bestTrack's index to bestTrack's variables
    std::vector<float> bestTrack_d0s(numElectrons, -99);
    std::vector<float> bestTrack_dzs(numElectrons, -99);
    std::vector<float> bestTrack_pts(numElectrons, -99);
    std::vector<float> bestTrack_etas(numElectrons, -99);
    std::vector<float> bestTrack_phis(numElectrons, -99);
    std::vector<float> bestTrack_pModes(numElectrons, -99);
    std::vector<float> bestTrack_etaModes(numElectrons, -99);
    std::vector<float> bestTrack_phiModes(numElectrons, -99);
    std::vector<float> bestTrack_qoverpModeErrors(numElectrons, -99);
    std::vector<float> bestTrack_chi2overndfs(numElectrons, -99);
    std::vector<int> bestTrack_charges(numElectrons, -99);

    for (size_t iElectron = 0; iElectron < numElectrons; iElectron++) {
      int bestTrack_index = bestTrack_indices[iElectron];
      if (bestTrack_index > -1) {
        const auto& electron = run3ScoutingElectronHandle->at(iElectron);
        bestTrack_d0s[iElectron] = electron.trkd0()[bestTrack_index];
        bestTrack_dzs[iElectron] = electron.trkdz()[bestTrack_index];
        bestTrack_pts[iElectron] = electron.trkpt()[bestTrack_index];
        bestTrack_etas[iElectron] = electron.trketa()[bestTrack_index];
        bestTrack_phis[iElectron] = electron.trkphi()[bestTrack_index];
        bestTrack_pModes[iElectron] = electron.trkpMode()[bestTrack_index];
        bestTrack_etaModes[iElectron] = electron.trketaMode()[bestTrack_index];
        bestTrack_phiModes[iElectron] = electron.trkphiMode()[bestTrack_index];
        bestTrack_qoverpModeErrors[iElectron] = electron.trkqoverpModeError()[bestTrack_index];
        bestTrack_chi2overndfs[iElectron] = electron.trkchi2overndf()[bestTrack_index];
        bestTrack_charges[iElectron] = electron.trkcharge()[bestTrack_index];
      }
    }

    // put into iEvent
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_d0s, "bestTrack-d0");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_dzs, "bestTrack-dz");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_pts, "bestTrack-pt");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_etas, "bestTrack-eta");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_phis, "bestTrack-phi");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_pModes, "bestTrack-pMode");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_etaModes, "bestTrack-etaMode");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_phiModes, "bestTrack-phiMode");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_qoverpModeErrors, "bestTrack-qoverpModeError");
    putValueMap<float>(iEvent, run3ScoutingElectronHandle, bestTrack_chi2overndfs, "bestTrack-chi2overndf");
    putValueMap<int>(iEvent, run3ScoutingElectronHandle, bestTrack_charges, "bestTrack-charge");
  }

}

template <typename T>
void ScoutingElectronBestTrackProducer::putValueMap(edm::Event& iEvent, edm::Handle<Run3ScoutingElectronCollection>& handle, const std::vector<T>& values, const std::string& label) {
  std::unique_ptr<edm::ValueMap<T>> valuemap(new edm::ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valuemap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valuemap), label);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingElectronBestTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>(("Run3ScoutingElectron"), edm::InputTag("hltScoutingEgammaPacker"));
  desc.add<std::vector<double>>(("TrackPtMin"), {0.0, 0.0});
  desc.add<std::vector<double>>(("TrackChi2OverNdofMax"), {9999.0, 9999.0});
  desc.add<std::vector<double>>(("RelativeEnergyDifferenceMax"), {9999.0, 9999.0});
  desc.add<std::vector<double>>(("DeltaPhiMax"), {9999.0, 9999.0});
  desc.add<bool>("produceBestTrackIndex", false);
  desc.add<bool>("produceBestTrackVariables", true);
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingElectronBestTrackProducer);
