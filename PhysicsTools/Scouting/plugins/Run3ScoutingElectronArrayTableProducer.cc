#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <vector>
#include <cassert>

class Run3ScoutingElectronArrayTableProducer: public edm::stream::EDProducer<>{
  public:
    Run3ScoutingElectronArrayTableProducer(const edm::ParameterSet &);
    ~Run3ScoutingElectronArrayTableProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    //void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    //void endStream() override {}

  private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> input_scouting_electron_token_;

    //std::string num_tracks_name_;
    std::string track_collection_name_;
};

Run3ScoutingElectronArrayTableProducer::Run3ScoutingElectronArrayTableProducer(const edm::ParameterSet &iConfig)
  : input_scouting_electron_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    //num_tracks_name_(iConfig.getParameter<std::string>("num_tracks_name")),
    track_collection_name_(iConfig.getParameter<std::string>("track_collection_name")){
  //produces<edm::ValueMap<uint>>(num_tracks_name_);
  produces<nanoaod::FlatTable>(track_collection_name_);
}

void Run3ScoutingElectronArrayTableProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup){
  edm::Handle<std::vector<Run3ScoutingElectron>> scouting_electron_handle;
  iEvent.getByToken(input_scouting_electron_token_, scouting_electron_handle);

  //auto num_tracks = std::make_unique<std::vector<uint>>();
  std::vector<float> trkd0, trkdz, trkpt, trketa, trkphi, trkchi2overndf;
  //std::vector<float> trkpMode, trketaMode, trkphiMode, trkqoverpModeError;
  std::vector<int> trkcharge;

  for(auto &electron : *scouting_electron_handle){
    size_t n_trks = electron.trkd0().size();
    assert((void("number of tracks mismatch"), 
      (electron.trkdz().size() == n_trks) 
      && (electron.trkpt().size() == n_trks)
      && (electron.trketa().size() == n_trks) 
      && (electron.trkphi().size() == n_trks)
      // && (electron.trkpMode().size() == n_trks) 
      // && (electron.trketaMode().size() == n_trks)
      // && (electron.trkphiMode().size() == n_trks) 
      // && (electron.trkqoverpModeError().size() == n_trks)
      && (electron.trkchi2overndf().size() == n_trks)
      && (electron.trkcharge().size() == n_trks)
    ));

    //num_tracks->push_back(n_trks);

    for(size_t itrk=0; itrk < n_trks; itrk++){
      trkd0.push_back(electron.trkd0()[itrk]);
      trkdz.push_back(electron.trkdz()[itrk]);
      trkpt.push_back(electron.trkpt()[itrk]);
      trketa.push_back(electron.trketa()[itrk]);
      trkphi.push_back(electron.trkphi()[itrk]);
      // trkpMode.push_back(electron.trkpMode()[itrk]);
      // trketaMode.push_back(electron.trketaMode()[itrk]);
      // trkphiMode.push_back(electron.trkphiMode()[itrk]);
      // trkqoverpModeError.push_back(electron.trkqoverpModeError()[itrk]);
      trkchi2overndf.push_back(electron.trkchi2overndf()[itrk]);
      trkcharge.push_back(electron.trkcharge()[itrk]);
    }
  }

  // std::unique_ptr<edm::ValueMap<uint>> num_tracks_VM(new edm::ValueMap<uint>());
  // edm::ValueMap<uint>::Filler filler_num_tracks(*num_tracks_VM);
  // filler_num_tracks.insert(scouting_electron_handle, num_tracks->begin(), num_tracks->end());
  // filler_num_tracks.fill();
  // iEvent.put(std::move(num_tracks_VM), num_tracks_name_);

  auto track_table = std::make_unique<nanoaod::FlatTable>(trkd0.size(), track_collection_name_, false);
  track_table->addColumn<float>("trkd0", trkd0, "flatten list of electron's associated track d0");
  track_table->addColumn<float>("trkdz", trkd0, "flatten list of electron's associated track dz");
  track_table->addColumn<float>("trkpt", trkd0, "flatten list of electron's associated track pt");
  track_table->addColumn<float>("trketa", trkd0, "flatten list of electron's associated track eta");
  track_table->addColumn<float>("trkphi", trkd0, "flatten list of electron's associated track phi");
  // track_table->addColumn<float>("trkpMode", trkd0, "flatten list of electron's associated track pMode");
  // track_table->addColumn<float>("trketaMode", trkd0, "flatten list of electron's associated track etaMode");
  // track_table->addColumn<float>("trkphiMode", trkd0, "flatten list of electron's associated track phiMode");
  // track_table->addColumn<float>("trkdqoverpModeError", trkd0, "flatten list of electron's associated track qoverpModeError");
  track_table->addColumn<float>("trkchi2overndf", trkd0, "flatten list of electron's associated track chi2overndf");
  track_table->addColumn<int>("trkcharge", trkd0, "flatten list of electron's associated track charge");

  iEvent.put(std::move(track_table), track_collection_name_);
}

void Run3ScoutingElectronArrayTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions){
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingEgammaPacker"));
  //desc.add<std::string>("num_tracks_name", "ScoutingElectronNumTracks");
  desc.add<std::string>("track_collection_name", "ScoutingElectronTrack");
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingElectronArrayTableProducer);