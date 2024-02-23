#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"

#include <vector>


class Run3ScoutingMuonVtxIndexProducer: public edm::stream::EDProducer<>{
  public:
    Run3ScoutingMuonVtxIndexProducer(const edm::ParameterSet &);
    ~Run3ScoutingMuonVtxIndexProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    void endStream() override {}

  private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scoutingmuon_token_;

    std::string num_vertices_name_;
    std::string vertex_index_collection_name_;
};

Run3ScoutingMuonVtxIndexProducer::Run3ScoutingMuonVtxIndexProducer(const edm::ParameterSet &iConfig)
  : input_scoutingmuon_token_(consume(iConfig.getParameter<edm::InputTag>("src"))),
    num_vertices_name_(consume(iConfig.getParameter<std::string>("num_vertices_name"))),
    vertex_index_collection_name_(consume(iConfig.getParameter<std::string>("vertex_index_collection_name"))){
    //produces<Run3ScoutingMuon>();
  produces<edm::ValueMap<uint>>(num_vertices_name_);
  produces<std::vector<int>>(vertex_index_collection_name_);
}

void Run3ScoutingMuonVtxIndexProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup){
  Handle<std::vector<Run3ScoutingMuon>> scouting_muon_handle;
  iEvent.getByToken(input_scoutingmuon_token_, scouting_particle_handle);

  auto num_vertices = std::make_unique<std::vector<uint>>();
  auto vertex_indices = std::make_unique<std::vector<int>>();

  for(auto &muon : *scouting_muon_handle){
    num_vertices -> push_back(muon.vtxIndx.size());
    for(auto &muon_vtx_idx: muon.vtxIndx) vertex_indices -> push_back(muon_vtx_idx);
  }

  edm::OrphanHandle<std::vector<Run3ScoutingMuon>> oh(scouting_muon_handle);

  std::unique_ptr<edm::ValueMap<uint>> num_vectices_VM(new edm::ValueMap<uint>());
  edm::ValueMap<uint>::Filler filler_num_vertices(*num_vectices_VM);
  filler_num_vertices.insert(oh, num_vertices -> begin(), num_vertices -> end());
  filler_num_vertices.fill();
  iEvent.put(std::move(num_vertices_VM), num_vertices_name);

  iEvent.put(std::move(vertex_indices));

  num_vertices_.clear();

}

void Run3ScoutingMuonVtxIndexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions){
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingMuonPacker"));
  desc.add<std::stringl>("num_vertices_name", "ScoutingMuonNumVertices");
  desc.add<std::string>("vertex_index_collection_name", "ScoutingMuonVertexIndexCollection");
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingMuonVtxIndexProducer);