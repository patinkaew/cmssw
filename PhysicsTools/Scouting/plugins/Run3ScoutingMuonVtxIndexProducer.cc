#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

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
    const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scouting_muon_token_;

    std::string num_vertices_name_;
    std::string vertex_index_collection_name_;
};

Run3ScoutingMuonVtxIndexProducer::Run3ScoutingMuonVtxIndexProducer(const edm::ParameterSet &iConfig)
  : input_scouting_muon_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    num_vertices_name_(iConfig.getParameter<std::string>("num_vertices_name")),
    vertex_index_collection_name_(iConfig.getParameter<std::string>("vertex_index_collection_name")){
    //produces<Run3ScoutingMuon>();
  produces<edm::ValueMap<uint>>(num_vertices_name_);
  produces<std::vector<int>>(vertex_index_collection_name_);
  produces<nanoaod::FlatTable>(vertex_index_collection_name_);
}

void Run3ScoutingMuonVtxIndexProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup){
  edm::Handle<std::vector<Run3ScoutingMuon>> scouting_muon_handle;
  iEvent.getByToken(input_scouting_muon_token_, scouting_muon_handle);

  auto num_vertices = std::make_unique<std::vector<uint>>();
  auto vertex_indices = std::make_unique<std::vector<int>>();
  //std::vector<int> vertex_indices;

  for(auto &muon : *scouting_muon_handle){
    num_vertices->push_back(muon.vtxIndx().size());
    std::cout << muon.vtxIndx().size() << std::endl; 
    for(auto &muon_vtx_idx: muon.vtxIndx()) vertex_indices->push_back(muon_vtx_idx);
    for(auto &vertex_index: *vertex_indices) std::cout << vertex_index << std::endl;  
  }

  //edm::OrphanHandle<std::vector<Run3ScoutingMuon>> oh(scouting_muon_handle);

  std::unique_ptr<edm::ValueMap<uint>> num_vertices_VM(new edm::ValueMap<uint>());
  edm::ValueMap<uint>::Filler filler_num_vertices(*num_vertices_VM);
  filler_num_vertices.insert(scouting_muon_handle, num_vertices->begin(), num_vertices->end());
  filler_num_vertices.fill();
  iEvent.put(std::move(num_vertices_VM), num_vertices_name_);

  iEvent.put(std::move(vertex_indices), vertex_index_collection_name_);
  
  auto vertex_index_table = std::make_unique<nanoaod::FlatTable>(vertex_indices->size(), vertex_index_collection_name_, false);
  vertex_index_table->addColumn<int>("VertexIndex", *vertex_indices, "flatten list of muon's associated vertex indices");

  iEvent.put(std::move(vertex_index_table), vertex_index_collection_name_);
}

void Run3ScoutingMuonVtxIndexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions){
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingMuonPacker"));
  desc.add<std::string>("num_vertices_name", "ScoutingMuonNumVertices");
  desc.add<std::string>("vertex_index_collection_name", "ScoutingMuonVertexIndexCollection");
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingMuonVtxIndexProducer);