#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <vector>

class Run3ScoutingMuonArrayTableProducer: public edm::stream::EDProducer<>{
  public:
    Run3ScoutingMuonArrayTableProducer(const edm::ParameterSet &);
    ~Run3ScoutingMuonArrayTableProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    void endStream() override {}

  private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scouting_muon_token_;

    std::string vertex_index_collection_name_;
    std::string hit_pattern_collection_name_;
};

Run3ScoutingMuonArrayTableProducer::Run3ScoutingMuonArrayTableProducer(const edm::ParameterSet &iConfig)
  : input_scouting_muon_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    vertex_index_collection_name_(iConfig.getParameter<std::string>("vertex_index_collection_name")),
    hit_pattern_collection_name_(iConfig.getParameter<std::string>("hit_pattern_collection_name")){
  produces<nanoaod::FlatTable>(vertex_index_collection_name_);
  produces<nanoaod::FlatTable>(hit_pattern_collection_name_);
}

void Run3ScoutingMuonArrayTableProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup){
  edm::Handle<std::vector<Run3ScoutingMuon>> scouting_muon_handle;
  iEvent.getByToken(input_scouting_muon_token_, scouting_muon_handle);

  //auto num_vertices = std::make_unique<std::vector<uint>>();
  std::vector<int> vertex_indices;
  std::vector<uint16_t> hit_patterns;

  for(auto &muon : *scouting_muon_handle){
    //num_vertices->push_back(muon.vtxIndx().size());
    for(auto &muon_vtx_idx: muon.vtxIndx()) vertex_indices.push_back(muon_vtx_idx);

    //num_hit_patterns->push_back(muon.trk_hitPattern().hitPattern.size());
    for(auto &muon_hit_pattern: muon.trk_hitPattern().hitPattern) hit_patterns.push_back(muon_hit_pattern);
  }

  auto vertex_index_table = std::make_unique<nanoaod::FlatTable>(vertex_indices.size(), vertex_index_collection_name_, false);
  vertex_index_table->addColumn<int>("vertexIndex", vertex_indices, "flatten list of muon's associated vertex indices");
  iEvent.put(std::move(vertex_index_table), vertex_index_collection_name_);

  auto hit_pattern_table = std::make_unique<nanoaod::FlatTable>(hit_patterns.size(), hit_pattern_collection_name_, false);
  hit_pattern_table->addColumn<int>("hitPattern", hit_patterns, "flatten list of muon's associated hit patterns");
  iEvent.put(std::move(hit_pattern_table), hit_pattern_collection_name_);

}

void Run3ScoutingMuonArrayTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions){
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingMuonPacker"));
  desc.add<std::string>("vertex_index_collection_name", "ScoutingMuonVertexIndex");
  desc.add<std::string>("hit_pattern_collection_name", "ScoutingMuonTrackHitPattern");
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingMuonArrayTableProducer);
