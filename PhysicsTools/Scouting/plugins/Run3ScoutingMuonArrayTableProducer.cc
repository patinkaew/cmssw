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
    //void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    //void endStream() override {}

  private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scouting_muon_token_;

    //std::string num_vertices_name_;
    std::string vertex_index_collection_name_;
    //std::string num_hit_patterns_name_;
    std::string hit_pattern_collection_name_;
};

Run3ScoutingMuonArrayTableProducer::Run3ScoutingMuonArrayTableProducer(const edm::ParameterSet &iConfig)
  : input_scouting_muon_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    //num_vertices_name_(iConfig.getParameter<std::string>("num_vertices_name")),
    vertex_index_collection_name_(iConfig.getParameter<std::string>("vertex_index_collection_name")),
    //num_hit_patterns_name_(iConfig.getParameter<std::string>("num_hit_patterns_name")),
    hit_pattern_collection_name_(iConfig.getParameter<std::string>("hit_pattern_collection_name")){
  //produces<edm::ValueMap<uint>>(num_vertices_name_);
  produces<nanoaod::FlatTable>(vertex_index_collection_name_);
  /*produces<edm::ValueMap<uint8_t>>("trk_hitPattern_hitCount");
  produces<edm::ValueMap<uint8_t>>("trk_hitPattern_beginTrackHits");
  produces<edm::ValueMap<uint8_t>>("trk_hitPattern_endTrackHits");
  produces<edm::ValueMap<uint8_t>>("trk_hitPattern_beginInner");
  produces<edm::ValueMap<uint8_t>>("trk_hitPattern_endInner");
  produces<edm::ValueMap<uint8_t>>("trk_hitPattern_beginOuter");
  produces<edm::ValueMap<uint8_t>>("trk_hitPattern_endOuter");*/
  //produces<edm::ValueMap<uint>>(num_hit_patterns_name_);
  produces<nanoaod::FlatTable>(hit_pattern_collection_name_);
}

void Run3ScoutingMuonArrayTableProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup){
  edm::Handle<std::vector<Run3ScoutingMuon>> scouting_muon_handle;
  iEvent.getByToken(input_scouting_muon_token_, scouting_muon_handle);

  //auto num_vertices = std::make_unique<std::vector<uint>>();
  std::vector<int> vertex_indices;
  /*
  auto trk_hitPattern_hitCount = std::make_unique<std::vector<uint8_t>>();
  auto trk_hitPattern_beginTrackHits = std::make_unique<std::vector<uint8_t>>();
  auto trk_hitPattern_endTrackHits = std::make_unique<std::vector<uint8_t>>();
  auto trk_hitPattern_beginInner = std::make_unique<std::vector<uint8_t>>();
  auto trk_hitPattern_endInner = std::make_unique<std::vector<uint8_t>>();
  auto trk_hitPattern_beginOuter = std::make_unique<std::vector<uint8_t>>();
  auto trk_hitPattern_endOuter = std::make_unique<std::vector<uint8_t>>();
  */
  //auto num_hit_patterns = std::make_unique<std::vector<uint>>();
  std::vector<uint16_t> hit_patterns;

  for(auto &muon : *scouting_muon_handle){
    //num_vertices->push_back(muon.vtxIndx().size());
    for(auto &muon_vtx_idx: muon.vtxIndx()) vertex_indices.push_back(muon_vtx_idx);
    /*
    trk_hitPattern_hitCount->push_back(muon.trk_hitPattern().hitCount);
    trk_hitPattern_beginTrackHits->push_back(muon.trk_hitPattern().beginTrackHits);
    trk_hitPattern_endTrackHits->push_back(muon.trk_hitPattern().endTrackHits);
    trk_hitPattern_beginInner->push_back(muon.trk_hitPattern().beginInner);
    trk_hitPattern_endInner->push_back(muon.trk_hitPattern().endInner);
    trk_hitPattern_beginOuter->push_back(muon.trk_hitPattern().beginOuter);
    trk_hitPattern_endOuter->push_back(muon.trk_hitPattern().endOuter);
    */

    //num_hit_patterns->push_back(muon.trk_hitPattern().hitPattern.size());
    for(auto &muon_hit_pattern: muon.trk_hitPattern().hitPattern) hit_patterns.push_back(muon_hit_pattern);
  }

  // vtx Index
  // std::unique_ptr<edm::ValueMap<uint>> num_vertices_VM(new edm::ValueMap<uint>());
  // edm::ValueMap<uint>::Filler filler_num_vertices(*num_vertices_VM);
  // filler_num_vertices.insert(scouting_muon_handle, num_vertices->begin(), num_vertices->end());
  // filler_num_vertices.fill();
  // iEvent.put(std::move(num_vertices_VM), num_vertices_name_);

  auto vertex_index_table = std::make_unique<nanoaod::FlatTable>(vertex_indices.size(), vertex_index_collection_name_, false);
  vertex_index_table->addColumn<int>("vertexIndex", vertex_indices, "flatten list of muon's associated vertex indices");
  iEvent.put(std::move(vertex_index_table), vertex_index_collection_name_);
  /*
  // hit Pattern (primitive)
  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_hitCount_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_hitCount(*trk_hitPattern_hitCount_VM);
  filler_trk_hitPattern_hitCount.insert(scouting_muon_handle, trk_hitPattern_hitCount->begin(), trk_hitPattern_hitCount->end());
  filler_trk_hitPattern_hitCount.fill();
  iEvent.put(std::move(trk_hitPattern_hitCount_VM), "trk_hitPattern_hitCount");

  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_beginTrackHits_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_beginTrackHits(*trk_hitPattern_beginTrackHits_VM);
  filler_trk_hitPattern_beginTrackHits.insert(scouting_muon_handle, trk_hitPattern_beginTrackHits->begin(), trk_hitPattern_beginTrackHits->end());
  filler_trk_hitPattern_beginTrackHits.fill();
  iEvent.put(std::move(trk_hitPattern_beginTrackHits_VM), "trk_hitPattern_beginTrackHits");

  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_endTrackHits_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_endTrackHits(*trk_hitPattern_beginTrackHits_VM);
  filler_trk_hitPattern_endTrackHits.insert(scouting_muon_handle, trk_hitPattern_endTrackHits->begin(), trk_hitPattern_endTrackHits->end());
  filler_trk_hitPattern_endTrackHits.fill();
  iEvent.put(std::move(trk_hitPattern_endTrackHits_VM), "trk_hitPattern_endTrackHits");

  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_beginInner_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_beginInner(*trk_hitPattern_beginInner_VM);
  filler_trk_hitPattern_beginInner.insert(scouting_muon_handle, trk_hitPattern_beginInner->begin(), trk_hitPattern_beginInner->end());
  filler_trk_hitPattern_beginInner.fill();
  iEvent.put(std::move(trk_hitPattern_beginInner_VM), "trk_hitPattern_beginInner");

  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_endInner_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_endInner(*trk_hitPattern_endInner_VM);
  filler_trk_hitPattern_endInner.insert(scouting_muon_handle, trk_hitPattern_endInner->begin(), trk_hitPattern_endInner->end());
  filler_trk_hitPattern_endInner.fill();
  iEvent.put(std::move(trk_hitPattern_endInner_VM), "trk_hitPattern_endInner");

  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_beginOuter_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_beginOuter(*trk_hitPattern_beginOuter_VM);
  filler_trk_hitPattern_beginOuter.insert(scouting_muon_handle, trk_hitPattern_beginOuter->begin(), trk_hitPattern_beginOuter->end());
  filler_trk_hitPattern_beginOuter.fill();
  iEvent.put(std::move(trk_hitPattern_beginOuter_VM), "trk_hitPattern_beginOuter");

  std::unique_ptr<edm::ValueMap<uint8_t>> trk_hitPattern_endOuter_VM(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_trk_hitPattern_endOuter(*trk_hitPattern_endOuter_VM);
  filler_trk_hitPattern_endOuter.insert(scouting_muon_handle, trk_hitPattern_endOuter->begin(), trk_hitPattern_endOuter->end());
  filler_trk_hitPattern_endOuter.fill();
  iEvent.put(std::move(trk_hitPattern_endOuter_VM), "trk_hitPattern_endOuter");
  */

  // hit pattern (vector)
  // std::unique_ptr<edm::ValueMap<uint>> num_hit_patterns_VM(new edm::ValueMap<uint>());
  // edm::ValueMap<uint>::Filler filler_num_hit_patterns(*num_hit_patterns_VM);
  // filler_num_hit_patterns.insert(scouting_muon_handle, num_hit_patterns->begin(), num_hit_patterns->end());
  // filler_num_hit_patterns.fill();
  // iEvent.put(std::move(num_hit_patterns_VM), num_hit_patterns_name_);

  auto hit_pattern_table = std::make_unique<nanoaod::FlatTable>(hit_patterns.size(), hit_pattern_collection_name_, false);
  hit_pattern_table->addColumn<int>("hitPattern", hit_patterns, "flatten list of muon's associated hit patterns");
  iEvent.put(std::move(hit_pattern_table), hit_pattern_collection_name_);

}

void Run3ScoutingMuonArrayTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions){
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingMuonPacker"));
  //desc.add<std::string>("num_vertices_name", "ScoutingMuonNumVertices");
  desc.add<std::string>("vertex_index_collection_name", "ScoutingMuonVertexIndex");
  //desc.add<std::string>("num_hit_patterns_name", "ScoutingMuonNumHitPatterns");
  desc.add<std::string>("hit_pattern_collection_name", "ScoutingMuonTrackHitPattern");
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingMuonArrayTableProducer);