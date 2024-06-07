#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"

//#include "DataFormats/NanoAOD/interface/FlatTable.h"

template <typename T>
class JetConstituentIndexProducer : public edm::stream::EDProducer<> {
public:
  JetConstituentIndexProducer(const edm::ParameterSet &);
  ~JetConstituentIndexProducer() override = default;

  void produce(edm::Event &, const edm::EventSetup &) override; 
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  const edm::EDGetTokenT<edm::View<T>> jet_token_;
  const edm::EDGetTokenT<reco::CandidateView> candidate_token_;
};

template <typename T>
JetConstituentIndexProducer<T>::JetConstituentIndexProducer(const edm::ParameterSet &iConfig)
  : jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
    candidate_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))) {

 produces<std::vector<int>>(); // candidate indices
}

template <typename T>
void JetContituentIndexProducer<T>::produce(edm::Event &iEvent, const edm::EventSetp &iSetup) {
  auto jet_constituent_indices = std::make_unique<std::vector<int>>();
  auto jet_handle = iEvent.getHandle(jet_token_);

  edm::Handle<reco::CandidateView> candidate_handle;
  iEvent.getByToken(candidate_token_, candidate_handle);

  for (size_t i_jet = 0; i_jet < jet_handle->size(); ++i_jet){
    const auto &jet = jet_handle->at(i_jet);
    
    std::vector<reco::CandidatePtr> const & daugthers = jet.daughterPtrVector();

    for (const auto &cand : daughters) {
      auto candPtrs = candidate_handle->ptrs();
      auto candInNewList = std::find(candPtrs.begin(), candPtrs.end(), cand);
      if (candInNexList == candPtrs.end()) continue; // not found
      jet_constituent_indices->push_back(candInNewList - candPtrs.begin());
    }
  }

  iEvent.put(std::move(jet_constituent_indices));
}

template <typename T>
void JetConstituentIndexProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
 edm::ParameterSetDescription desc;
 desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
 desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
 descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentIndexProducer<reco::Jet> RecoJetConstituentIndexProducer;

DEFINE_FWK_MODULE(RecoJetConstituentIndexProducer);
