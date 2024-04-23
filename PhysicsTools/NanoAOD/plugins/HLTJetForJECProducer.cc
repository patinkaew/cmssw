// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

class HLTJetForJECProducer : public edm::stream::EDProducer<> {
public:
  explicit HLTJetForJECProducer(const edm::ParameterSet &);
  ~HLTJetForJECProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
  void endStream() override {}

private:
  const edm::EDGetTokenT<std::vector<reco::PFJet>> in_corrected_token_;
  const edm::EDGetTokenT<std::vector<reco::PFJet>> in_raw_token_;
  const double ptMin_;
  const double maxDeltaR_;
  std::vector<float> rawFactor_;
};

//
// constructors and destructor
//
HLTJetForJECProducer::HLTJetForJECProducer(
    const edm::ParameterSet &iConfig)
    : in_corrected_token_(consumes(iConfig.getParameter<edm::InputTag>("corrected"))),
      in_raw_token_(consumes(iConfig.getParameter<edm::InputTag>("raw"))), 
      ptMin_(iConfig.getParameter<double>("ptMin")),
      maxDeltaR_(iConfig.getParameter<double>("maxDeltaR"))
      {
  //register products
  produces<reco::PFJetCollection>("corrected");
  produces<reco::PFJetCollection>("raw");
  produces<edm::ValueMap<float>>("rawFactor");
}

HLTJetForJECProducer::~HLTJetForJECProducer() = default;

// ------------ method called to produce the data  ------------
void HLTJetForJECProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup) {
  auto out_corrected = std::make_unique<reco::PFJetCollection>();
  auto out_raw = std::make_unique<reco::PFJetCollection>();

  //try{
    edm::Handle<std::vector<reco::PFJet>> in_corrected_Handle;
    edm::Handle<std::vector<reco::PFJet>> in_raw_Handle;
    iEvent.getByToken(in_corrected_token_, in_corrected_Handle);
    iEvent.getByToken(in_raw_token_, in_raw_Handle);

    //for (unsigned int i_corrected = 0; i_corrected < in_corrected_Handle->size(); ++i_corrected){
    for (const auto& in_corrected : *in_corrected_Handle){
      //auto &in_corrected = (*in_corrected_Handle)[i_corrected];
      if (in_corrected.pt() >= ptMin_){
        for (const auto& in_raw : *in_raw_Handle){
        //for (unsigned int i_raw = 0; i_raw < in_raw_Handle->size(); ++i_raw){
          //auto &in_raw = (*in_raw_Handle)[i_raw];
          if (deltaR(in_corrected.p4(), in_raw.p4()) < maxDeltaR_){
            out_corrected->push_back(in_corrected);
            out_raw->push_back(in_raw);
            rawFactor_.push_back(1 - (in_raw.pt()/in_corrected.pt()));
          }
        }
      } 
    }
  // } catch (const cms::Exception& ex) {
  //     if (ex.category() != "ProductNotFound") throw;
  // }

  // save output
  iEvent.put(std::move(out_raw), "raw");
  edm::OrphanHandle<reco::PFJetCollection> oh = iEvent.put(std::move(out_corrected), "corrected");

  std::unique_ptr<edm::ValueMap<float>> rawFactor_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_rawFactor(*rawFactor_VM);
  filler_rawFactor.insert(oh, rawFactor_.begin(), rawFactor_.end());
  filler_rawFactor.fill();
  iEvent.put(std::move(rawFactor_VM), "rawFactor");

  // clear
  rawFactor_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HLTJetForJECProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("corrected", edm::InputTag("hltPFJetsCorrectedMatchedToCaloJets10"));
  desc.add<edm::InputTag>("raw", edm::InputTag("hltAK4PFJets"));
  desc.add<double>("ptMin", 0);
  desc.add<double>("maxDeltaR", 0.2);
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(HLTJetForJECProducer);
