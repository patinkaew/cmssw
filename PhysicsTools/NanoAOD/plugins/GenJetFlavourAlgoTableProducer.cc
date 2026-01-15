// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourAlgoInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourAlgoInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

class GenJetFlavourAlgoTableProducer : public edm::stream::EDProducer<> {
public:
  explicit GenJetFlavourAlgoTableProducer(const edm::ParameterSet& iConfig)
      : name_(iConfig.getParameter<std::string>("name")),
        src_(consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("src"))),
        cut_(iConfig.getParameter<std::string>("cut"), true),
        //findByMatchingDeltaR_(iConfig.getParameter<bool>("findByMatchingDeltaR")),
        deltaR_(iConfig.getParameter<double>("deltaR")) {
    /*const auto& jetFlavourAlgoVPSet = iConfig.getParameter<std::vector<edm::ParameterSet>>("algos");
    for (const auto& jetFlavourAlgoPSet : jetFlavourAlgoVPSet) {
      algoNames_.push_back(jetFlavourAlgoPSet.getParameter<std::string>("name"));
      algoTokens_.push_back(consumes(iConfig.getParameter<edm::InputTag>("tag"));
    }*/
    const auto& algosPSet = iConfig.getParameter<edm::ParameterSet>("algos");
    for (const std::string& algoName : algosPSet.getParameterNamesForType<edm::InputTag>()) {
      algoNames_.push_back(algoName);
      algoTokens_.push_back(consumes<reco::JetFlavourAlgoInfoMatchingCollection>(algosPSet.getParameter<edm::InputTag>(algoName)));
    }
    produces<nanoaod::FlatTable>();
  }

  ~GenJetFlavourAlgoTableProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("input genJet collection");
    desc.add<std::string>("name")->setComment("name of the genJet FlatTable we are extending with flavour information");
    desc.add<std::string>("cut")->setComment("cut on input genJet collection");
    desc.add<double>("deltaR")->setComment("deltaR to match genjets");
    
    edm::ParameterSetDescription algos;
    algos.setComment("input flavour algo info collections");
    algos.addNode(
      edm::ParameterWildcard<edm::InputTag>("*", edm::RequireAtLeastOne, true));
    desc.add<edm::ParameterSetDescription>("algos", algos);

    descriptions.add("genJetFlavourAlgoTable", desc);
  }

private:
  void produce(edm::Event&, edm::EventSetup const&) override;

  std::string name_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > src_;
  const StringCutObjectSelector<reco::GenJet> cut_;
  //const bool findByMatchingDeltaR_;
  const double deltaR_;
  std::vector<std::string> algoNames_;
  std::vector<edm::EDGetTokenT<reco::JetFlavourAlgoInfoMatchingCollection>> algoTokens_;
};

// ------------ method called to produce the data  ------------
void GenJetFlavourAlgoTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& jets = iEvent.getHandle(src_);
  std::vector<edm::Handle<reco::JetFlavourAlgoInfoMatchingCollection>> algos;
  algos.reserve(algoTokens_.size());
  for (const auto& algoToken : algoTokens_) {
    algos.push_back(iEvent.getHandle(algoToken));
  }

  unsigned int ncand = 0;
  unsigned int nalgos = algoNames_.size();
  std::vector<std::vector<int>> flags(nalgos);
  std::vector<std::vector<int>> ds(nalgos);
  std::vector<std::vector<int>> us(nalgos);
  std::vector<std::vector<int>> ss(nalgos);
  std::vector<std::vector<int>> cs(nalgos);
  std::vector<std::vector<int>> bs(nalgos);
  std::vector<std::vector<int>> ts(nalgos);
  std::vector<std::vector<int16_t>> heaviestFlavours(nalgos);
  std::vector<std::vector<int16_t>> heaviestFlavourMod2s(nalgos);

  for (const reco::GenJet& jet : *jets) {
    if (!cut_(jet))
      continue;
    ++ncand;
    for (unsigned int ialgo = 0; ialgo < nalgos; ialgo++) {
      bool matched = false;
      for (const reco::JetFlavourAlgoInfoMatching& jetFlavourAlgoInfoMatching : *(algos[ialgo])) {
        if (deltaR(jet.p4(), jetFlavourAlgoInfoMatching.first->p4()) < deltaR_) {
          flags[ialgo].push_back(jetFlavourAlgoInfoMatching.second.flag()); 
          ds[ialgo].push_back(jetFlavourAlgoInfoMatching.second.d()); 
          us[ialgo].push_back(jetFlavourAlgoInfoMatching.second.u()); 
          ss[ialgo].push_back(jetFlavourAlgoInfoMatching.second.s()); 
          cs[ialgo].push_back(jetFlavourAlgoInfoMatching.second.c()); 
          bs[ialgo].push_back(jetFlavourAlgoInfoMatching.second.b()); 
          ts[ialgo].push_back(jetFlavourAlgoInfoMatching.second.t()); 
          heaviestFlavours[ialgo].push_back(jetFlavourAlgoInfoMatching.second.heaviestFlavour());
          heaviestFlavourMod2s[ialgo].push_back(jetFlavourAlgoInfoMatching.second.heaviestFlavourMod2());
          matched = true;
          break;
        }
      }
      if (!matched) {
        flags[ialgo].push_back(-999);
        ds[ialgo].push_back(-999);
        us[ialgo].push_back(-999);
        ss[ialgo].push_back(-999);
        cs[ialgo].push_back(-999);
        bs[ialgo].push_back(-999);
        ts[ialgo].push_back(-999);
        heaviestFlavours[ialgo].push_back(-999);
        heaviestFlavourMod2s[ialgo].push_back(-999);
      }
    }
  }

  auto tab = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);
  for (unsigned int ialgo = 0; ialgo < nalgos; ialgo++) {
    tab->addColumn<int>(algoNames_[ialgo] + "flag", flags[ialgo], "flag");
    tab->addColumn<int>(algoNames_[ialgo] + "d", ds[ialgo], "d");
    tab->addColumn<int>(algoNames_[ialgo] + "u", us[ialgo], "u");
    tab->addColumn<int>(algoNames_[ialgo] + "s", ss[ialgo], "s");
    tab->addColumn<int>(algoNames_[ialgo] + "c", cs[ialgo], "c");
    tab->addColumn<int>(algoNames_[ialgo] + "b", bs[ialgo], "b");
    tab->addColumn<int>(algoNames_[ialgo] + "t", ts[ialgo], "t");
    tab->addColumn<int16_t>(algoNames_[ialgo] + "heaviestFlavour", heaviestFlavours[ialgo], "heaviestFlavour");
    tab->addColumn<int16_t>(algoNames_[ialgo] + "heaviestFlavourMod2", heaviestFlavourMod2s[ialgo], "heaviestFlavourMod2");
  }

  iEvent.put(std::move(tab));
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJetFlavourAlgoTableProducer);
