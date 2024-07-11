// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

class Run3ScoutingPFJetToRecoPFJetProducer : public edm::stream::EDProducer<> {
public:
    explicit Run3ScoutingPFJetToRecoPFJetProducer(const edm::ParameterSet &);
    ~Run3ScoutingPFJetToRecoPFJetProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    void endStream() override {}

    reco::PFJet createPFJet(Run3ScoutingPFJet scoutingPFJet);

private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>> input_scoutingPFJet_token_; 
};

Run3ScoutingPFJetToRecoPFJetProducer::Run3ScoutingPFJetToRecoPFJetProducer(
    const edm::ParameterSet &iConfig)
    : input_scoutingPFJet_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingPFJet"))){
    
    produces<reco::PFJetCollection>();
}

reco::PFJet Run3ScoutingPFJetToRecoPFJetProducer::createPFJet(Run3ScoutingPFJet scoutingPFJet){
    // fill LorentzVector P4
    float px = scoutingPFJet.pt() * cos(scoutingPFJet.phi());
    float py = scoutingPFJet.pt() * sin(scoutingPFJet.phi());
    float pz = scoutingPFJet.pt() * sinh(scoutingPFJet.eta());
    float p = scoutingPFJet.pt() * cosh(scoutingPFJet.eta());
    float energy = std::sqrt(p * p + scoutingPFJet.m() * scoutingPFJet.m());
    reco::Particle::LorentzVector p4(px, py, pz, energy);
    
    // fill vertex with default (0, 0, 0)
    reco::Particle::Point vertex(0, 0, 0);
    
    // fill specific
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
    
    // create reco::PFJet output
    reco::PFJet recoPFJet(p4, vertex, specific);

    // set jetArea
    recoPFJet.setJetArea(scoutingPFJet.jetArea());
    
    return recoPFJet;
}

void Run3ScoutingPFJetToRecoPFJetProducer::produce(edm::Event &iEvent, edm::EventSetup const &setep){
    edm::Handle<std::vector<Run3ScoutingPFJet>> scoutingPFJetHandle;
    iEvent.getByToken(input_scoutingPFJet_token_, scoutingPFJetHandle);

    auto outPFJetCollection = std::make_unique<reco::PFJetCollection>();

    for (unsigned int ijet = 0; ijet < scoutingPFJetHandle->size(); ++ijet) {
        auto &scoutingPFJet = (*scoutingPFJetHandle)[ijet];
        outPFJetCollection->push_back(createPFJet(scoutingPFJet));
    }

    iEvent.put(std::move(outPFJetCollection));
}

void Run3ScoutingPFJetToRecoPFJetProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingPFJet", edm::InputTag("hltScoutingPFPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingPFJetToRecoPFJetProducer);
