#include <memory>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

//#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
//#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
//#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "fastjet/contrib/SoftKiller.hh"

/*
 * HLTScoutingUnpackProducer unpacks Run3Scouting data formats to their reco format counterparts
 * and cross link collections (currently PFCandidate - Track)
*/

class HLTScoutingUnpackProducer : public edm::stream::EDProducer<> {
public:
    //using Run3ScoutingElectron = std::vector<Run3ScoutingElectron>;
    //using Run3ScoutingMuon = std::vector<Run3ScoutingMuon>;
    using Run3ScoutingPFJetCollection = std::vector<Run3ScoutingPFJet>;
    using Run3ScoutingParticleCollection = std::vector<Run3ScoutingParticle>;
    //using Run3ScoutingPhotonCollection = std::vector<Run3ScoutingPhoton>;
    using Run3ScoutingTrackCollection = std::vector<Run3ScoutingTrack>;
    using Run3ScoutingVertexCollection = std::vector<Run3ScoutingVertex>;
    
    template <typename T> using RefCollection = std::vector<edm::Ref<std::vector<T>>>;
    template <typename T> using RefMap = edm::ValueMap<edm::Ref<std::vector<T>>>;

    explicit HLTScoutingUnpackProducer(edm::ParameterSet const& params);
    ~HLTScoutingUnpackProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
private: 
    void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;
    
    // private helper functions
    reco::PFJet createPFJet(Run3ScoutingPFJet const& scoutingPFJet);
    reco::Vertex createVertex(Run3ScoutingVertex const& scoutingVertex);
    reco::Track createTrack(Run3ScoutingTrack const& scoutingTrack);
    reco::PFCandidate createPFCandidate(Run3ScoutingParticle const& scoutingPFCandidate);
    std::unique_ptr<reco::Track> createPFTrack(Run3ScoutingParticle const& scoutingPFCandidate, edm::Handle<Run3ScoutingVertexCollection> const& scoutingPrimaryVertex_collection_handle);
    
    // register products for reco object and corresponding Value to Ref to original scouting objects
    template<typename RecoObjectType, typename ScoutingObjectType>
    void produceWithRef(std::string const& name);

    // put reco object together with ValueMap to Ref to original scouting object
    template<typename RecoObjectType, typename ScoutingObjectType>
    void putWithRef(edm::Event& iEvent, std::string const& name, std::unique_ptr<std::vector<RecoObjectType>>& recoObject_collection_ptr, std::unique_ptr<RefCollection<ScoutingObjectType>>& scoutingObjectRef_collection_ptr);

    //edm::EDGetTokenT<Run3ScoutingElectronCollection> scoutingElectron_collection_token_;
    //edm::EDGetTokenT<Run3ScoutingMuonCollection> scoutingMuon_collection_token_;
    edm::EDGetTokenT<Run3ScoutingPFJetCollection> scoutingPFJet_collection_token_;
    edm::EDGetTokenT<Run3ScoutingParticleCollection> scoutingPFCandidate_collection_token_;
    //edm::EDGetTokenT<Run3ScoutingPhotonCollection> scoutingPhoton_collection_token_;
    edm::EDGetTokenT<Run3ScoutingTrackCollection> scoutingTrack_collection_token_;
    edm::EDGetTokenT<Run3ScoutingVertexCollection> scoutingPrimaryVertex_collection_token_;

    edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particle_data_table_token_;

    bool produce_PFCandidate_;
    bool produce_PFCHSCandidate_; // CHS = charged hadron subtraction
    bool produce_PFSKCandidate_; // SK = soft killer
    HepPDT::ParticleDataTable const *particle_data_table_;

    //edm::EDPutTokenT<reco::TrackCollection> recoTrack_collection_token_;
    //edm::EDPutTokenT<reco::TrackCollection> recoPFTrack_collection_token_;
    
    inline static const std::string REF_TO_SCOUTING_LABEL_SUFFIX_ = "-RefToScouting"; 
};

// constructor
HLTScoutingUnpackProducer::HLTScoutingUnpackProducer(edm::ParameterSet const& params)
    : //scoutingElectron_token_(consumes(params.getParameter<edm::InputTag>("scoutingElectron"))),
      //scoutingMuon_token_(consumes(params.getParameter<edm::InputTag>("scoutingMuon"))),
      scoutingPFJet_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingPFJet"))),
      scoutingPFCandidate_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingPFCandidate"))),
      //scoutingPhoton_token_(consumes(params.getParameter<edm::InputTag>("scoutingPhoton"))),
      scoutingTrack_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingTrack"))),
      scoutingPrimaryVertex_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingPrimaryVertex"))),
      particle_data_table_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()),
      produce_PFCandidate_(params.getParameter<bool>("producePFCandidate")),
      produce_PFCHSCandidate_(params.getParameter<bool>("producePFCHSCandidate")),
      produce_PFSKCandidate_(params.getParameter<bool>("producePFSKCandidate")){
    
    produces<reco::PFJetCollection>("PFJet");
    //produces<edm::ValueMap<edm::Ref<Run3ScoutingPFJetCollection>>>("PFJet_ScoutingRef");
    
    //recoTrack_collection_token_ = produces<reco::TrackCollection>("Track");
    //produces<reco::TrackCollection>("Track");
    //produces<edm::ValueMap<edm::Ref<Run3ScoutingTrackCollection>>>("Track_ScoutingRef");
    produceWithRef<reco::Track, Run3ScoutingTrack>("Track");

    //produces<reco::VertexCollection>("PrimaryVertex");
    //produces<Run3ScoutingVertexRefMap>("PrimaryVertex_RefToScouting");
    produceWithRef<reco::Vertex, Run3ScoutingVertex>("PrimaryVertex");
    
    if (produce_PFCandidate_){
        produces<reco::PFCandidateCollection>("PFCandidate");
        //produces<edm::ValueMap<edm::Ref<Run3ScoutingParticleCollection>>>("PFCandidate_ScoutingRef");
        //produces<reco::ValueMap<int>>("PFCandidate_vertexIndex");
        //producePFCandidateCollection("PFCandidate");
    }
    if (produce_PFCHSCandidate_){
        produces<reco::PFCandidateCollection>("PFCHSCandidate");
    }
    /*if (produce_PFSKCandidate_){
        produces<reco::PFCandidateCollection>("PFSKCandidate");
    }*/
    if (produce_PFCandidate_ || produce_PFCHSCandidate_ || produce_PFSKCandidate_) {
        //recoPFTrack_collection_token_ = produces<reco::TrackCollection>("PFTrack"); // reco::track collection from best track parameters of PF
        produces<reco::TrackCollection>("PFTrack"); // reco::track collection from best track parameters of PF
    }
}

template <typename RecoObjectType, typename ScoutingObjectType>
void HLTScoutingUnpackProducer::produceWithRef(std::string const& name) {
    produces<std::vector<RecoObjectType>>(name);
    produces<RefMap<ScoutingObjectType>>(name + REF_TO_SCOUTING_LABEL_SUFFIX_);
}

void HLTScoutingUnpackProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
    // produce reco::PFJet
    edm::Handle<Run3ScoutingPFJetCollection> scoutingPFJet_collection_handle = iEvent.getHandle(scoutingPFJet_collection_token_);
    auto recoPFJet_collection_ptr = std::make_unique<reco::PFJetCollection>();
    //auto scoutingPFJetRef_collection_ptr = std::make_unique<RefCollection<Run3ScoutingPFJet>>();
    if (scoutingPFJet_collection_handle.isValid()) {
        //for (size_t scoutingPFJet_index=0; scoutingPFJet_index < scoutingPFJet_collection_handle->size(); scoutingPFJet_index++) {
        for (auto const& scoutingPFJet: *scoutingPFJet_collection_handle) {
            recoPFJet_collection_ptr->push_back(createPFJet(scoutingPFJet));
            //edm::Ref<Run3ScoutingPFJetCollection>(scoutingPFJet_collection_handle, i);
        }
    }

    // produce reco::Vertex
    edm::Handle<Run3ScoutingVertexCollection> scoutingPrimaryVertex_collection_handle = iEvent.getHandle(scoutingPrimaryVertex_collection_token_);
    auto recoPrimaryVertex_collection_ptr = std::make_unique<reco::VertexCollection>();
    auto scoutingPrimaryVertexRef_collection_ptr = std::make_unique<RefCollection<Run3ScoutingVertex>>();
    if (scoutingPrimaryVertex_collection_handle.isValid()) {
        for (size_t scoutingPrimaryVertex_index=0; scoutingPrimaryVertex_index < scoutingPrimaryVertex_collection_handle->size(); scoutingPrimaryVertex_index++) {
        //for (auto const& scoutingPrimaryVertex: *scoutingPrimaryVertex_collection_handle) {
            auto &scoutingPrimaryVertex = scoutingPrimaryVertex_collection_handle->at(scoutingPrimaryVertex_index);
            recoPrimaryVertex_collection_ptr->push_back(createVertex(scoutingPrimaryVertex));
            scoutingPrimaryVertexRef_collection_ptr->push_back(edm::Ref<Run3ScoutingVertexCollection>(scoutingPrimaryVertex_collection_handle, scoutingPrimaryVertex_index));
        }
    }

    // produce reco::Track
    edm::Handle<Run3ScoutingTrackCollection> scoutingTrack_collection_handle = iEvent.getHandle(scoutingTrack_collection_token_);
    auto recoTrack_collection_ptr = std::make_unique<reco::TrackCollection>();
    auto scoutingTrackRef_collection_ptr = std::make_unique<RefCollection<Run3ScoutingTrack>>();
    if (scoutingTrack_collection_handle.isValid()) {
        for (size_t scoutingTrack_index=0; scoutingTrack_index < scoutingTrack_collection_handle->size(); scoutingTrack_index++) {
            auto &scoutingTrack = scoutingTrack_collection_handle->at(scoutingTrack_index);
        //for (auto const& scoutingTrack: *scoutingTrack_collection_handle) {
            recoTrack_collection_ptr->push_back(createTrack(scoutingTrack));
            scoutingTrackRef_collection_ptr->push_back(edm::Ref<Run3ScoutingTrackCollection>(scoutingTrack_collection_handle, scoutingTrack_index));
        }
    }

    // produce reco::PFCandidate
    edm::Handle<Run3ScoutingParticleCollection> scoutingPFCandidate_collection_handle = iEvent.getHandle(scoutingPFCandidate_collection_token_);
    auto recoPFCandidate_collection_ptr = std::make_unique<reco::PFCandidateCollection>();
    auto recoPFCHSCandidate_collection_ptr = std::make_unique<reco::PFCandidateCollection>();
    auto recoPFTrack_collection_ptr = std::make_unique<reco::TrackCollection>();
    //reco::TrackRefProd recoTrack_collection_refprod = iEvent.getRefBeforePut(recoTrack_collection_token_);
    //reco::TrackRefProd recoPFTrack_collection_refprod = iEvent.getRefBeforePut(recoPFTrack_collection_token_);
    if (scoutingPFCandidate_collection_handle.isValid() && (produce_PFCandidate_ || produce_PFCHSCandidate_ || produce_PFSKCandidate_)) {
        reco::TrackRefProd recoTrack_collection_refProd = iEvent.getRefBeforePut<reco::TrackCollection>("Track");
        reco::TrackRefProd recoPFTrack_collection_refProd = iEvent.getRefBeforePut<reco::TrackCollection>("PFTrack");
        particle_data_table_ = iSetup.getHandle(particle_data_table_token_).product();
        for (auto const& scoutingPFCandidate: *scoutingPFCandidate_collection_handle) {    
            reco::PFCandidate recoPFCandidate = createPFCandidate(scoutingPFCandidate);
            std::unique_ptr<reco::Track> recoPFTrack_ptr = createPFTrack(scoutingPFCandidate, scoutingPrimaryVertex_collection_handle);
            // build track and vertex for output reco::PFCandidate
            if (recoPFTrack_ptr != nullptr) {
                // try to search for track built from ScoutingTrack containing more information
                int recoTrack_index = -1; //findCompatibleScoutingTrack(recoPFTrack_ptr, recoTrack_collection_ptr);
                if (recoTrack_index >= 0) { // found
                    reco::TrackRef trackRef(recoTrack_collection_refProd, recoTrack_index);
                    recoPFCandidate.setTrackRef(trackRef);
                    recoPFCandidate.setVertex((recoTrack_collection_ptr->at(recoTrack_index)).vertex());
                } else { // not found 
                    recoPFTrack_collection_ptr->push_back(*recoPFTrack_ptr);
                    reco::TrackRef trackRef(recoPFTrack_collection_refProd, recoPFTrack_collection_ptr->size() - 1);
                    recoPFCandidate.setTrackRef(trackRef);
                    recoPFCandidate.setVertex(recoPFTrack_ptr->vertex());
                }
            }
            recoPFCandidate_collection_ptr->push_back(recoPFCandidate);
            if (produce_PFCHSCandidate_ && (scoutingPFCandidate.vertex() >= 0)) {
                recoPFCHSCandidate_collection_ptr->push_back(recoPFCandidate);
            }
        }
    }

    // put products back in Event
    iEvent.put(std::move(recoPFJet_collection_ptr), "PFJet");

    //iEvent.put(std::move(recoPrimaryVertex_collection_ptr), "PrimaryVertex");
    putWithRef<reco::Vertex, Run3ScoutingVertex>(iEvent, "PrimaryVertex", recoPrimaryVertex_collection_ptr, scoutingPrimaryVertexRef_collection_ptr);
    
    //iEvent.put(std::move(recoTrack_collection_ptr), "Track");
    putWithRef<reco::Track, Run3ScoutingTrack>(iEvent, "Track", recoTrack_collection_ptr, scoutingTrackRef_collection_ptr);
    
    if (produce_PFCandidate_ || produce_PFCHSCandidate_ || produce_PFSKCandidate_) {
        //std::cout << "PFTrack put " << produce_PFCandidate_ << std::endl;
        iEvent.put(std::move(recoPFTrack_collection_ptr), "PFTrack");
    }
    if (produce_PFCandidate_) {
        iEvent.put(std::move(recoPFCandidate_collection_ptr), "PFCandidate");
    }
    if (produce_PFCHSCandidate_) {
        iEvent.put(std::move(recoPFCHSCandidate_collection_ptr), "PFCHSCandidate");
    }
    
}

reco::PFJet HLTScoutingUnpackProducer::createPFJet(Run3ScoutingPFJet const& scoutingPFJet) {
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

reco::Vertex HLTScoutingUnpackProducer::createVertex(Run3ScoutingVertex const& scoutingVertex) {
    // fill point coordinate
    reco::Vertex::Point point(scoutingVertex.x(), scoutingVertex.y(), scoutingVertex.z());
    
    // fill error
    std::vector<float> error_vec(6);
    error_vec[0] = scoutingVertex.xError() * scoutingVertex.xError(); // cov(0, 0)
    error_vec[1] = 0; // cov(0, 1)
    error_vec[2] = 0; // cov(0, 2)
    error_vec[3] = scoutingVertex.yError() * scoutingVertex.yError(); // cov(1, 1)
    error_vec[4] = 0; // cov(1, 2)
    error_vec[5] = scoutingVertex.zError() * scoutingVertex.zError(); // cov(2, 2)

    // off-diagonal errors are added in the begining of 2024
    // see https://github.com/cms-sw/cmssw/pull/43758
    try { 
      error_vec[1] = scoutingVertex.xyCov();
      error_vec[2] = scoutingVertex.xzCov();
      error_vec[4] = scoutingVertex.yzCov();
    }
    catch (...) { // do nothing
    }

    reco::Vertex::Error error(error_vec.begin(), error_vec.end());

    return scoutingVertex.isValidVtx() ? reco::Vertex(point, error, scoutingVertex.chi2(), scoutingVertex.ndof(), scoutingVertex.tracksSize()): reco::Vertex(point, error);

}

reco::Track HLTScoutingUnpackProducer::createTrack(Run3ScoutingTrack const& scoutingTrack){
    float chi2 = scoutingTrack.tk_chi2();
    float ndof = scoutingTrack.tk_ndof();
    reco::TrackBase::Point referencePoint(scoutingTrack.tk_vx(), scoutingTrack.tk_vy(), scoutingTrack.tk_vz());
    
    float px = scoutingTrack.tk_pt() * cos(scoutingTrack.tk_phi());
    float py = scoutingTrack.tk_pt() * sin(scoutingTrack.tk_phi());
    float pz = scoutingTrack.tk_pt() * sinh(scoutingTrack.tk_eta());
    reco::TrackBase::Vector momentum(px, py, pz);
    
    int charge = scoutingTrack.tk_charge();
    
    std::vector<float> cov_vec(15); // 5*(5+1)/2 = 15
    cov_vec[0] = scoutingTrack.tk_qoverp_Error() * scoutingTrack.tk_qoverp_Error(); // cov(0, 0)
    cov_vec[1] = scoutingTrack.tk_qoverp_lambda_cov(); // cov(0, 1)
    cov_vec[2] = scoutingTrack.tk_qoverp_phi_cov(); // cov(0, 2)
    cov_vec[3] = scoutingTrack.tk_qoverp_dxy_cov(); // cov(0, 3)
    cov_vec[4] = scoutingTrack.tk_qoverp_dsz_cov(); // cov(0, 4)
    cov_vec[5] = scoutingTrack.tk_lambda_Error() * scoutingTrack.tk_lambda_Error(); // cov(1, 1)
    cov_vec[6] = scoutingTrack.tk_lambda_phi_cov(); // cov(1, 2)
    cov_vec[7] = scoutingTrack.tk_lambda_dxy_cov(); // cov(1, 3)
    cov_vec[8] = scoutingTrack.tk_lambda_dsz_cov(); // cov(1, 4)
    cov_vec[9] = scoutingTrack.tk_phi_Error() * scoutingTrack.tk_phi_Error(); // cov(2, 2) 
    cov_vec[10] = scoutingTrack.tk_phi_dxy_cov(); // cov(2, 3)
    cov_vec[11] = scoutingTrack.tk_phi_dsz_cov(); // cov(2, 4)
    cov_vec[12] = scoutingTrack.tk_dxy_Error() * scoutingTrack.tk_dxy_Error(); // cov(3, 3)
    cov_vec[13] = scoutingTrack.tk_dxy_dsz_cov(); // cov(3, 4)
    cov_vec[14] = scoutingTrack.tk_dsz_Error() * scoutingTrack.tk_dsz_Error(); // cov(4, 4)
    reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

    reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm); // undefined
    reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed); // confirmed
    
    // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1 
  
    reco::Track recoTrack(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);

    return recoTrack;
}

reco::PFCandidate HLTScoutingUnpackProducer::createPFCandidate(Run3ScoutingParticle const& scoutingPFCandidate) {
    int pdgId = scoutingPFCandidate.pdgId();
    // fill mass and charge from pdgId
    float m = 0.;
    int charge = 0;
    auto particle_data_ptr = particle_data_table_->particle(HepPDT::ParticleID(pdgId)); // particle data
    // fake pdgId: h_HF (1), egamma_HF (2), X(0)
    if (!(pdgId == 1 || pdgId == 2 || pdgId == 0) && (particle_data_ptr != nullptr)) {
        m = particle_data_ptr->mass();
        charge = particle_data_ptr->charge();
    } else {
        return reco::PFCandidate();
    }
    
    // fill 4-momentum
    float px = scoutingPFCandidate.pt() * cos(scoutingPFCandidate.phi());
    float py = scoutingPFCandidate.pt() * sin(scoutingPFCandidate.phi());
    float pz = scoutingPFCandidate.pt() * sinh(scoutingPFCandidate.eta());
    float p = scoutingPFCandidate.pt() * cosh(scoutingPFCandidate.eta());
    float energy = std::sqrt(p * p + m * m);
    reco::Particle::LorentzVector p4(px, py, pz, energy);
    
    static const reco::PFCandidate dummy;
    return reco::PFCandidate(charge, p4, dummy.translatePdgIdToType(pdgId));
}

std::unique_ptr<reco::Track> HLTScoutingUnpackProducer::createPFTrack(Run3ScoutingParticle const& scoutingPFCandidate,
                                                                      edm::Handle<Run3ScoutingVertexCollection> const& scoutingPrimaryVertex_collection_handle) {
    // retrive pdgId
    int pdgId = scoutingPFCandidate.pdgId();
    // no track for photon(22), neutral hadron(130), h_HF (1), egamma_HF (2), X(0)
    // see https://github.com/cms-sw/cmssw/blob/master/DataFormats/ParticleFlowCandidate/src/PFCandidate.cc#L231-L276
    if (pdgId == 22 || pdgId == 130 || pdgId == 1 || pdgId == 2 || pdgId == 0) {
      return nullptr;
    }

    // set track charge from PFCandidate's charge
    int charge;
    auto particle_data_ptr = particle_data_table_->particle(HepPDT::ParticleID(pdgId)); // particle data
    if ( particle_data_ptr != nullptr) {
        charge = particle_data_ptr->charge();
    } else {
        return nullptr;
    }

    // fill chi2 and ndof
    // only save normalized chi2 = chi2/ndof is saved, assume ndof = 1
    float chi2 = scoutingPFCandidate.normchi2();
    int ndof = 1;

    if (chi2 == 999) {
        return nullptr;
    }

    // retrieve track parameters
    float pt = scoutingPFCandidate.trk_pt();
    float eta = scoutingPFCandidate.trk_eta();
    float phi = scoutingPFCandidate.trk_phi();
    if (scoutingPFCandidate.relative_trk_vars()){
        pt += scoutingPFCandidate.pt();
        eta += scoutingPFCandidate.eta();
        phi += scoutingPFCandidate.phi();
    }
  
    // fill momentum
    float px = pt * cos(phi);
    float py = pt * sin(phi);
    float pz = pt * sinh(eta);
    float p = pt * cosh(eta);
    reco::TrackBase::Vector momentum(px, py, pz);
    
    // fill 5D curvilinear covariance matrix
    // only dxysig and dzsig are stored
    std::vector<float> cov_vec(15, 0); // number of element is 5*(5+1)/2 = 15, default value to 0
    float dxy = scoutingPFCandidate.dxy();
    float dz = scoutingPFCandidate.dz();
    float dxyError = dxy / scoutingPFCandidate.dxysig();
    float dzError = dz / scoutingPFCandidate.dzsig();
    float dszError = dzError * pt/p;
    cov_vec[12] = dxyError * dxyError;
    cov_vec[14] = dszError * dszError;
    reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

    reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm); // undefined
    reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed); // confirmed FIXME: change to use scoutingPFCandidate.quality()

    // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1

    // fill reference point
    // reference point is not stored, try PFCandidate's vertex, otherwise default to closest point to PV0 
    Run3ScoutingVertex pv0 = scoutingPrimaryVertex_collection_handle->at(0);
    reco::TrackBase::Point pv0_point(pv0.x(), pv0.y(), pv0.z());
    int vertex_index = scoutingPFCandidate.vertex();
    // try using PFCandidate's vertex as referencePoint
    if (vertex_index >= 0) {
        Run3ScoutingVertex pv = scoutingPrimaryVertex_collection_handle->at(vertex_index);
        reco::TrackBase::Point referencePoint(pv.x(), pv.y(), pv.z());
        auto recoTrack_ptr = std::make_unique<reco::Track>(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);
        float epsilon = 0.0001;
        if (fabs(recoTrack_ptr->dxy(pv0_point) - scoutingPFCandidate.dxy()) < epsilon
            && fabs(recoTrack_ptr->dxy(pv0_point)/recoTrack_ptr->dxyError() - scoutingPFCandidate.dxysig()) < epsilon
            && fabs(recoTrack_ptr->dz(pv0_point) - scoutingPFCandidate.dz()) < epsilon
            && fabs(recoTrack_ptr->dz(pv0_point)/recoTrack_ptr->dzError() - scoutingPFCandidate.dzsig()) < epsilon) {
            return recoTrack_ptr;
        }
    }
    // estimate to closest point to PV0
    float vx = pv0.x() - fabs(dxy)*sin(phi);
    float vy = pv0.y() + fabs(dxy)*cos(phi);
    float vz = pv0.z() + dz;
    reco::TrackBase::Point referencePoint(vx, vy, vz);
    return std::make_unique<reco::Track>(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);
}

template <typename RecoObjectType, typename ScoutingObjectType>
void HLTScoutingUnpackProducer::putWithRef(edm::Event& iEvent, std::string const& name,
                                           std::unique_ptr<std::vector<RecoObjectType>>& recoObject_collection_ptr, std::unique_ptr<RefCollection<ScoutingObjectType>>& scoutingObjectRef_collection_ptr) {
    auto recoObject_collection_handle = iEvent.put(std::move(recoObject_collection_ptr), name);

    std::unique_ptr<RefMap<ScoutingObjectType>> refmap_to_scouting(new RefMap<ScoutingObjectType>());
    typename RefMap<ScoutingObjectType>::Filler filler_refmap_to_scouting(*refmap_to_scouting);
    filler_refmap_to_scouting.insert(recoObject_collection_handle, scoutingObjectRef_collection_ptr->begin(), scoutingObjectRef_collection_ptr->end());
    filler_refmap_to_scouting.fill();
    iEvent.put(std::move(refmap_to_scouting), name + REF_TO_SCOUTING_LABEL_SUFFIX_);
}

void HLTScoutingUnpackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("scoutingPFJet", edm::InputTag("hltScoutingPFPacker"));
    desc.add<edm::InputTag>("scoutingPFCandidate", edm::InputTag("hltScoutingPFPacker"));
    desc.add<edm::InputTag>("scoutingTrack", edm::InputTag("hltScoutingTrack"));
    desc.add<edm::InputTag>("scoutingPrimaryVertex", edm::InputTag("hltScoutingPrimaryVertex"));
    
    desc.add<bool>("producePFCandidate", true);
    desc.add<bool>("producePFCHSCandidate", false);
    desc.add<bool>("producePFSKCandidate", false);

    descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HLTScoutingUnpackProducer);
