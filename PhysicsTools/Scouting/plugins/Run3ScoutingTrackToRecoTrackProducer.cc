// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class Run3ScoutingTrackToRecoTrackProducer : public edm::stream::EDProducer<> {
public:
    explicit Run3ScoutingTrackToRecoTrackProducer(const edm::ParameterSet &);
    ~Run3ScoutingTrackToRecoTrackProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    void endStream() override {}

    reco::Track createTrack(Run3ScoutingTrack scoutingTrack);

private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingTrack>> input_scoutingTrack_token_; 
};

Run3ScoutingTrackToRecoTrackProducer::Run3ScoutingTrackToRecoTrackProducer(
    const edm::ParameterSet &iConfig)
    : input_scoutingTrack_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingTrack"))){
    
    produces<reco::TrackCollection>();
}

reco::Track Run3ScoutingTrackToRecoTrackProducer::createTrack(Run3ScoutingTrack scoutingTrack){
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

void Run3ScoutingTrackToRecoTrackProducer::produce(edm::Event &iEvent, edm::EventSetup const &setep){
    edm::Handle<std::vector<Run3ScoutingTrack>> scoutingTrackHandle;
    iEvent.getByToken(input_scoutingTrack_token_, scoutingTrackHandle);

    auto recoTrackCollection = std::make_unique<reco::TrackCollection>();

    for (unsigned int itrack = 0; itrack < scoutingTrackHandle->size(); ++itrack) {
        auto &scoutingTrack = (*scoutingTrackHandle)[itrack];
        recoTrackCollection->push_back(createTrack(scoutingTrack));
    }

    iEvent.put(std::move(recoTrackCollection));
}

void Run3ScoutingTrackToRecoTrackProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingTrack", edm::InputTag("hltScoutingTrackPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingTrackToRecoTrackProducer);
