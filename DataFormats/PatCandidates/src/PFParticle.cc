//
//

#include "DataFormats/PatCandidates/interface/PFParticle.h"

using namespace pat;

/// constructor from reco::PFCandidate
PFParticle::PFParticle(const reco::PFCandidate& aPFParticle)
    : PATObject<reco::PFCandidate>(aPFParticle) {}

/// constructor from PFParticleType
PFParticle::PFParticle(const edm::RefToBase<reco::PFCandidate>& aPFParticle)
    : PATObject<reco::PFCandidate>(aPFParticle) {}
