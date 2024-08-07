#ifndef RecoParticleFlow_PFClusterProducer_PFClusterComparator_
#define RecoParticleFlow_PFClusterProducer_PFClusterComparator_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

/**\class PFClusterAnalyzer 
\brief test analyzer for PFClusters
*/
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

class PFClusterComparator : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit PFClusterComparator(const edm::ParameterSet &);

  ~PFClusterComparator() override = default;

  void analyze(const edm::Event &, const edm::EventSetup &) override;

  void beginRun(const edm::Run &r, const edm::EventSetup &c) override {}
  void endRun(const edm::Run &r, const edm::EventSetup &c) override {}

private:
  void fetchCandidateCollection(edm::Handle<reco::PFClusterCollection> &c,
                                const edm::EDGetTokenT<reco::PFClusterCollection> &token,
                                const edm::Event &iSetup) const;

  /*   void printElementsInBlocks(const reco::PFCluster& cluster, */
  /* 			     std::ostream& out=std::cout) const; */

  /// PFClusters in which we'll look for pile up particles
  const edm::EDGetTokenT<reco::PFClusterCollection> inputTokenPFClusters_;
  const edm::EDGetTokenT<reco::PFClusterCollection> inputTokenPFClustersCompare_;

  /// verbose ?
  const bool verbose_;

  /// print the blocks associated to a given candidate ?
  const bool printBlocks_;

  TH1F *log10E_old, *log10E_new, *deltaEnergy;
  TH1F *posX_old, *posX_new, *deltaX;
  TH1F *posY_old, *posY_new, *deltaY;
  TH1F *posZ_old, *posZ_new, *deltaZ;
};

#endif
