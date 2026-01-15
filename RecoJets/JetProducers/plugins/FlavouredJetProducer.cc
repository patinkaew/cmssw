#include "RecoJets/JetProducers/plugins/FlavouredJetProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <memory>

FlavouredJetProducerBase::FlavouredJetProducerBase(const edm::ParameterSet& iConfig)
    : FastjetJetProducer(iConfig),
      recombination_(iConfig.getParameter<std::string>("recombination")) {
  // this only works with GenJet
  if (jetType_ != "GenJet") {
    throw cms::Exception("FlavouredJetProducer", "Flavoured jet algorithms only works with GenJet, but JetType is " + jetType_);
  }
  
  // set get flavour recombiner
  if (recombination_ == "net") {
    fjRecombiner_ = std::make_shared<fastjet::contrib::FlavRecombiner>(fastjet::contrib::FlavRecombiner::FlavSummation::net);
  } else if (recombination_ == "mod2" || recombination_ == "modulo_2") {
    fjRecombiner_ = std::make_shared<fastjet::contrib::FlavRecombiner>(fastjet::contrib::FlavRecombiner::FlavSummation::modulo_2);
  } else if (recombination_ == "any" || recombination_ == "any_abs") {
    fjRecombiner_ = std::make_shared<fastjet::contrib::FlavRecombiner>(fastjet::contrib::FlavRecombiner::FlavSummation::any_abs);
  } else {
    throw cms::Exception("Configuration", "Unknown flavour recombination scheme: " + recombination_);
  }
  fjJetDefinition_->set_recombiner(fjRecombiner_.get());

  produces<reco::JetFlavourAlgoInfoMatchingCollection>(jetCollInstanceName_);
}

void FlavouredJetProducerBase::fillDescriptionsFromFlavouredJetProducerBase(edm::ParameterSetDescription& desc) {
  // From FastjetJetProducer
  FastjetJetProducer::fillDescriptionsFromFastJetProducer(desc);

  // From VirtualJetProducer
  VirtualJetProducer::fillDescriptionsFromVirtualJetProducer(desc);

  desc.add<std::string>("recombination", "net");
}

void FlavouredJetProducerBase::inputTowers() {
  // first call inputTowerts from FastjetJetProducer (actually from VirtualJetProducer)
  FastjetJetProducer::inputTowers();

  // add pdgId
  for (auto& fjInput : fjInputs_) {
    auto index = fjInput.user_index();
    auto pdgId = inputs_[index]->pdgId();
    fjInput.set_user_info(new fastjet::contrib::FlavHistory(pdgId));
  }
}

void FlavouredJetProducerBase::output(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  FastjetJetProducer::output(iEvent, iSetup);
  
  auto jetFlavourAlgoInfos = std::make_unique<reco::JetFlavourAlgoInfoMatchingCollection>(reco::JetRefBaseProd(genJetHandle_));

  for (unsigned int ijet = 0; ijet < fjJets_.size(); ijet++) {
    auto flav = fastjet::contrib::FlavHistory::current_flavour_of(fjJets_[ijet]);
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<std::vector<reco::GenJet>>(genJetHandle_, ijet));
    (*jetFlavourAlgoInfos)[jetRef] = reco::JetFlavourAlgoInfo(flav._flav_content);
  }

  iEvent.put(std::move(jetFlavourAlgoInfos), jetCollInstanceName_);
}

////////////////////////////
// SDF (SoftDrop Flavour) //
////////////////////////////

SDFJetProducer::SDFJetProducer(const edm::ParameterSet& iConfig)
    : FlavouredJetProducerBase(iConfig) {
  /*
  bool mod2 = false;

  if (recombination_ == "net") {
    mod2 = false;
  } else if (recombination_ == "mod2" || recombination_ == "modulo_2") {
    mod2 = true;
  } else if (recombination_ == "any" || recombination_ == "any_abs") {
    throw cms::Exception("Configuration", "SDF does not support any flavour recombination scheme");
  } else {
    throw cms::Exception("Configuration", "Unknown flavour recombination scheme: " + recombination_);
  }
  */

  sdFlavCalc_ = fastjet::contrib::SDFlavourCalc(
      iConfig.getParameter<double>("beta"),
      iConfig.getParameter<double>("zcut"),
      iConfig.getParameter<double>("R0"));
}

void SDFJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  fillDescriptionsFromFlavouredJetProducerBase(desc);
  
  desc.add<std::string>("jetCollInstanceName", "");
  descriptions.add("SDFJetProducer", desc);
}

void SDFJetProducer::runAlgorithm(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  FastjetJetProducer::runAlgorithm(iEvent, iSetup);
  
  sdFlavCalc_(fjJets_);
}

DEFINE_FWK_MODULE(SDFJetProducer);

/////////////////////////////
// CMP (Flavoured anti-kt) //
/////////////////////////////

CMPJetProducer::CMPJetProducer(const edm::ParameterSet& iConfig)
    : FlavouredJetProducerBase(iConfig),
      a_(iConfig.getParameter<double>("a")),
      correction_type_(iConfig.getParameter<std::string>("correctionType")),
      clustering_type_(iConfig.getParameter<std::string>("clusteringType")) {
  
  fastjet::CMPPlugin::CorrectionType corr;
  if (correction_type_ == "NoCorrection" ) {
      corr = fastjet::CMPPlugin::CorrectionType::NoCorrection;
  } else if (correction_type_ == "SqrtCoshyCosPhiArgument") {
      corr = fastjet::CMPPlugin::CorrectionType::SqrtCoshyCosPhiArgument;
  } else if (correction_type_ == "SqrtCoshyCosPhiArgument_a2") {
      corr = fastjet::CMPPlugin::CorrectionType::SqrtCoshyCosPhiArgument_a2;
  } else if (correction_type_ == "CoshyCosPhi") {
      corr = fastjet::CMPPlugin::CorrectionType::CoshyCosPhi;
  } else if (correction_type_ == "OverAllCoshyCosPhi") {
      corr = fastjet::CMPPlugin::CorrectionType::OverAllCoshyCosPhi;
  } else if (correction_type_ == "OverAllCoshyCosPhi_a2") {
      corr = fastjet::CMPPlugin::CorrectionType::OverAllCoshyCosPhi_a2;
  } else {
      throw cms::Exception("Configuration", "Unknown correction type for CMP: " + correction_type_);
  }

  fastjet::CMPPlugin::ClusteringType clust;
  if (clustering_type_ == "DynamicKtMax") {
    clust = fastjet::CMPPlugin::ClusteringType::DynamicKtMax;
  } else if (clustering_type_ == "FixedKtMax") {
      clust = fastjet::CMPPlugin::ClusteringType::FixedKtMax;
  } else {
      throw cms::Exception("Configuration", "Unknown clustering type for CMP: " + clustering_type_);
  }

  fjCMPPlugin_ = std::make_shared<fastjet::CMPPlugin>(rParam_, a_, corr, clust);
  fjJetDefinition_ = std::make_shared<fastjet::JetDefinition>(fjCMPPlugin_.get());
  fjJetDefinition_->set_recombiner(fjRecombiner_.get());
}

void CMPJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  fillDescriptionsFromFlavouredJetProducerBase(desc);
  desc.add<double>("a", 0.1);
  desc.add<std::string>("correctionType", "SqrtCoshyCosPhiArgument_a2");
  desc.add<std::string>("clusteringType", "DynamicKtMax");

  desc.add<std::string>("jetCollInstanceName", "");
  descriptions.add("CMPJetProducer", desc);
}

DEFINE_FWK_MODULE(CMPJetProducer);

//////////////////////////////////////////////
// IFN (Interleaved Flavour Neutralisation) //
//////////////////////////////////////////////

IFNJetProducer::IFNJetProducer(const edm::ParameterSet& iConfig)
    : FlavouredJetProducerBase(iConfig),
      alpha_(iConfig.getParameter<double>("alpha")),
      omega_(iConfig.getParameter<double>("omega")) {
  // save original fjJetDefinition_as fjBaseJetDefinition
  fjBaseJetDefinition_ = fjJetDefinition_;

  // set IFNPlugin
  if (recombination_ == "net") {
    fjIFNPlugin_ = std::make_shared<fastjet::contrib::IFNPlugin>(*fjBaseJetDefinition_, alpha_, omega_, fastjet::contrib::FlavRecombiner::FlavSummation::net);
  } else if (recombination_ == "mod2" || recombination_ == "modulo_2") {
    fjIFNPlugin_ = std::make_shared<fastjet::contrib::IFNPlugin>(*fjBaseJetDefinition_, alpha_, omega_, fastjet::contrib::FlavRecombiner::FlavSummation::modulo_2);
  } else if (recombination_ == "any" || recombination_ == "any_abs") {
    fjIFNPlugin_ = std::make_shared<fastjet::contrib::IFNPlugin>(*fjBaseJetDefinition_, alpha_, omega_, fastjet::contrib::FlavRecombiner::FlavSummation::any_abs);
  } else {
    throw cms::Exception("Configuration", "Unknown flavour recombination scheme: " + recombination_);
  }
  
  // redefine fjJetDefinition_
  fjJetDefinition_ = std::make_shared<fastjet::JetDefinition>(fjIFNPlugin_.get());
}

void IFNJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  
  fillDescriptionsFromFlavouredJetProducerBase(desc);
  
  desc.add<double>("alpha", 1.0);
  desc.add<double>("omega", 2.0);

  desc.add<std::string>("jetCollInstanceName", "");
  descriptions.add("IFNJetProducer", desc);
}

DEFINE_FWK_MODULE(IFNJetProducer);

////////////////////////////
// GHS (Flavour Dressing) //
///////////////////////////

GHSJetProducer::GHSJetProducer(const edm::ParameterSet& iConfig)
    : FlavouredJetProducerBase(iConfig),
      ptcut_(iConfig.getParameter<double>("ptcut")),
      alpha_(iConfig.getParameter<double>("alpha")),
      omega_(iConfig.getParameter<double>("omega")) {
}

void GHSJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  fillDescriptionsFromFlavouredJetProducerBase(desc);
  
  desc.add<double>("ptcut", 15.0);
  desc.add<double>("alpha", 1.0);
  desc.add<double>("omega", 2.0);
  
  desc.add<std::string>("jetCollInstanceName", "");
  descriptions.add("GHSJetProducer", desc);
}

void GHSJetProducer::runAlgorithm(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  FastjetJetProducer::runAlgorithm(iEvent, iSetup);
  
  fjBaseJets_ = fjJets_;

  fjJets_ = fastjet::contrib::run_GHS(fjBaseJets_, ptcut_, alpha_, omega_, *fjRecombiner_);
}

DEFINE_FWK_MODULE(GHSJetProducer);
