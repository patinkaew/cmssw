#ifndef RecoJets_JetProducers_plugins_FlavouredJetProducer_h
#define RecoJets_JetProducers_plugins_FlavouredJetProducer_h

#include "RecoJets/JetProducers/interface/JetSpecific.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetMatching/interface/JetFlavourAlgoInfoMatching.h"
#include "RecoJets/JetProducers/plugins/FastjetJetProducer.h"

#include "fastjet/contrib/FlavInfo.hh"
#include "fastjet/contrib/SDFPlugin.hh"
#include "fastjet/contrib/CMPPlugin.hh"
#include "fastjet/contrib/IFNPlugin.hh"
#include "fastjet/contrib/GHSAlgo.hh"

class FlavouredJetProducerBase : public FastjetJetProducer {
public:
  explicit FlavouredJetProducerBase(const edm::ParameterSet& iConfig);
  ~FlavouredJetProducerBase() override = default;

  static void fillDescriptionsFromFlavouredJetProducerBase(edm::ParameterSetDescription& descriptions);

protected:
  void inputTowers() override; 
  void output(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  //void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  std::string recombination_;
  std::shared_ptr<fastjet::contrib::FlavRecombiner> fjRecombiner_;
};

////////////////////////////
// SDF (SoftDrop Flavour) //
////////////////////////////
class SDFJetProducer : public FlavouredJetProducerBase {
public:
  explicit SDFJetProducer(const edm::ParameterSet& iConfig);
  ~SDFJetProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void runAlgorithm(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

private:
  fastjet::contrib::SDFlavourCalc sdFlavCalc_;
};

/////////////////////////////
// CMP (Flavoured anti-kt) //
/////////////////////////////
class CMPJetProducer : public FlavouredJetProducerBase {
public:
  explicit CMPJetProducer(const edm::ParameterSet& iConfig);
  ~CMPJetProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  double a_;
  std::string correction_type_;
  std::string clustering_type_;
  std::shared_ptr<fastjet::CMPPlugin> fjCMPPlugin_;
};

/////////////////////////////////////////////
// IFN (Interleaved Flavour Neutralisation //
/////////////////////////////////////////////
class IFNJetProducer : public FlavouredJetProducerBase {
public: 
  explicit IFNJetProducer(const edm::ParameterSet& iConfig);
  ~IFNJetProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  double alpha_;
  double omega_;
  std::shared_ptr<fastjet::JetDefinition> fjBaseJetDefinition_;
  std::shared_ptr<fastjet::contrib::IFNPlugin> fjIFNPlugin_;
};

////////////////////////////
// GHS (Flavour Dressing) //
////////////////////////////
class GHSJetProducer : public FlavouredJetProducerBase {
public:
  explicit GHSJetProducer(const edm::ParameterSet& iConfig);
  ~GHSJetProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void runAlgorithm(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

private:
  double ptcut_;
  double alpha_;
  double omega_;
  std::vector<fastjet::PseudoJet> fjBaseJets_;
};

#endif // RecoJets_JetProducers_plugins_FlavouredJetProducer
