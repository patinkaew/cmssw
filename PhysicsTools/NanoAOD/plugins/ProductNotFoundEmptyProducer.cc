// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


template<typename T>
class ProductNotFoundEmptyProducer: public edm::stream::EDProducer<> {
public:
    explicit ProductNotFoundEmptyProducer(const edm::ParameterSet&);
    ~ProductNotFoundEmptyProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    void endStream() override {}

private:
     const edm::EDGetTokenT<T> input_token_;
};

template<typename T>
ProductNotFoundEmptyProducer<T>::ProductNotFoundEmptyProducer(const edm::ParameterSet& iConfig)
  :input_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))){
    produces<T>();
}

template<typename T>
ProductNotFoundEmptyProducer<T>::~ProductNotFoundEmptyProducer() = default;

// ------------ method called to produce the data  ------------
template<typename T>
void ProductNotFoundEmptyProducer<T>::produce(edm::Event &iEvent, edm::EventSetup const &setup){
    try{
        edm::Handle<T> input_handle;
        iEvent.getByToken(input_token_, input_handle);
        //auto out = std::make_unique<T> (new T(iEvent.get(input_token_)));
        //std::unique_ptr<T> out;
        //out = std::move();
        iEvent.put(std::move(std::make_unique<T>(std::move(*input_handle))));
    } catch (const cms::Exception& ex) {
      if (ex.category() == "ProductNotFound") {
        auto empty = std::make_unique<T> ();
        iEvent.put(std::move(empty));
      }else{
        throw;
      }
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<typename T>
void ProductNotFoundEmptyProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
typedef ProductNotFoundEmptyProducer<double> DoubleProductNotFoundEmptyProducer;
DEFINE_FWK_MODULE(DoubleProductNotFoundEmptyProducer);

#include "DataFormats/JetReco/interface/PFJetCollection.h"
typedef ProductNotFoundEmptyProducer<reco::PFJetCollection> PFJetCollectionProductNotFoundEmptyProducer;
DEFINE_FWK_MODULE(PFJetCollectionProductNotFoundEmptyProducer);
