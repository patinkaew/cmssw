#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiHostCollection.h"
#include "DataFormats/HcalDigi/interface/alpaka/HcalDigiDeviceCollection.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    class HcalDigisProducerPortable : public stream::EDProducer<> {
    public:
      explicit HcalDigisProducerPortable(edm::ParameterSet const& ps);
      ~HcalDigisProducerPortable() override = default;
      static void fillDescriptions(edm::ConfigurationDescriptions&);
    
    private:
      void produce(device::Event&, device::EventSetup const&) override;
    
    private:
      // input product tokens
      edm::EDGetTokenT<HBHEDigiCollection> hbheDigiToken_;
      edm::EDGetTokenT<QIE11DigiCollection> qie11DigiToken_;
    
      // type aliases
      using HostCollectionPhase1 = hcal::Phase1DigiHostCollection;
      using HostCollectionPhase0 = hcal::Phase0DigiHostCollection;

      using DeviceCollectionPhase1 = hcal::Phase1DigiDeviceCollection;
      using DeviceCollectionPhase0 = hcal::Phase0DigiDeviceCollection;
    
      // output product tokens
      device::EDPutToken<DeviceCollectionPhase1> digisF01HEToken_;
      device::EDPutToken<DeviceCollectionPhase0> digisF5HBToken_;
      device::EDPutToken<DeviceCollectionPhase1> digisF3HBToken_;
    
    
      struct ConfigParameters {
        uint32_t maxChannelsF01HE, maxChannelsF5HB, maxChannelsF3HB;
      };
      ConfigParameters config_;
    
      // per event host buffers
      //HostCollectionPhase1 hf01_;
      //HostCollectionPhase1 hf3_;
      //HostCollectionPhase0 hf5_;
    
      //// device products: product owns memory (i.e. not the module)
      //DeviceCollectionPhase1 df01_;
      //DeviceCollectionPhase1 df3_;
      //DeviceCollectionPhase0 df5_;
    };
    
    void HcalDigisProducerPortable::fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
      edm::ParameterSetDescription desc;
    
      // FIXME
      desc.add<edm::InputTag>("hbheDigisLabel", edm::InputTag("hcalDigis"));
      desc.add<edm::InputTag>("qie11DigiLabel", edm::InputTag("hcalDigis"));
      desc.add<std::string>("digisLabelF01HE", std::string{"f01HEDigisGPU"});
      desc.add<std::string>("digisLabelF5HB", std::string{"f5HBDigisGPU"});
      desc.add<std::string>("digisLabelF3HB", std::string{"f3HBDigisGPU"});
      desc.add<uint32_t>("maxChannelsF01HE", 10000u);
      desc.add<uint32_t>("maxChannelsF5HB", 10000u);
      desc.add<uint32_t>("maxChannelsF3HB", 10000u);
    
      confDesc.addWithDefaultLabel(desc);
    }
    
    HcalDigisProducerPortable::HcalDigisProducerPortable(const edm::ParameterSet& ps)
        : hbheDigiToken_{consumes<HBHEDigiCollection>(ps.getParameter<edm::InputTag>("hbheDigisLabel"))},
          qie11DigiToken_{consumes<QIE11DigiCollection>(ps.getParameter<edm::InputTag>("qie11DigiLabel"))},
          digisF01HEToken_{produces(ps.getParameter<std::string>("digisLabelF01HE"))},
          digisF5HBToken_{produces(ps.getParameter<std::string>("digisLabelF5HB"))},
          digisF3HBToken_{produces(ps.getParameter<std::string>("digisLabelF3HB"))} {
      config_.maxChannelsF01HE = ps.getParameter<uint32_t>("maxChannelsF01HE");
      config_.maxChannelsF5HB = ps.getParameter<uint32_t>("maxChannelsF5HB");
      config_.maxChannelsF3HB = ps.getParameter<uint32_t>("maxChannelsF3HB");
    
      // this is a preallocation for the max statically known number of time samples
      // actual stride/nsamples will be inferred from data
      //hf01_.stride = hcal::compute_stride<hcal::Flavor1>(QIE11DigiCollection::MAXSAMPLES);
      //hf5_.stride = hcal::compute_stride<hcal::Flavor5>(HBHEDataFrame::MAXSAMPLES);
      //hf3_.stride = hcal::compute_stride<hcal::Flavor3>(QIE11DigiCollection::MAXSAMPLES);
    
    }
    
    void HcalDigisProducerPortable::produce(device::Event& event, device::EventSetup const& setup) {
      
      const auto hbheDigis = event.getHandle(hbheDigiToken_);
      const auto qie11Digis = event.getHandle(qie11DigiToken_);
   
      auto const nsamples = (*hbheDigis)[0].size();
      //stride = nsamples * WORDS_PER_SAMPLE + Flavor::HEADER_WORDS; 
      //TODO:: get HEADER_WORDS/WORDS_PER_SAMPLE from DataFormat
      auto const stride = nsamples * 0.5 + 1;
      auto const size = hbheDigis->size() * stride;           // number of channels * stride

      // stack host memory in the queue
      HostCollectionPhase0 hf5_(size,event.queue());
      // device product
      auto df5_ = std::make_unique<DeviceCollectionPhase0>(size, alpaka::getDev(event.queue()));
      // set SoA_Scalar;
      hf5_.view().stride() = stride;

      for (unsigned int i =0 ; i< hbheDigis->size();++i){
        auto const& hbhe = (*hbheDigis)[i];
        auto const id = hbhe.id().rawId();
        auto const presamples = hbhe.presamples();
        uint16_t   header_word = (1 << 15) | (0x5 << 12) | (0 << 10) | ((hbhe.sample(0).capid() & 0x3) << 8);

        auto hf5_vi          = hf5_.view()[i];  
        hf5_vi.ids()         = id;
        hf5_vi.npresamples() = presamples;
        hf5_vi.data()[0]     = header_word;
        //TODO:: get HEADER_WORDS/WORDS_PER_SAMPLE from DataFormat
        for (unsigned int i = 0; i < hf5_.view().stride() - 1 ; i++) {
          uint16_t s0 = (0 << 7) | (static_cast<uint8_t>(hbhe.sample(2 * i).adc()) & 0x7f);
          uint16_t s1 = (0 << 7) | (static_cast<uint8_t>(hbhe.sample(2 * i + 1).adc()) & 0x7f);
          uint16_t sample = (s1 << 8) | s0;
          hf5_vi.data()[i+1] = sample;
        }
      }
      // copy hf5 to df5
      alpaka::memcpy(event.queue(), (*df5_).buffer(), hf5_.const_buffer());
      
      // put to event 
      event.put(digisF5HBToken_, std::move(df5_)); 

    }
}
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HcalDigisProducerPortable);
