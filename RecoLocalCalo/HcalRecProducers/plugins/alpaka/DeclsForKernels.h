#ifndef RecoLocalCalo_HcalRecProducers_plugins_alpaka_DeclsForKernels_h
#define RecoLocalCalo_HcalRecProducers_plugins_alpaka_DeclsForKernels_h

#include <functional>
#include <optional>

#include "CondFormats/HcalObjects/interface/HcalChannelStatus.h"
#include "DataFormats/HcalDigi/interface/alpaka/HcalDigiDeviceCollection.h"
#include "DataFormats/HcalRecHit/interface/alpaka/HcalRecHitDeviceCollection.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
//#include "CondFormats/DataRecord/interface/HcalCombinedRecordsGPU.h"
//#include "CondFormats/DataRecord/interface/HcalGainWidthsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalGainsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalLUTCorrsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalQIEDataRcd.h"
//#include "CondFormats/DataRecord/interface/HcalQIETypesRcd.h"
//#include "CondFormats/DataRecord/interface/HcalRecoParamsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalSiPMCharacteristicsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalSiPMParametersRcd.h"
//#include "CondFormats/DataRecord/interface/HcalTimeCorrsRcd.h"
//#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

//#include "CondFormats/HcalObjects/interface/HcalConvertedEffectivePedestalWidthsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalConvertedEffectivePedestalsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalGainWidthsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalGainsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalLUTCorrsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalQIECodersGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalQIETypesGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalRecoParamsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalRespCorrsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalSiPMCharacteristicsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalSiPMParametersGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalTimeCorrsGPU.h"
//#include "CondFormats/HcalObjects/interface/HcalChannelQualityGPU.h"
#include "CondFormats/HcalObjects/interface/alpaka/HcalMahiConditionsDevice.h"
#include "CondFormats/HcalObjects/interface/alpaka/HcalSiPMCharacteristicsDevice.h"
#include "CondFormats/HcalObjects/interface/alpaka/HcalRecoParamWithPulseShapeDevice.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace hcal {
    namespace reconstruction {

      struct ConditionsProducts {
        HcalMahiConditionsPortableDevice const& mahi;
        HcalSiPMCharacteristicsPortableDevice const& sipmCharacteristics;
        HcalRecoParamWithPulseShapeDevice const& recoParams;
      };

      struct ConfigParameters {
        uint32_t maxTimeSamples;
        uint32_t kprep1dChannelsPerBlock;
        int sipmQTSShift;
        int sipmQNTStoSum;
        int firstSampleShift;
        bool useEffectivePedestals;

        float meanTime;
        float timeSigmaSiPM, timeSigmaHPD;
        float ts4Thresh;

        std::array<uint32_t, 3> kernelMinimizeThreads;

        // FIXME:
        //   - add "getters" to HcalTimeSlew calib formats
        //   - add ES Producer to consume what is produced above not to replicate.
        //   which ones to use is hardcoded, therefore no need to send those to the device
        bool applyTimeSlew;
        float tzeroTimeSlew, slopeTimeSlew, tmaxTimeSlew;
      };

    }  // namespace reconstruction
  }    // namespace hcal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_HcalRecProducers_plugins_DeclsForKernels_h
