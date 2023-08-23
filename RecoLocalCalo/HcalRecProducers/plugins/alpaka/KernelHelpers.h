#ifndef RecoLocalCalo_HcalRecProducers_plugins_alpaka_KernelHelpers_h
#define RecoLocalCalo_HcalRecProducers_plugins_alpaka_KernelHelpers_h

#include <alpaka/alpaka.hpp>
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConstants.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::hcal::reconstruction {

  ALPAKA_FN_ACC ALPAKA_FN_INLINE float compute_time_slew_delay(float const fC,
                                                               float const tzero,
                                                               float const slope,
                                                               float const tmax);

  // this is from HcalTimeSlew.
  // HcalTimeSlew are values that come in from ESProducer that takes them
  // from a python config. see DeclsForKernels for more explanation
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float compute_time_slew_delay(float const fC,
                                                               float const tzero,
                                                               float const slope,
                                                               float const tmax);

  // Conditions are transferred once per IOV
  // Access is performed based on the det id which is converted to a linear index
  // 2 funcs below are taken from HcalTopology (reimplemented here).
  // Inputs are constants that are also taken from HcalTopology
  // but passed to the kernel as arguments using the HclaTopology itself

  ALPAKA_FN_ACC uint32_t did2linearIndexHB(
      uint32_t const didraw, int const maxDepthHB, int const firstHBRing, int const lastHBRing, int const nEtaHB);

  ALPAKA_FN_ACC uint32_t did2linearIndexHE(uint32_t const didraw,
                                           int const maxDepthHE,
                                           int const maxPhiHE,
                                           int const firstHERing,
                                           int const lastHERing,
                                           int const nEtaHE);

  ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t get_qiecoder_index(uint32_t const capid, uint32_t const range);

  ALPAKA_FN_ACC float compute_reco_correction_factor(float const par1,
                                                     float const par2,
                                                     float const par3,
                                                     float const x);

  // compute the charge using the adc, qie type and the appropriate qie shape array
  ALPAKA_FN_ACC float compute_coder_charge(
      int const qieType, uint8_t const adc, uint8_t const capid, float const* qieOffsets, float const* qieSlopes);

  ALPAKA_FN_ACC float compute_diff_charge_gain(int const qieType,
                                               uint8_t adc,
                                               uint8_t const capid,
                                               float const* qieOffsets,
                                               float const* qieSlopes,
                                               bool const isqie11);

  // TODO: remove what's not needed
  // originally from from RecoLocalCalo/HcalRecAlgos/src/PulseShapeFunctor.cc
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float compute_pulse_shape_value(float const pulse_time,
                                                                 int const sample,
                                                                 int const shift,
                                                                 float const* acc25nsVec,
                                                                 float const* diff25nsItvlVec,
                                                                 float const* accVarLenIdxMinusOneVec,
                                                                 float const* diffVarItvlIdxMinusOneVec,
                                                                 float const* accVarLenIdxZeroVec,
                                                                 float const* diffVarItvlIdxZeroVec);
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::hcal::reconstruction

#endif  // RecoLocalCalo_HcalRecProducers_plugins_alpaka_KernelHelpers_h
