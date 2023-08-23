#ifndef HeterogeneousCore_AlpakaCore_interface_atomicMaxPair_h
#define HeterogeneousCore_AlpakaCore_interface_atomicMaxPair_h
#include <alpaka/alpaka.hpp>

#include "FWCore/Utilities/interface/bit_cast.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#if defined(__CUDA_ARCH__) or defined(__HIP_DEVICE_COMPILE__)
template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>, typename F>
ALPAKA_FN_ACC ALPAKA_FN_INLINE void atomicMaxPair(const TAcc& acc,
                                                  unsigned long long int* address,
                                                  std::pair<unsigned int, float> value,
                                                  F comparator) {
  unsigned long long int val = (static_cast<unsigned long long int>(value.first) << 32) + __float_as_uint(value.second);
  unsigned long long int expected;
  unsigned long long int old = *address;
  do {
    expected = old;
    if (comparator(std::pair<unsigned int, float>{static_cast<unsigned int>(val >> 32 & 0xffffffff),
                                                  __uint_as_float(val & 0xffffffff)},
                   value)) {
      old = atomicCAS(address, expected, val);
    } else {
      break;
    }
  } while (expected != old);
}

auto comparator = [](std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) {
  return a.second < b.second;
};
#else
template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>, typename F>
ALPAKA_FN_ACC ALPAKA_FN_INLINE void atomicMaxPair(const TAcc& acc,
                                                  unsigned long long int* address,
                                                  std::pair<unsigned int, float> value,
                                                  F comparator) {
  unsigned long long int val =
      (static_cast<unsigned long long int>(value.first) << 32) + edm::bit_cast<unsigned int>(value.second);
  unsigned long long int expected;
  unsigned long long int old = *address;
  do {
    expected = old;
    if (comparator(std::pair{static_cast<unsigned int>(val >> 32 & 0xffffffff),
                             edm::bit_cast<float>(static_cast<unsigned int>(val & 0xffffffff))},
                   value)) {
      old = alpaka::atomicCas(acc, address, expected, val);
    } else {
      break;
    }
  } while (expected != old);
}

auto comparator = [](std::pair<unsigned int, float> a, std::pair<unsigned int, float> b) {
  return a.second < b.second;
};
#endif  // __CUDA_ARCH__ or __HIP_DEVICE_COMPILE__

#endif  // HeterogeneousCore_AlpakaCore_interface_atomicMaxPair_h
