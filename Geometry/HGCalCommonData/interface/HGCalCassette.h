#ifndef Geometry_HGCalCommonData_HGCalCassette_h
#define Geometry_HGCalCommonData_HGCalCassette_h

#include <cmath>
#include <cstdint>
#include <vector>

class HGCalCassette {
public:
  HGCalCassette(int cassette, const std::vector<double>& shifts) { setParameter(cassette, shifts, true); }
  HGCalCassette() {}

  void setParameter(int cassette, const std::vector<double>& shifts, bool both = true);
  void setParameterScint(const std::vector<double>& shifts);
  void setParameterRetract(const std::vector<double>& shifts);
  std::pair<double, double> getShift(int layer, int zside, int cassette, bool scnt = false) const;
  std::pair<double, double> getShiftScnt(int layer, int zside, double phi) const;
  static int cassetteIndex(int det, int layer, int zside, int cassette);
  static int cassetteType(int det, int zside, int cassette);

private:
  static constexpr int positEE_[6] = {2, 1, 0, 5, 4, 3};
  static constexpr int positHE_[12] = {5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  int cassette_;
  bool typeHE_;
  std::vector<double> shifts_, shiftsScnt_, retractScnt_;
  static constexpr int32_t factor_ = 100;
};

#endif
