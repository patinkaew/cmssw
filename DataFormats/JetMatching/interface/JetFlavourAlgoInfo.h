#ifndef DataFormats_JetMatching_JetFlavourAlgoInfo_H
#define DataFormats_JetMatching_JetFlavourAlgoInfo_H

#include <array>
#include <cmath>

namespace reco {
  /**\class JetFlavourAlgoInfo JetFlavourAlgoInfo.h DataFormats/JetMatching/interface/JetFlavourAlgoInfo.h
 * \brief Class storing the jet flavour information from flavoured jet algorithms
 */

  class JetFlavourAlgoInfo {
  public:
    typedef std::array<int, 7> FlavourContent;

    enum {i_flag = 0, i_down, i_up, i_strange, i_charm, i_bottom, i_top};

    JetFlavourAlgoInfo() 
      : m_flavourContent{0, 0, 0, 0, 0, 0, 0}, 
        m_flavourContentMod2(m_flavourContent) {}

    JetFlavourAlgoInfo(int* flavourContent) {
      m_flavourContent[0] = flavourContent[0];
      m_flavourContentMod2[0] = flavourContent[0];
      for (unsigned int i = 1; i <= 6; i++) {
        m_flavourContent[i] = flavourContent[i];
        m_flavourContentMod2[i] = std::abs(flavourContent[i]) % 2;
      }
    }

    JetFlavourAlgoInfo(FlavourContent flavourContent) : JetFlavourAlgoInfo(flavourContent.data()) {}

    FlavourContent flavourContent() { return m_flavourContent; }

    FlavourContent flavourContentMod2() { return m_flavourContentMod2; }

    int flag() const { return m_flavourContent[i_flag]; }

    int flavour(int index) const { return m_flavourContent[index]; }
    int d() const { return flavour(i_down); }
    int u() const { return flavour(i_up); }
    int s() const { return flavour(i_strange); }
    int c() const { return flavour(i_charm); }
    int b() const { return flavour(i_bottom); }
    int t() const { return flavour(i_top); }

    int flavourMod2(int index) const { return m_flavourContentMod2[index]; }
    int dMod2() const { return flavourMod2(i_down); }
    int uMod2() const { return flavourMod2(i_up); }
    int sMod2() const { return flavourMod2(i_strange); }
    int cMod2() const { return flavourMod2(i_charm); }
    int bMod2() const { return flavourMod2(i_bottom); }
    int tMod2() const { return flavourMod2(i_top); }

    int iHeaviestFlavour() const {
      for (unsigned int i = 6; i > 0; i--) {
        if (m_flavourContent[i] != 0) {
          return i;
        }
      }
      return 0;
    }

    int iHeaviestFlavourMod2() const {
      for (unsigned int i = 6; i > 0; i--) {
        if (m_flavourContentMod2[i] != 0) {
          return i;
        }
      }
      return 0;
    }

    int heaviestFlavour() const {
      auto h = iHeaviestFlavour();
      switch(h) {
        case 0: return 21;
        case 1: return d() > 0 ? 2 : -2;
        case 2: return u() > 0 ? 1 : -1;
        default: return m_flavourContent[h] > 0 ? h : -h;
      };
    }

    int heaviestFlavourMod2() const {
      auto h = iHeaviestFlavourMod2();
      switch(h) {
        case 0: return 21;
        case 1: return 2;
        case 2: return 1;
        default: return h;
      };
    }

  private:
    FlavourContent m_flavourContent;
    FlavourContent m_flavourContentMod2;
    std::string m_recombiner;
  };
  
  typedef std::vector<JetFlavourAlgoInfo> JetFlavourAlgoInfoCollection;
}  // namespace reco
#endif
