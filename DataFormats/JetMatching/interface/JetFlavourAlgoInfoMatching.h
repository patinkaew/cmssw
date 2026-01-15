#ifndef DataFormats_JetMatching_JetFlavourAlgoInfoMatching_h
#define DataFormats_JetMatching_JetFlavourAlgoInfoMatching_h

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/RefToBaseProd.h"
#include "DataFormats/JetMatching/interface/JetFlavourAlgoInfo.h"
#include <vector>

namespace reco {

  typedef edm::AssociationVector<edm::RefToBaseProd<reco::Jet>, std::vector<reco::JetFlavourAlgoInfo> >
      JetFlavourAlgoInfoMatchingCollectionBase;

  class JetFlavourAlgoInfoMatchingCollection : public JetFlavourAlgoInfoMatchingCollectionBase {
  public:
    JetFlavourAlgoInfoMatchingCollection() : JetFlavourAlgoInfoMatchingCollectionBase() {}

    JetFlavourAlgoInfoMatchingCollection(const JetFlavourAlgoInfoMatchingCollectionBase &v)
        : JetFlavourAlgoInfoMatchingCollectionBase(v) {}
  };

  typedef JetFlavourAlgoInfoMatchingCollection::value_type JetFlavourAlgoInfoMatching;

  typedef edm::Ref<JetFlavourAlgoInfoMatchingCollection> JetFlavourAlgoInfoMatchingRef;

  typedef edm::RefProd<JetFlavourAlgoInfoMatchingCollection> JetFlavourAlgoInfoMatchingRefProd;

  typedef edm::RefVector<JetFlavourAlgoInfoMatchingCollection> JetFlavourAlgoInfoMatchingRefVector;

}  // namespace reco

#endif  // DataFormats_JetMatching_JetFlavourAlgoInfoMatching_h
