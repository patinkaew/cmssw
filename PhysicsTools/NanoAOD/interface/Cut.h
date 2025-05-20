#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/NanoAOD/interface/Cut.h"

template <typename ObjType>
class CutBase {
public:
  virtual ~CutBase() = default;
  virtual void beginEvent(const edm::Event &iEvent) = 0;
  virtual bool operator()(const edm::Ptr<ObjType> &) = 0;
};

template <typename ObjType>
class Cut : public CutBase<ObjType> {
public:
  Cut(const edm::ParameterSet &cfg)
      : cut_(cfg.getParameter<std::string>("cut"),
             cfg.existsAs<bool>("lazyEval") ? cfg.getParameter<bool>("lazyEval") : false) {}
  void beginEvent(const edm::Event &iEvent) {}
  bool operator()(const edm::Ptr<ObjType> &obj) {
    return cut_(*obj);
  }

private:
  StringCutObjectSelector<ObjType> cut_;
};

template <typename ObjType, typename ValType>
class ValueMapCut : public CutBase<ObjType> {
public:
  ValueMapCut(const edm::ParameterSet &cfg, edm::ConsumesCollector &cc)
      : cut_(cfg.getParameter<std::string>("cut"),
             cfg.existsAs<bool>("lazyEval") ? cfg.getParameter<bool>("lazyEval") : false),
        token_(cc.consumes<edm::ValueMap<ValType>>(cfg.getParameter<edm::InputTag>("src"))),
        value_(nanoaod::cut::Value<ValType>()) {}

  void beginEvent(const edm::Event &iEvent) {
    vmap_ = iEvent.getHandle(token_);
  }
  
  bool operator()(const edm::Ptr<ObjType> &obj) {
    value_.v = (*vmap_)[obj]; // set value
    return cut_(value_);
  }

private:
  StringCutObjectSelector<nanoaod::cut::Value<ValType>> cut_;
  edm::EDGetTokenT<edm::ValueMap<ValType>> token_;
  nanoaod::cut::Value<ValType> value_;
  edm::Handle<edm::ValueMap<ValType>> vmap_;
};

template <typename ObjType>
class Cuts : public CutBase<ObjType> {
public:
  Cuts(const edm::ParameterSet &cfg, edm::ConsumesCollector &&cc) 
      : mode_(cfg.existsAs<std::string>("mode") && !cfg.existsAs<std::string>("expr") 
            ? cfg.getParameter<std::string>("mode") : ""),
        expr_(cfg.existsAs<std::string>("expr") && !cfg.existsAs<std::string>("mode")
            ? cfg.getParameter<std::string>("expr") : ""),
        cutResultsEvaluator_(expr_) {
    // default to mode = "ALL"
    if (mode_.empty() && expr_.empty()) {
      mode_ = "ALL";
    }
    // build cuts
    if (cfg.existsAs<std::vector<edm::ParameterSet>>("cuts")) {
      const auto &cutsPSet = cfg.getParameter<std::vector<edm::ParameterSet>>("cuts");
      for (const auto &cutPSet : cutsPSet) {
        if (!cutPSet.existsAs<edm::InputTag>("src")) {
          cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new Cut<ObjType>(cutPSet)));
        } else {
          const std::string &type = cutPSet.getParameter<std::string>("type");
          if (type == "int")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new IntValueMapCut(cutPSet, cc)));
          else if (type == "uint")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new UIntValueMapCut(cutPSet, cc)));
          else if (type == "int64")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new Int64ValueMapCut(cutPSet, cc)));
          else if (type == "uint64")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new UInt64ValueMapCut(cutPSet, cc)));
          else if (type == "float")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new FloatValueMapCut(cutPSet, cc)));
          else if (type == "double")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new DoubleValueMapCut(cutPSet, cc)));
          else if (type == "uint8")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new UInt8ValueMapCut(cutPSet, cc)));
          else if (type == "int16")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new Int16ValueMapCut(cutPSet, cc)));
          else if (type == "uint16")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new UInt16ValueMapCut(cutPSet, cc)));
          else if (type == "bool")
            cuts_.push_back(std::unique_ptr<CutBase<ObjType>>(new BoolValueMapCut(cutPSet, cc)));
          else
            throw cms::Exception("Configuration", "unsupported type " + type);
        }
      }
      if (!expr_.empty()) {
        // allocate memory for cut results
        cutResults_ = nanoaod::cut::CutResults(cuts_.size());
      }
    }
  }

  void beginEvent(const edm::Event& iEvent) {
    for (const auto &cut : cuts_)
      cut->beginEvent(iEvent);
  }

  static edm::ParameterSetDescription cutsDescription() {
    edm::ParameterSetDescription cut;
    cut.add<std::string>("cut", "")->setComment("a logical expression to define the cut");
    cut.addUntracked<bool>("lazyEval", false);
    cut.addOptionalNode(
        edm::ParameterDescription<edm::InputTag>("src", true,
            edm::Comment("input tag of the ValueMap used in the cut")) and 
        edm::ParameterSwitch<std::string>(
            edm::ParameterDescription<std::string>("type", "float", true),
            std::move(edm::allowedValues<std::string>(
                "int", "uint", "int64", "uint64", "float", "double", "uint8", "int16", "uint16", "bool"))),
        false);
    edm::ParameterSetDescription cuts;
    cuts.setComment("a parameters set to define cuts");
    cuts.addNode(
        edm::ParameterSwitch<std::string>(
            edm::ParameterDescription<std::string>("mode", "ALL", true,
                edm::Comment("a mode (ALL or ANY) to combine all cuts to evaluate a final decision")),
            std::move(edm::allowedValues<std::string>("ALL", "ANY"))) xor
        edm::ParameterDescription<std::string>("expr", "", true,
            edm::Comment("a logical statement to combine all cuts to evaluate a final decision")));
    cuts.addVPSetOptional("cuts", cut, {edm::ParameterSet()});

    return cuts;
  }
  
  bool operator()(const edm::Ptr<ObjType> &obj) {
    if (mode_ == "ALL") {
      for (const auto &cut : cuts_) {
        if (!(*cut)(obj)) return false;
      }
      return true;
    } else if (mode_ == "ANY") {
      for (const auto &cut : cuts_) {
        if ((*cut)(obj)) return true;
      }
      return false;
    } else { // use expr
      // evaluate each cut and fill in the results
      for (unsigned int i = 0; i < cuts_.size(); i++ ) {
        const auto &cut = cuts_[i];
        cutResults_.c[i] = (*cut)(obj);
      }
      return cutResultsEvaluator_(cutResults_);
    }
  }

private:
  std::string mode_;
  std::string expr_;
  std::vector<std::unique_ptr<CutBase<ObjType>>> cuts_;
  nanoaod::cut::CutResults cutResults_;
  StringCutObjectSelector<nanoaod::cut::CutResults> cutResultsEvaluator_;

  typedef ValueMapCut<ObjType, int32_t> IntValueMapCut;
  typedef ValueMapCut<ObjType, uint32_t> UIntValueMapCut;
  typedef ValueMapCut<ObjType, int64_t> Int64ValueMapCut;
  typedef ValueMapCut<ObjType, uint64_t> UInt64ValueMapCut;
  typedef ValueMapCut<ObjType, float> FloatValueMapCut;
  typedef ValueMapCut<ObjType, double> DoubleValueMapCut;
  typedef ValueMapCut<ObjType, bool> BoolValueMapCut;
  typedef ValueMapCut<ObjType, uint8_t> UInt8ValueMapCut;
  typedef ValueMapCut<ObjType, int16_t> Int16ValueMapCut;
  typedef ValueMapCut<ObjType, uint16_t> UInt16ValueMapCut;
};
