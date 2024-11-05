#ifndef CommonTools_Utils_TypedStringObjectMethodCaller_h
#define CommonTools_Utils_TypedStringObjectMethodCaller_h
/* \class TypedStringObjectMethodCaller
 *
 * \author Patin Inkaew, HIP
 *
 */

#include "boost/spirit/include/classic_core.hpp"
#include "boost/spirit/include/classic_grammar_def.hpp"
#include "boost/spirit/include/classic_chset.hpp"

#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/Utils/interface/parser/Exception.h"
#include "CommonTools/Utils/interface/parser/MethodChain.h"
#include "CommonTools/Utils/interface/parser/MethodChainGrammar.h"
#include "FWCore/Reflection/interface/ObjectWithDict.h"
#include "CommonTools/Utils/interface/parser/MethodSetter.h"
#include "CommonTools/Utils/interface/parser/MethodInvoker.h"
#include "CommonTools/Utils/interface/parser/MethodStack.h"
#include "CommonTools/Utils/interface/parser/TypeStack.h"
#include "CommonTools/Utils/interface/parser/MethodArgumentStack.h"
#include "CommonTools/Utils/interface/parser/AnyMethodArgument.h"

/*
template <typename T, typename R>
struct TypedStringObjectMethodCaller {
  TypedStringObjectMethodCaller(const std::string &name)
    : type_(typeid(T)), singleInvoker_(type_, name, std::vector<reco::parser::AnyMethodArgument>()) {
  }
  R operator()(const T &t) const {
    edm::ObjectWithDict o(type_, const_cast<T *>(&t));
    std::vector<reco::parser::StorageManager> storage;
    storage.reserve(1);
    std::pair<edm::ObjectWithDict, bool> ret = singleInvoker_.invoke(o, storage);
    void* addr = ret.first.address();
    return *static_cast<R*>(addr);
  }

private:
  edm::TypeWithDict type_;
  reco::parser::SingleInvoker singleInvoker_;
}; */

template <typename T, typename R>
struct TypedStringObjectMethodCaller {
  TypedStringObjectMethodCaller(const std::string expr) : type_(typeid(T)) {
    using namespace boost::spirit::classic;
    reco::parser::MethodChainGrammar grammar(methodchain_, type_, false);
    const char* startingFrom = expr.c_str();
    //parse(startingFrom, grammar.use_parser<1>() >> end_p, space_p).full;
    std::cout << "MethodChainGrammar starts parsing " << expr << std::endl;
    std::cout << &methodchain_ << std::endl;
    try {
      parse(startingFrom, grammar >> end_p, space_p).full;
      //parse(startingFrom, real_p >> *(',' >> real_p), space_p).full;
      std::cout << "MethodChainGrammer parse successfully" << std::endl;
    } catch (reco::parser::BaseException& e) {
      throw edm::Exception(edm::errors::Configuration)
        << "MethodChainGrammer parse error:" << reco::parser::baseExceptionWhat(e) << " (char " << e.where - startingFrom << ")\n";
    }
  }

  R operator()(const T &t) const {
    edm::ObjectWithDict o(type_, const_cast<T *>(&t));
    edm::ObjectWithDict ret = methodchain_->value(o);
    return *static_cast<R*>(ret.address());
  }

private:
  reco::parser::MethodChainPtr methodchain_;
  edm::TypeWithDict type_;
};

#endif
