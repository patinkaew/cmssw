#include "CommonTools/Utils/interface/parser/MethodChainSetter.h"
#include "CommonTools/Utils/interface/parser/MethodChain.h"
#include "CommonTools/Utils/interface/returnType.h"
#include "CommonTools/Utils/interface/parser/Exception.h"
#include <string>
#include <iostream>

using namespace reco::parser;
using namespace std;

void MethodChainSetter::operator()(const char *begin, const char *end) const {
  std::cout << "MethodChainSetter: Pushed [" << std::string(begin,end) << "]" << std::endl;
  if (!methStack_.empty())
    push(begin, end);
  /*else if (!lazyMethStack_.empty())
    lazyPush(begin, end);*/
  else
    throw Exception(begin) << " Expression didn't parse neither hastily nor lazyly. This must not happen.\n";
}

void MethodChainSetter::push(const char *begin, const char *end) const {
  edm::TypeWithDict type = typeStack_.back();
  std::cout << "MethodChainSetter push start" << std::endl;
  methchain_ = std::make_shared<MethodChain>(methStack_);
  //methchain_ = std::shared_ptr<MethodChain>(new MethodChain(methStack_));
  std::cout << "MethodChainSetter push done" << std::endl;
  methStack_.clear();
  typeStack_.resize(1);
}

/*
// not support for now
void MethodChainSetter::lazyPush(const char *begin, const char *end) const {
  //exprStack_.push_back(std::shared_ptr<ExpressionBase>(new ExpressionLazyVar(lazyMethStack_)));
  lazyMethStack_.clear();
  typeStack_.resize(1);
}
*/
