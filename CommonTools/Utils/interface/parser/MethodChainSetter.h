#ifndef CommonTools_Utils_MethodChainSetter_h
#define CommonTools_Utils_MethodChainSetter_h
/* \class reco::parser::MethodChainSetter
 *
 * Method Chain setter
 *
 *
 *
 */
#include "CommonTools/Utils/interface/parser/MethodChain.h"
#include "CommonTools/Utils/interface/parser/MethodStack.h"
#include "CommonTools/Utils/interface/parser/TypeStack.h"

#include <iostream>

namespace reco {
  namespace parser {
    struct MethodChainSetter {
      MethodChainSetter(MethodChainPtr &methchain,
                        MethodStack &methStack,
                        LazyMethodStack &lazyMethStack,
                        TypeStack &typeStack)
          : methchain_(methchain), methStack_(methStack), lazyMethStack_(lazyMethStack), typeStack_(typeStack) {
        std::cout << "MethodChainSetter from constructor" << std::endl;
        std::cout << methchain << ", " << methchain_ << std::endl;
        std::cout << &methStack << ", " << &methStack_ << std::endl;
        std::cout << &typeStack << ", " << &typeStack_ << std::endl;
        }
      void operator()(const char *, const char *) const;

    private:
      void push(const char *, const char *) const;
      //void lazyPush(const char *, const char *) const;

      MethodChainPtr &methchain_;
      MethodStack &methStack_;
      LazyMethodStack &lazyMethStack_;
      TypeStack &typeStack_;
    };
  }  // namespace parser
}  // namespace reco

#endif
