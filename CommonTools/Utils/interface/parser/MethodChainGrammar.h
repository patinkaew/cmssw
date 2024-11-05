#ifndef CommonTools_Utils_MethodChainGrammar_h
#define CommonTools_Utils_MethodChainGrammar_h
/* \class MethodChainGrammer
 *
 * subset grammar of the full Grammer (CommonTools/Utils/interface/parser/Grammar.h),
 * allowing only a chain of methods 
 *
 */

#include "boost/spirit/include/classic_core.hpp"
#include "boost/spirit/include/classic_grammar_def.hpp"
#include "boost/spirit/include/classic_chset.hpp"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Reflection/interface/ObjectWithDict.h"
#include "CommonTools/Utils/interface/parser/MethodChain.h"
#include "CommonTools/Utils/interface/parser/MethodChainSetter.h"
#include "CommonTools/Utils/interface/parser/MethodSetter.h"
#include "CommonTools/Utils/interface/parser/MethodArgumentSetter.h"
#include "CommonTools/Utils/interface/parser/MethodInvoker.h"
#include "CommonTools/Utils/interface/parser/MethodStack.h"
#include "CommonTools/Utils/interface/parser/TypeStack.h"
#include "CommonTools/Utils/interface/parser/MethodArgumentStack.h"
#include "CommonTools/Utils/interface/parser/AnyMethodArgument.h"
#include "CommonTools/Utils/interface/parser/Exception.h"

#include <iostream>

namespace reco {
  namespace parser {
    struct MethodChainGrammar : public boost::spirit::classic::grammar<MethodChainGrammar> {
      MethodChainPtr *methchain_;
      bool lazy_; 
      mutable MethodStack methStack;
      mutable LazyMethodStack lazyMethStack;
      mutable MethodArgumentStack methArgStack;
      mutable TypeStack typeStack;

      MethodChainGrammar(MethodChainPtr &methchain, const edm::TypeWithDict& iType, bool lazy = false)
          : methchain_(&methchain), lazy_(lazy) {
        std::cout << "calling from MethodChainGrammar constructor" << std::endl;
        std::cout << methchain << std::endl;
        std::cout << &methchain << std::endl;
        std::cout << methchain_ << std::endl;
        std::cout << iType.name() << std::endl;
        typeStack.push_back(iType);
      }

      template <typename ScannerT>
      //struct definition : public boost::spirit::classic::grammar_def<boost::spirit::classic::rule<ScannerT>,
      //                                                               boost::spirit::classic::same,
      //                                                               boost::spirit::classic::same> {
      struct definition {
        typedef boost::spirit::classic::rule<ScannerT> rule;
        rule metharg, method, arrayAccess, methodchain;
        definition(const MethodChainGrammar& self) {
          using namespace boost::spirit::classic;

          std::cout << "calling from MethodChainGrammar" << std::endl;
          std::cout << *self.methchain_ << std::endl;
          std::cout << &self.methStack << std::endl;

          MethodArgumentSetter methodArg_s(self.methArgStack);
          MethodSetter method_s(self.methStack, self.lazyMethStack, self.typeStack, self.methArgStack, self.lazy_);
          MethodChainSetter methodchain_s(*self.methchain_, self.methStack, self.lazyMethStack, self.typeStack);

          BOOST_SPIRIT_DEBUG_RULE(methodchain);
          BOOST_SPIRIT_DEBUG_RULE(arrayAccess);
          BOOST_SPIRIT_DEBUG_RULE(method);
          BOOST_SPIRIT_DEBUG_RULE(metharg);

          std::cout << "calling from MethodChainGrammar after debug rules" << std::endl;

          boost::spirit::classic::assertion<SyntaxErrors> expectParenthesis(kMissingClosingParenthesis);
          boost::spirit::classic::assertion<SyntaxErrors> expect(kSyntaxError);

          std::cout << "calling from MethodChainGrammar after assertion" << std::endl;

          metharg = (strict_real_p[methodArg_s]) | (int_p[methodArg_s]) |
                    (ch_p('"') >> *(~ch_p('"')) >> ch_p('"'))[methodArg_s] |
                    (ch_p('\'') >> *(~ch_p('\'')) >> ch_p('\''))[methodArg_s];
          method = // alnum_p doesn't accept underscores, so we use chset<>; lexeme_d needed to avoid whitespace skipping within method names
                   (lexeme_d[alpha_p >> *chset<>("a-zA-Z0-9_")] >> ch_p('(') >> metharg >> *(ch_p(',') >> metharg) >>
                    expectParenthesis(ch_p(')')))[method_s] |
                   ((lexeme_d[alpha_p >> *chset<>("a-zA-Z0-9_")])[method_s] >> !(ch_p('(') >> ch_p(')')));
          arrayAccess = (ch_p('[') >> metharg >> *(ch_p(',') >> metharg) >> expectParenthesis(ch_p(']')))[method_s];
          methodchain = (method >> *(arrayAccess | (ch_p('.') >> expect(method))))[methodchain_s];
          //methodchain = (method >> *((ch_p('.') >> expect(method))));
         
          //this->start_parsers(methodchain);

          std::cout << "calling from MethodChainGrammar after this->start_parsers" << std::endl;
        }

        rule const& start() const { return methodchain; }
      };
    };
  }
}

#endif
