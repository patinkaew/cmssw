#ifndef CommonTools_Utils_MethodChain_h
#define CommonTools_Utils_MethodChain_h

/* \class reco::parser::MethodChain
 *
 * Chain of methods
 *
 * \author original version: Patin Inkaew, HIP
 *
 */

#include "CommonTools/Utils/interface/parser/MethodInvoker.h"
#include "CommonTools/Utils/interface/TypeCode.h"

#include <vector>
#include <oneapi/tbb/concurrent_queue.h>

namespace reco {
  namespace parser {

    /// Evaluate an object's method or datamember (or chain of them)
    class MethodChain {
    private:  // Private Data Members
      std::vector<MethodInvoker> methods_;
      using Objects = std::vector<std::pair<edm::ObjectWithDict, bool>>;
      mutable oneapi::tbb::concurrent_queue<Objects> objectsCache_;

    private:  // Private Methods
      [[nodiscard]] Objects initObjects_() const;

      Objects borrowObjects() const;
      void returnObjects(Objects&&) const;

    public:  // Public Static Methods

      /// allocate an object to hold the result of a given member (if needed)
      /// this method is used also from the LazyInvoker code
      /// returns true if objects returned from this will require a destructor
      static bool makeStorage(edm::ObjectWithDict& obj, const edm::TypeWithDict& retType);

      /// delete an objecty, if needed
      /// this method is used also from the LazyInvoker code
      static void delStorage(edm::ObjectWithDict&);

    public:  // Public Methods
      MethodChain(const std::vector<MethodInvoker>& methods);
      MethodChain(const MethodChain&);
      ~MethodChain();
      edm::ObjectWithDict value(const edm::ObjectWithDict&) const;
    };
/*
    /// Same as ExpressionVar but with lazy resolution of object methods
    /// using the dynamic type of the object, and not the one fixed at compile time
    class ExpressionLazyVar : public ExpressionBase {
    private:  // Private Data Members
      std::vector<LazyInvoker> methods_;

    public:
      ExpressionLazyVar(const std::vector<LazyInvoker>& methods);
      ~ExpressionLazyVar() override;
      double value(const edm::ObjectWithDict&) const override;
    };
*/
    typedef std::shared_ptr<MethodChain> MethodChainPtr;
    //typedef std::vector<MethodChain> MethodChainStack;
  }  // namespace parser
}  // namespace reco

#endif  // CommonTools_Utils_MethodChain_h
