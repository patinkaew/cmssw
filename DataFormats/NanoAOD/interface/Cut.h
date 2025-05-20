#ifndef DataFormats_NanoAOD_Cut_h
#define DataFormats_NanoAOD_Cut_h

namespace nanoaod {
  namespace cut{

    template <typename T>
    struct Value {
      T v;
      Value () {}
      Value(T value) : v(value) {}
    };
    
    struct CutResults {
      std::vector<uint8_t> c;
      CutResults() {}
      CutResults(int size) : c(std::vector<uint8_t>(size)) {}
    };

  } // namespace cut
} // namespace nanoaod

#endif
