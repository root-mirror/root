#ifndef ROOT_REP
#define ROOT_REP

#include <string>

template <int D1, int D2>
struct RepStd {
   typedef typename ROOT::Math::MatRepStd<double, D1, D2> R_t;
   static std::string Name()
   {
      return "ROOT::Math::MatRepStd<double," + ROOT::Math::Util::ToString(D1) + "," + ROOT::Math::Util::ToString(D2) +
             "> ";
   }
   static std::string Name32()
   {
      return "ROOT::Math::MatRepStd<Double32_t," + ROOT::Math::Util::ToString(D1) + "," +
             ROOT::Math::Util::ToString(D2) + "> ";
   }
   static std::string SName() { return ""; }
};
template <int D1>
struct RepSym {
   typedef typename ROOT::Math::MatRepSym<double, D1> R_t;
   static std::string Name() { return "ROOT::Math::MatRepSym<double," + ROOT::Math::Util::ToString(D1) + "> "; }
   static std::string Name32() { return "ROOT::Math::MatRepSym<Double32_t," + ROOT::Math::Util::ToString(D1) + "> "; }
   static std::string SName() { return "MatRepSym"; }
};

#endif // ROOT_REP
