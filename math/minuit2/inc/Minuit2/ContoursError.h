// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ContoursError
#define ROOT_Minuit2_ContoursError

#include "Minuit2/MnConfig.h"
#include "Minuit2/MinosError.h"

#include <vector>
#include <utility>

namespace ROOT {

namespace Minuit2 {

class ContoursError {

public:
   ContoursError(unsigned int parx, unsigned int pary, const std::vector<std::pair<double, double>> &points,
                 const MinosError &xmnos, const MinosError &ymnos, unsigned int nfcn)
      : fParX(parx), fParY(pary), fPoints(points), fXMinos(xmnos), fYMinos(ymnos), fNFcn(nfcn)
   {
   }

   ~ContoursError() {}

   ContoursError(const ContoursError &cont)
      : fParX(cont.fParX), fParY(cont.fParY), fPoints(cont.fPoints), fXMinos(cont.fXMinos), fYMinos(cont.fYMinos),
        fNFcn(cont.fNFcn)
   {
   }

   ContoursError &operator()(const ContoursError &cont)
   {
      fParX = cont.fParX;
      fParY = cont.fParY;
      fPoints = cont.fPoints;
      fXMinos = cont.fXMinos;
      fYMinos = cont.fYMinos;
      fNFcn = cont.fNFcn;
      return *this;
   }

   const std::vector<std::pair<double, double>> &operator()() const { return fPoints; }

   std::pair<double, double> XMinos() const { return fXMinos(); }

   std::pair<double, double> YMinos() const { return fYMinos(); }

   unsigned int Xpar() const { return fParX; }
   unsigned int Ypar() const { return fParY; }

   const MinosError &XMinosError() const { return fXMinos; }

   const MinosError &YMinosError() const { return fYMinos; }

   unsigned int NFcn() const { return fNFcn; }
   double XMin() const { return fXMinos.Min(); }
   double YMin() const { return fYMinos.Min(); }

private:
   unsigned int fParX;
   unsigned int fParY;
   std::vector<std::pair<double, double>> fPoints;
   MinosError fXMinos;
   MinosError fYMinos;
   unsigned int fNFcn;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ContoursError
