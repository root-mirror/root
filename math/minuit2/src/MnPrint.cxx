// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnPrint.h"
#include "Minuit2/LAVector.h"
#include "Minuit2/LASymMatrix.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnGlobalCorrelationCoeff.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MinuitParameter.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include "Minuit2/MnPlot.h"

#include <iomanip>
#include <utility>

constexpr int PRECISION = 10;
constexpr int WIDTH = PRECISION + 7;

namespace ROOT {
namespace Minuit2 {

int gPrintLevel = 0;

int MnPrint::SetGlobalLevel(int level)
{
   std::swap(level, gPrintLevel);
   return level;
}

int MnPrint::GlobalLevel()
{
   return gPrintLevel;
}

MnPrint::MnPrint(const char *prefix, int level) : fPrefix{prefix}, fLevel{level} {}

int MnPrint::SetLevel(int level)
{
   std::swap(level, fLevel);
   return level;
}

int MnPrint::Level() const
{
   return fLevel;
}

MnPrint::Oneline::Oneline(double fcn, double edm, int ncalls, int iter)
   : fFcn(fcn), fEdm(edm), fNcalls(ncalls), fIter(iter)
{
}

MnPrint::Oneline::Oneline(const MinimumState &state, int iter)
   : MnPrint::Oneline(state.Fval(), state.Edm(), state.NFcn(), iter)
{
}

MnPrint::Oneline::Oneline(const FunctionMinimum &fmin, int iter) : MnPrint::Oneline(fmin.State(), iter) {}

std::ostream &operator<<(std::ostream &os, const MnPrint::Oneline &x)
{
   // print iteration, function value, edm and ncalls in one single line
   if (x.fIter >= 0)
      os << std::setw(4) << x.fIter << " - ";
   const int pr = os.precision(PRECISION);
   os << "FCN = " << std::setw(WIDTH) << x.fFcn << " Edm = " << std::setw(WIDTH) << x.fEdm << " NCalls = " << std::setw(6)
      << x.fNcalls;
   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const LAVector &vec)
{
   // print a vector
   const int pr = os.precision(PRECISION);
   const int nrow = vec.size();
   for (int i = 0; i < nrow; i++) {
      os << '\n';
      os.width(WIDTH);
      os << vec(i);
   }
   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const LASymMatrix &matrix)
{
   // print a matrix
   const int pr = os.precision(8);
   const int n = matrix.Nrow();
   for (int i = 0; i < n; i++) {
      os << '\n';
      for (int j = 0; j < n; j++) {
         os.width(15);
         os << matrix(i, j);
      }
   }
   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const MnUserParameters &par)
{
   // print the MnUserParameter object
   os << "\n  Pos |    Name    |  type   |      Value       |    Error +/-";

   int pr = os.precision();

   const double eps2 = par.Precision().Eps2();
   for (auto &&p : par.Parameters()) {
      os << "\n" << std::setw(5) << p.Number() << " | " << std::setw(10) << p.Name() << " |";
      if (p.IsConst())
         os << "  const  |";
      else if (p.IsFixed())
         os << "  fixed  |";
      else if (p.HasLimits())
         os << " limited |";
      else
         os << "  free   |";
      os.precision(PRECISION);
      os.width(WIDTH);
      os << p.Value() << " | " << std::setw(12);
      if (p.Error() > 0) {
        os << p.Error();
        if (p.HasLimits()) {
           if (std::fabs(p.Value() - p.LowerLimit()) < eps2) {
              os << " (at lower limit)";
           } else if (std::fabs(p.Value() - p.UpperLimit()) < eps2) {
              os << " (at upper limit)";
           }
        }
      }
   }
   os.precision(pr);

   return os;
}

std::ostream &operator<<(std::ostream &os, const MnUserCovariance &matrix)
{
   // print the MnUserCovariance
   const int pr = os.precision(6);

   unsigned int n = matrix.Nrow();
   for (unsigned int i = 0; i < n; i++) {
      os << '\n';
      for (unsigned int j = 0; j < n; j++) {
         os.width(13);
         os << matrix(i, j);
      }
      os << " | ";
      double di = matrix(i, i);
      for (unsigned int j = 0; j < n; j++) {
         double dj = matrix(j, j);
         os.width(13);
         os << matrix(i, j) / std::sqrt(std::fabs(di * dj));
      }
   }
   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const MnGlobalCorrelationCoeff &coeff)
{
   // print the global correlation coefficient
   const int pr = os.precision(6);
   for (unsigned int i = 0; i < coeff.GlobalCC().size(); i++) {
      os << '\n';
      os.width(6 + 7);
      os << coeff.GlobalCC()[i];
   }
   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const MnUserParameterState &state)
{
   // print the MnUserParameterState
   const int pr = os.precision(PRECISION);
   os << "\n  Valid         : " << (state.IsValid() ? "yes" : "NO") << "\n  Function calls: " << state.NFcn()
      << "\n  Minimum value : " << state.Fval() << "\n  Edm           : " << state.Edm()
      << "\n  Parameters    : " << state.Parameters() << "\n  CovarianceStatus: " << state.CovarianceStatus()
      << "\n  Covariance and correlation matrix: ";
   if (state.HasCovariance())
      os << state.Covariance();
   else
      os << "matrix is not present or not valid";
   if (state.HasGlobalCC())
      os << "\n  Global correlation coefficients: " << state.GlobalCC();

   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const FunctionMinimum &min)
{
   // print the FunctionMinimum
   const int pr = os.precision(PRECISION);
   os << "\n  Valid         : " << (min.IsValid() ? "yes" : "NO") << "\n  Function calls: " << min.NFcn()
      << "\n  Minimum value : " << min.Fval() << "\n  Edm           : " << min.Edm()
      << "\n  Internal parameters: " << min.Parameters().Vec();
   if (min.HasValidCovariance())
      os << "\n  Internal covariance matrix: " << min.Error().Matrix();
   os << "\n  External parameters: " << min.UserParameters();
   // os << min.UserCovariance() << '\n';
   // os << min.UserState().GlobalCC() << '\n';

   if (!min.IsValid()) {
      os << "\n  FunctionMinimum is invalid:";
      if (!min.State().IsValid())
         os << "\n    State is invalid";
      if (min.IsAboveMaxEdm())
         os << "\n    Edm is above max";
      if (min.HasReachedCallLimit())
         os << "\n    Reached call limit";
   }

   os.precision(pr);

   return os;
}

std::ostream &operator<<(std::ostream &os, const MinimumState &min)
{
   const int pr = os.precision(PRECISION);
   os << "\n  Minimum value : " << min.Fval() << "\n  Edm           : " << min.Edm()
      << "\n  Internal parameters:" << min.Vec() << "\n  Internal gradient  :" << min.Gradient().Vec();
   if (min.HasCovariance())
      os << "\n  Internal covariance matrix:" << min.Error().Matrix();
   os.precision(pr);
   return os;
}

std::ostream &operator<<(std::ostream &os, const MnMachinePrecision &prec)
{
   // print the Precision
   int pr = os.precision(PRECISION);
   os << "MnMachinePrecision " << prec.Eps() << '\n';
   os.precision(pr);

   return os;
}

std::ostream &operator<<(std::ostream &os, const MinosError &me)
{
   // print the Minos Error
   os << "Minos # of function calls: " << me.NFcn() << '\n';

   if (!me.IsValid())
      os << "Minos Error is not valid." << '\n';
   if (!me.LowerValid())
      os << "lower Minos Error is not valid." << '\n';
   if (!me.UpperValid())
      os << "upper Minos Error is not valid." << '\n';
   if (me.AtLowerLimit())
      os << "Minos Error is Lower limit of Parameter " << me.Parameter() << "." << '\n';
   if (me.AtUpperLimit())
      os << "Minos Error is Upper limit of Parameter " << me.Parameter() << "." << '\n';
   if (me.AtLowerMaxFcn())
      os << "Minos number of function calls for Lower Error exhausted." << '\n';
   if (me.AtUpperMaxFcn())
      os << "Minos number of function calls for Upper Error exhausted." << '\n';
   if (me.LowerNewMin()) {
      os << "Minos found a new Minimum in negative direction." << '\n';
      os << me.LowerState() << '\n';
   }
   if (me.UpperNewMin()) {
      os << "Minos found a new Minimum in positive direction." << '\n';
      os << me.UpperState() << '\n';
   }

   int pr = os.precision();

   os << "No  |"
      << "|   Name    |"
      << "|   Value@min   |"
      << "|    negative   |"
      << "|   positive  " << '\n';
   os << std::setw(4) << me.Parameter() << std::setw(5) << "||";
   os << std::setw(10) << me.LowerState().Name(me.Parameter()) << std::setw(3) << "||";
   os << std::setprecision(PRECISION) << std::setw(WIDTH) << me.Min() << " ||" << std::setprecision(PRECISION)
      << std::setw(WIDTH) << me.Lower() << " ||" << std::setw(WIDTH) << me.Upper() << '\n';

   os << '\n';
   os.precision(pr);

   return os;
}

std::ostream &operator<<(std::ostream &os, const ContoursError &ce)
{
   // print the ContoursError
   os << "Contours # of function calls: " << ce.NFcn() << '\n';
   os << "MinosError in x: " << '\n';
   os << ce.XMinosError() << '\n';
   os << "MinosError in y: " << '\n';
   os << ce.YMinosError() << '\n';
   MnPlot plot;
   plot(ce.XMin(), ce.YMin(), ce());
   for (auto ipar = ce().begin(); ipar != ce().end(); ++ipar) {
      os << ipar - ce().begin() << "  " << (*ipar).first << "  " << (*ipar).second << '\n';
   }
   os << '\n';

   return os;
}

std::ostream &operator<<(std::ostream &os, const std::pair<double, double> &point)
{
   os << "\t x = " << point.first << "  y = " << point.second << std::endl;
   return os;
}

} // namespace Minuit2
} // namespace ROOT
