/// \file
/// \ingroup tutorial_r
/// \notebook -nodraw
/// More Information for R interpolation in
/// http://stat.ethz.ch/R-manual/R-patched/library/stats/html/approxfun.html
/// NOTE: this example illustrates an interpolation with random points given from ROOT
/// and procedures made in R's environment.
///
/// \macro_code
///
/// \author Omar Zapata

#include<TRInterface.h>
#include<TRandom.h>
#include<vector>

void Interpolation()
{
   ROOT::R::TRInterface &r=ROOT::R::TRInterface::Instance();
   //Creating points
   TRandom rg;
   std::vector<Double_t> x(10),y(10);
   for(int i=0;i<10;i++)
   {
      x[i]=i;
      y[i]=rg.Gaus();
   }

   r["x"]=x;
   r["y"]=y;


   // do plotting only in non-batch mode
   if (!gROOT->IsBatch() )  {

      r<<"dev.new()";//Required to activate new window for plot
      //Plot parameter. Plotting using two rows and one column
      r<<"par(mfrow = c(2,1))";

      //plotting the points
      r<<"plot(x, y, main = 'approx(.) and approxfun(.)')";

      //The function "approx" returns a list with components x and y
      //containing n coordinates which interpolate the given data points according to the method (and rule) desired.
      r<<"points(approx(x, y), col = 2, pch = '*')";
      r<<"points(approx(x, y, method = 'constant'), col = 4, pch = '*')";
   }
   else {
      r << "print('Interpolated points')";
      r << "print(approx(x,y,n=20))";
   }

   //The function "approxfun" returns a function performing (linear or constant)
   //interpolation of the given data.
   //For a given set of x values, this function will return the corresponding interpolated values.
   r<<"f <- approxfun(x, y)";
   //using approxfun with const method
   r<<"fc <- approxfun(x, y, method = 'const')";

   if (!gROOT->IsBatch() ) {
      r<<"curve(f(x), 0, 11, col = 'green2')";
      r<<"points(x, y)";

      r<<"curve(fc(x), 0, 10, col = 'darkblue', add = TRUE)";
      // different interpolation on left and right side :
      r<<"plot(approxfun(x, y, rule = 2:1), 0, 11,col = 'tomato', add = TRUE, lty = 3, lwd = 2)";
   }
   else {
      r << "x2=x+0.5";
      r << "print('Result of approxfun with default method')";
      r << "print(paste('x = ',x,'  f(x) = ',f(x2)))";
      r << "print('Result of approxfun with const method')";
      r << "print(paste('x = ',x,'  f(x) = ',fc(x2)))";
   }
}
