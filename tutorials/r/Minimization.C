/// \file
/// \ingroup tutorial_r
/// \notebook -nodraw
/// Example based  in
/// http://root.cern/root/html/tutorials/fit/NumericalMinimization.C.html
/// http://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
///
/// \macro_code
///
/// \author Omar Zapata

#include<TRInterface.h>

//in the next function the *double pointer must be changed by a TVectorD,
//because the pointer has no meaning in R enviroment.
Double_t RosenBrock(const TVectorD xx )
{
   const Double_t x = xx[0];
   const Double_t y = xx[1];
   const Double_t tmp1 = y-x*x;
   const Double_t tmp2 = 1-x;
   return 100*tmp1*tmp1+tmp2*tmp2;
}

TVectorD RosenBrockGrad(const TVectorD xx )
{
   const Double_t x = xx[0];
   const Double_t y = xx[1];
   TVectorD grad(2);
   grad[0]=-400 * x * (y - x * x) - 2 * (1 - x);
   grad[1]=200 * (y - x * x);
   return grad;
}

void Minimization()
{
   ROOT::R::TRInterface &r=ROOT::R::TRInterface::Instance();

   //passsing RosenBrock function to R
   r["RosenBrock"]=ROOT::R::TRFunctionExport(RosenBrock);

   //passsing RosenBrockGrad function to R
   r["RosenBrockGrad"]=ROOT::R::TRFunctionExport(RosenBrockGrad);
   //the option "method" could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

   //the option "control" lets you put some constraints like
   //"maxit" The maximum number of iterations.
   //"abstol" The absolute convergence tolerance.
   r.Execute("result <- optim( c(0.01,0.01), RosenBrock,method='BFGS',control = list(maxit = 1000000) )");
   //"reltol" Relative convergence tolerance.

   //Getting results from R
   TVectorD  min=r.Eval("result$par");

   std::cout.precision(8);
   //printing results
   std::cout<<"-----------------------------------------"<<std::endl;
   std::cout<<"Minimum x="<<min[0]<<" y="<<min[1]<<std::endl;
   std::cout<<"Value at minimum ="<<RosenBrock(min)<<std::endl;

   //using the gradient
   r.Execute("optimHess(result$par, RosenBrock, RosenBrockGrad)");
   r.Execute("hresult <- optim(c(-1.2,1), RosenBrock, NULL, method = 'BFGS', hessian = TRUE)");
   //getting the min calculated with the gradient
   TVectorD  hmin=r.Eval("hresult$par");

   //printing results
   std::cout<<"-----------------------------------------"<<std::endl;
   std::cout<<"Minimization with the Gradient"<<std::endl;
   std::cout<<"Minimum x="<<hmin[0]<<" y="<<hmin[1]<<std::endl;
   std::cout<<"Value at minimum ="<<RosenBrock(hmin)<<std::endl;
}
