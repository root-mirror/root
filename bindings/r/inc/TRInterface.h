// @(#)root/r:$Id$
// Author: Omar Zapata   29/05/2013


/*************************************************************************
 * Copyright (C) 2013-2014, Omar Andres Zapata Mesa                      *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_R_TRInterface
#define ROOT_R_TRInterface

#include <TRObject.h>

#include <TRDataFrame.h>

#include <TRFunctionExport.h>

#include <TRFunctionImport.h>

#include <TThread.h>

/**
 @namespace ROOT::R
 namespace associated R package for ROOT.
 @defgroup R R Interface for Statistical Computing
 */
namespace ROOT {
   namespace R {
      /**
      \class TRInterface
         ROOT R was implemented using the
         <A HREF="http://www.r-project.org/">R Project</A> library and the modules
         <A HREF="http://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</A> and
         <A HREF="http://cran.r-project.org/web/packages/RInside/index.html">RInside</A>
         <h2>Users Guide </h2>
         <a href="http://oproject.org/tiki-index.php?page=ROOT+R+Users+Guide"> http://oproject.org/tiki-index.php?page=ROOT+R+Users+Guide</a><br>
         <a href="https://root.cern.ch/drupal/content/how-use-r-root-root-r-interface"> https://root.cern.ch/drupal/content/how-use-r-root-root-r-interface</a>

         \ingroup R
       */


      /**
       <center><h2>TRInterface class</h2></center>

      </p>
      The TRInterface class lets you procces R code from ROOT.<br>
      You can call R libraries and their functions, plot results in R or ROOT,<br>
      and use the power of ROOT and R at the same time.<br>
      It also lets you pass scalars, vectors and matrices from ROOT to R<br>
      and from R to ROOT using TRObject; but you can to use overloaded opetarors [],<< and >> <br>
      to work with ROOTR like work with streams of data.<br>

      TRInterface class can not be instantiated directly, but you can create objects using the static methods
      TRInterface& Instance() and TRInterface* InstancePtr() to create your own objects.<br>
      <br>
      </p>
      Show an example below:
      Create an exponential fit, the idea is to create a set of numbers x,y with noise from ROOT,
      pass them to R and fit the data to \f$ x^3 \f$, get the fitted coefficient(power) and plot the data,
      the known function and the fitted function.
      \code{.cpp}

         TCanvas *c1 = new TCanvas("c1","Curve Fit",700,500);
         c1->SetGrid();

         // draw a frame for multiples graphs
         TMultiGraph *mg = new TMultiGraph();

         // create the first graph (points with gaussian noise)
         const Int_t n = 24;
         Double_t x[n] ;
         Double_t y[n] ;
         //Generate points along a X^3 with noise
         TRandom rg;
         rg.SetSeed(520);
         for (Int_t i = 0; i < n; i++) {
            x[i] = rg.Uniform(0, 1);
            y[i] = TMath::Power(x[i], 3) + rg.Gaus() * 0.06;
         }

         TGraph *gr1 = new TGraph(n,x,y);
         gr1->SetMarkerColor(kBlue);
         gr1->SetMarkerStyle(8);
         gr1->SetMarkerSize(1);
         mg->Add(gr1);

            // create second graph
         TF1 *f_known=new TF1("f_known","pow(x,3)",0,1);
         TGraph *gr2 = new TGraph(f_known);
         gr2->SetMarkerColor(kRed);
         gr2->SetMarkerStyle(8);
         gr2->SetMarkerSize(1);
         mg->Add(gr2);

         //passing x and y values to R for fitting
         ROOT::R::TRInterface &r=ROOT::R::TRInterface::Instance();
         r["x"]<<TVectorD(n, x);
         r["y"]<<TVectorD(n, y);
         //creating a R data frame
         r<<"ds<-data.frame(x=x,y=y)";
         //fitting x and y to X^power using Nonlinear Least Squares
         r<<"m <- nls(y ~ I(x^power),data = ds, start = list(power = 1),trace = T)";
         //getting the fitted value (power)
         Double_t power;
         r["summary(m)$coefficients[1]"]>>power;

         TF1 *f_fitted=new TF1("f_fitted","pow(x,[0])",0,1);
         f_fitted->SetParameter(0,power);
         //plotting the fitted function
         TGraph *gr3 = new TGraph(f_fitted);
         gr3->SetMarkerColor(kGreen);
         gr3->SetMarkerStyle(8);
         gr3->SetMarkerSize(1);

         mg->Add(gr3);
         mg->Draw("ap");

         //displaying basic results
         TPaveText *pt = new TPaveText(0.1,0.6,0.5,0.9,"brNDC");
         pt->SetFillColor(18);
         pt->SetTextAlign(12);
         pt->AddText("Fitting x^power ");
         pt->AddText(" \"Blue\"   Points with gaussian noise to be fitted");
         pt->AddText(" \"Red\"    Known function x^3");
         TString fmsg;
         fmsg.Form(" \"Green\"  Fitted function with power=%.4lf",power);
         pt->AddText(fmsg);
         pt->Draw();
         c1->Update();
      \endcode
      @ingroup R
       */
      class TRInterface: public TObject {
      protected:
         RInside *fR;
         TThread *th;
      public:
         //Proxy class to use operators for assignation Ex: r["name"]=object
         class Binding {
         public:
            Binding(TRInterface *rnt, TString name): fInterface(rnt), fName(name) {}
            Binding &operator=(const Binding &obj)
            {
               fInterface = obj.fInterface;
               fName = obj.fName;
               return *this;
            }
            template <class T> Binding &operator=(const T &data)
            {
               fInterface->Assign<T>(data, fName);
               return *this;
            }
            Binding &operator=(const TRFunctionExport &fun)
            {
               //The method assign is not a template for a function
               fInterface->Assign(fun, fName);
               return *this;
            }

            Binding &operator<<(const TRFunctionExport &fun)
            {
               //The method assign is not a template for a function
               fInterface->Assign(fun, fName);
               return *this;
            }

            Binding &operator=(const TRDataFrame &df)
            {
               fInterface->Assign(df, fName);
               return *this;
            }

            Binding &operator<<(const TRDataFrame &df)
            {
               fInterface->Assign(df, fName);
               return *this;
            }

            template <class T> Binding &operator >>(T &var)
            {
               var = fInterface->Eval(fName).As<T>();
               return *this;
            }

            template <class T> Binding &operator <<(T var)
            {
               fInterface->Assign<T>(var, fName);
               return *this;
            }
#include<TRInterface_Binding.h>
            template <class T> operator T()
            {
               return fInterface->Eval(fName);
            }

         private:
            TRInterface *fInterface;
            TString fName;
         };
      private:
         /**
         The command line arguments are by deafult argc=0 and argv=NULL,
         The verbose mode is by default disabled but you can enable it to show procedures information in stdout/stderr         \note some time can produce so much noise in the output
         \param argc default 0
         \param argv default null
         \param loadRcpp default true
         \param verbose default false
         \param interactive default true
         */
         TRInterface(const Int_t argc = 0, const Char_t *argv[] = NULL, const Bool_t loadRcpp = true,
                     const Bool_t verbose = false, const Bool_t interactive = true);

      public:
         ~TRInterface();

         /**
         Method to set verbose mode, that produce extra output
         \note some time can produce so much noise in the output
         \param status boolean to enable of disable
         */
         void SetVerbose(Bool_t status);
         /**
         Method to eval R code and you get the result in a reference to TRObject
         \param code R code
         \param ans reference to TRObject
         \return an true or false if the execution was sucessful or not.
         */
         Int_t Eval(const TString &code, TRObject  &ans); // parse line, returns in ans; error code rc
         /**
         Method to eval R code
         \param code R code
         */
         void  Execute(const TString &code);

         // "unhide" TObject::Execute methods.
         using TObject::Execute;

         /**
         Method to eval R code and you get the result in a  TRObject
         \param code R code
         \return a TRObject with result
         */
         TRObject Eval(const TString &code);


         /**
         Template method to assign C++ variables into R enviroment
         \param var any R wrappable datatype
         \param name  name of the variable in R's enviroment
         */
         template<typename T >void Assign(const T &var, const TString &name)
         {
            // This method lets you pass variables from ROOT to R.
            // The template T should be a supported ROOT datatype and
            // the TString's name is the name of the variable in the R enviroment.
            fR->assign<T>(var, name.Data());
         }
         /**
         Method to assign TRFunctionExport in R's enviroment
         \param fun TRFunctionExport
         \param name  name of the variable in R's enviroment
         */
         void Assign(const TRFunctionExport &fun, const TString &name);
         /**
         Method to assign TRDataFrame in R's enviroment
         \param df TRDataFrame
         \param name  name of the variable in R's enviroment
         */
         void Assign(const TRDataFrame &df, const TString &name);

         /**
         Method to get a R prompt to work interactively with tab completation support
         */
         void Interactive();

         /**
         Init event loop in a thread to support actions in windows from R graphics system
         */
         void ProcessEventsLoop();

         /**
         Method to verify if a package is installed
         \param pkg R's pkg name
         \return true or false if the package is installed or not
         */
         Bool_t IsInstalled(TString pkg);
         /**
         Method to load an R's package
         \param pkg R's pkg name
         \return true or false if the package was loaded or not
         */
         Bool_t Require(TString pkg);
         /**
         Method to install an R's package
         \param pkg R's pkg name
         \param repos url for R's package repository
         \return true or false if the package was installed or not
         */
         Bool_t Install(TString pkg, TString repos = "http://cran.r-project.org");
         Binding operator[](const TString &name);

         /**
         static method to get an TRInterface instance reference
         \return TRInterface instance reference
         */
         static TRInterface &Instance();
         /**
         static method to get an TRInterface instance pointer
         \return TRInterface instance pointer
         */
         static TRInterface *InstancePtr();

         ClassDef(TRInterface, 0)
      };
   }
}

inline ROOT::R::TRInterface &operator<<(ROOT::R::TRInterface &r, TString code)
{
   r.Execute(code);
   return r;
}

#endif
