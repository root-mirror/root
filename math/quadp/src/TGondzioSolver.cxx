// @(#)root/quadp:$Id$
// Author: Eddy Offermann   May 2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*************************************************************************
 * Parts of this file are copied from the OOQP distribution and          *
 * are subject to the following license:                                 *
 *                                                                       *
 * COPYRIGHT 2001 UNIVERSITY OF CHICAGO                                  *
 *                                                                       *
 * The copyright holder hereby grants you royalty-free rights to use,    *
 * reproduce, prepare derivative works, and to redistribute this software*
 * to others, provided that any changes are clearly documented. This     *
 * software was authored by:                                             *
 *                                                                       *
 *   E. MICHAEL GERTZ      gertz@mcs.anl.gov                             *
 *   Mathematics and Computer Science Division                           *
 *   Argonne National Laboratory                                         *
 *   9700 S. Cass Avenue                                                 *
 *   Argonne, IL 60439-4844                                              *
 *                                                                       *
 *   STEPHEN J. WRIGHT     swright@cs.wisc.edu                           *
 *   Computer Sciences Department                                        *
 *   University of Wisconsin                                             *
 *   1210 West Dayton Street                                             *
 *   Madison, WI 53706   FAX: (608)262-9777                              *
 *                                                                       *
 * Any questions or comments may be directed to one of the authors.      *
 *                                                                       *
 * ARGONNE NATIONAL LABORATORY (ANL), WITH FACILITIES IN THE STATES OF   *
 * ILLINOIS AND IDAHO, IS OWNED BY THE UNITED STATES GOVERNMENT, AND     *
 * OPERATED BY THE UNIVERSITY OF CHICAGO UNDER PROVISION OF A CONTRACT   *
 * WITH THE DEPARTMENT OF ENERGY.                                        *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
///
/// \class TGondzioSolver
///
/// Derived class of TQpSolverBase implementing Gondzio-correction
/// version of Mehrotra's original predictor-corrector algorithm.
///
////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "TMath.h"
#include "TGondzioSolver.h"
#include "TQpLinSolverDens.h"

ClassImp(TGondzioSolver);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGondzioSolver::TGondzioSolver()
{
   fPrintlevel               = 0;
   fTsig                     = 0.0;
   fMaximum_correctors       = 0;
   fNumberGondzioCorrections = 0;

   fStepFactor0 = 0.0;
   fStepFactor1 = 0.0;
   fAcceptTol   = 0.0;
   fBeta_min    = 0.0;
   fBeta_max    = 0.0;

   fCorrector_step  = 0;
   fStep            = 0;
   fCorrector_resid = 0;
   fFactory         = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGondzioSolver::TGondzioSolver(TQpProbBase *of,TQpDataBase *prob,Int_t verbose)
{
   fFactory = of;
   fStep            = fFactory->MakeVariables(prob);
   fCorrector_step  = fFactory->MakeVariables(prob);
   fCorrector_resid = fFactory->MakeResiduals(prob);

   fPrintlevel = verbose;
   fTsig       = 3.0;            // the usual value for the centering exponent (tau)

   fMaximum_correctors = 3;      // maximum number of Gondzio correctors

   fNumberGondzioCorrections = 0;

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   fStepFactor0 = 0.08;
   fStepFactor1 = 1.08;

   // accept the enhanced step if it produces a small improvement in
   // the step length
   fAcceptTol = 0.005;

   //define the Gondzio correction box
   fBeta_min = 0.1;
   fBeta_max = 10.0;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TGondzioSolver::TGondzioSolver(const TGondzioSolver &another) : TQpSolverBase(another)
{
   *this = another;
}


////////////////////////////////////////////////////////////////////////////////
/// Solve the quadratic programming problem as formulated through prob, store
/// the final solution in iterate->fX . Monitor the residuals during the iterations
/// through resid . The status is returned as defined in TQpSolverBase::ETerminationCode .

Int_t TGondzioSolver::Solve(TQpDataBase *prob,TQpVar *iterate,TQpResidual *resid)
{
   Int_t status_code;
   Double_t alpha = 1;
   Double_t sigma = 1;

   fDnorm = prob->DataNorm();

   // initialization of (x,y,z) and factorization routine.
   fSys = fFactory->MakeLinSys(prob);
   this->Start(fFactory,iterate,prob,resid,fStep);

   fIter = 0;
   fNumberGondzioCorrections = 0;
   Double_t mu = iterate->GetMu();

   Int_t done = 0;
   do {
      fIter++;
      // evaluate residuals and update algorithm status:
      resid->CalcResids(prob,iterate);

      //  termination test:
      status_code = this->DoStatus(prob,iterate,resid,fIter,mu,0);
      if (status_code != kNOT_FINISHED ) break;
      if (fPrintlevel >= 10)
         this->DoMonitor(prob,iterate,resid,alpha,sigma,fIter,mu,status_code,0);

      // *** Predictor step ***

      resid->Set_r3_xz_alpha(iterate,0.0);

      fSys->Factor(prob,iterate);
      fSys->Solve(prob,iterate,resid,fStep);
      fStep->Negate();

      alpha = iterate->StepBound(fStep);

      // calculate centering parameter
      Double_t muaff = iterate->MuStep(fStep,alpha);
      sigma = TMath::Power(muaff/mu,fTsig);

      if (fPrintlevel >= 10)
         this->DoMonitor(prob,iterate,resid,alpha,sigma,fIter,mu,status_code,2);

      // *** Corrector step ***

      // form right hand side of linear system:
      resid->Add_r3_xz_alpha(fStep,-sigma*mu);

      fSys->Solve(prob,iterate,resid,fStep);
      fStep->Negate();

      // calculate distance to boundary along the Mehrotra
      // predictor-corrector step:
      alpha = iterate->StepBound(fStep);

      // prepare for Gondzio corrector loop: zero out the
      // corrector_resid structure:
      fCorrector_resid->Clear_r1r2();

      // calculate the target box:
      Double_t rmin = sigma*mu*fBeta_min;
      Double_t rmax = sigma*mu*fBeta_max;

      Int_t stopCorrections     = 0;
      fNumberGondzioCorrections = 0;

      // enter the Gondzio correction loop:
      if (fPrintlevel >= 10)
         std::cout << "**** Entering the correction loop ****" << std::endl;

      while (fNumberGondzioCorrections < fMaximum_correctors  &&
      alpha < 1.0 && !stopCorrections) {

         // copy current variables into fcorrector_step
         *fCorrector_step = *iterate;

         // calculate target steplength
         Double_t alpha_target = fStepFactor1*alpha+fStepFactor0;
         if (alpha_target > 1.0) alpha_target = 1.0;

         // add a step of this length to corrector_step
         fCorrector_step->Saxpy(fStep,alpha_target);

         // place XZ into the r3 component of corrector_resids
         fCorrector_resid->Set_r3_xz_alpha(fCorrector_step,0.0);

         // do the projection operation
         fCorrector_resid->Project_r3(rmin,rmax);

         // solve for corrector direction
         fSys->Solve(prob,iterate,fCorrector_resid,fCorrector_step);

         // add the current step to corrector_step, and calculate the
         // step to boundary along the resulting direction
         fCorrector_step->Saxpy(fStep,1.0);
         Double_t alpha_enhanced = iterate->StepBound(fCorrector_step);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if (alpha_enhanced == 1.0) {
            *fStep = *fCorrector_step;
            alpha = alpha_enhanced;
            fNumberGondzioCorrections++;
            stopCorrections = 1;
         }
         else if(alpha_enhanced >= (1.0+fAcceptTol)*alpha) {
            // if enhanced step length is significantly better than the
            // current alpha, make the enhanced step official, but maybe
            // keep correcting
            *fStep = *fCorrector_step;
            alpha = alpha_enhanced;
            fNumberGondzioCorrections++;
            stopCorrections = 0;
         }
         else {
            // otherwise quit the correction loop
            stopCorrections = 1;
         }
      }

      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.x
      alpha = this->FinalStepLength(iterate,fStep);

      // alternatively, just use a crude step scaling factor.
      // alpha = 0.995 * iterate->StepBound(fStep);

      // actually take the step (at last!) and calculate the new mu

      iterate->Saxpy(fStep,alpha);
      mu = iterate->GetMu();
   } while (!done);

   resid->CalcResids(prob,iterate);
   if (fPrintlevel >= 10)
      this->DoMonitor(prob,iterate,resid,alpha,sigma,fIter,mu,status_code,1);

   return status_code;
}


////////////////////////////////////////////////////////////////////////////////
/// Print information about the optimization process and monitor the convergence
/// status of thye algorithm

void TGondzioSolver::DefMonitor(TQpDataBase* /* data */,TQpVar* /* vars */,
                                TQpResidual *resid,
                                Double_t alpha,Double_t sigma,Int_t i,Double_t mu,
                                Int_t status_code,Int_t level)
{
   switch (level) {
      case 0 : case 1:
      {
         std::cout << std::endl << "Duality Gap: " << resid->GetDualityGap() << std::endl;
         if (i > 1) {
            std::cout << " Number of Corrections = " << fNumberGondzioCorrections
               << " alpha = " << alpha << std::endl;
         }
         std::cout << " *** Iteration " << i << " *** " << std::endl;
         std::cout << " mu = " << mu << " relative residual norm = "
            << resid->GetResidualNorm()/fDnorm << std::endl;

         if (level == 1) {
            // Termination has been detected by the status check; print
            // appropriate message
            if (status_code == kSUCCESSFUL_TERMINATION) {
               std::cout << std::endl
                  << " *** SUCCESSFUL TERMINATION ***"
                  << std::endl;
            }
            else if (status_code == kMAX_ITS_EXCEEDED) {
               std::cout << std::endl
                  << " *** MAXIMUM ITERATIONS REACHED *** " << std::endl;
            }
            else if (status_code == kINFEASIBLE) {
               std::cout << std::endl
                  << " *** TERMINATION: PROBABLY INFEASIBLE *** "
                  << std::endl;
            }
            else if (status_code == kUNKNOWN) {
               std::cout << std::endl
                  << " *** TERMINATION: STATUS UNKNOWN *** " << std::endl;
            }
         }
      } break;
      case 2:
         std::cout << " *** sigma = " << sigma << std::endl;
         break;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Deconstructor

TGondzioSolver::~TGondzioSolver()
{
   if (fCorrector_step)  { delete fCorrector_step;  fCorrector_step  = 0; }
   if (fStep)            { delete fStep;            fStep            = 0; }
   if (fCorrector_resid) { delete fCorrector_resid; fCorrector_resid = 0; }
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TGondzioSolver &TGondzioSolver::operator=(const TGondzioSolver &source)
{
   if (this != &source) {
      TQpSolverBase::operator=(source);

      fPrintlevel               = source.fPrintlevel;
      fTsig                     = source.fTsig ;
      fMaximum_correctors       = source.fMaximum_correctors;
      fNumberGondzioCorrections = source.fNumberGondzioCorrections;

      fStepFactor0 = source.fStepFactor0;
      fStepFactor1 = source.fStepFactor1;
      fAcceptTol   = source.fAcceptTol;
      fBeta_min    = source.fBeta_min;
      fBeta_max    = source.fBeta_max;

      if (fCorrector_step)  delete fCorrector_step;
      if (fStep)            delete fStep;
      if (fCorrector_resid) delete fCorrector_resid;

      fCorrector_step  = new TQpVar(*source.fCorrector_step);
      fStep            = new TQpVar(*source.fStep);
      fCorrector_resid = new TQpResidual(*source.fCorrector_resid);
      fFactory         = source.fFactory;
   }
   return *this;
}
