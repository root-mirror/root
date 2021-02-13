// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RuleFitAPI                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Interface to Friedman's RuleFit method                                    *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker    <Andreas.Hocker@cern.ch>     - CERN, Switzerland       *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *      Helge Voss         <Helge.Voss@cern.ch>         - MPI-KP Heidelberg, Ger. *
 *      Kai Voss           <Kai.Voss@cern.ch>           - U. of Victoria, Canada  *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-KP Heidelberg, Germany                                                *
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RuleFitAPI
#define ROOT_TMVA_RuleFitAPI

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RuleFitAPI                                                           //
//                                                                      //
// J Friedman's RuleFit method                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <vector>

#include "TMVA/MsgLogger.h"

namespace TMVA {

   class MethodRuleFit;
   class RuleFit;

   class RuleFitAPI {

   public:

      RuleFitAPI( const TMVA::MethodRuleFit *rfbase, TMVA::RuleFit *rulefit, EMsgType minType );

      virtual ~RuleFitAPI();

      // welcome message
      void WelcomeMessage();

      // message on howto get the binary
      void HowtoSetupRF();

      // Set RuleFit working directory
      void SetRFWorkDir(const char * wdir);

      // Check RF work dir - aborts if it fails
      void CheckRFWorkDir();

      // run rf_go.exe in various modes
      inline void TrainRuleFit();
      inline void TestRuleFit();
      inline void VarImp();

      // read result into MethodRuleFit
      Bool_t ReadModelSum();

      // Get working directory
      const TString GetRFWorkDir() const { return fRFWorkDir; }

   protected:

      enum ERFMode    { kRfRegress=1, kRfClass=2 };          // RuleFit modes, default=Class
      enum EModel     { kRfLinear=0, kRfRules=1, kRfBoth=2 }; // models, default=Both (rules+linear)
      enum ERFProgram { kRfTrain=0, kRfPredict, kRfVarimp };    // rf_go.exe running mode

      // integer parameters
      typedef struct {
         Int_t mode;
         Int_t lmode;
         Int_t n;
         Int_t p;
         Int_t max_rules;
         Int_t tree_size;
         Int_t path_speed;
         Int_t path_xval;
         Int_t path_steps;
         Int_t path_testfreq;
         Int_t tree_store;
         Int_t cat_store;
      } IntParms;

      // float parameters
      typedef struct {
         Float_t  xmiss;
         Float_t  trim_qntl;
         Float_t  huber;
         Float_t  inter_supp;
         Float_t  memory_par;
         Float_t  samp_fract;
         Float_t  path_inc;
         Float_t  conv_fac;
      } RealParms;

      // setup
      void InitRuleFit();
      void FillRealParmsDef();
      void FillIntParmsDef();
      void ImportSetup();
      void SetTrainParms();
      void SetTestParms();

      // run
      Int_t  RunRuleFit();

      // set rf_go.exe running mode
      void SetRFTrain()   { fRFProgram = kRfTrain; }
      void SetRFPredict() { fRFProgram = kRfPredict; }
      void SetRFVarimp()  { fRFProgram = kRfVarimp; }

      // handle rulefit files
      inline TString GetRFName(TString name);
      inline Bool_t  OpenRFile(TString name, std::ofstream & f);
      inline Bool_t  OpenRFile(TString name, std::ifstream & f);

      // read/write binary files
      inline Bool_t WriteInt(std::ofstream &   f, const Int_t   *v, Int_t n=1);
      inline Bool_t WriteFloat(std::ofstream & f, const Float_t *v, Int_t n=1);
      inline Int_t  ReadInt(std::ifstream & f,   Int_t *v, Int_t n=1) const;
      inline Int_t  ReadFloat(std::ifstream & f, Float_t *v, Int_t n=1) const;

      // write rf_go.exe i/o files
      Bool_t WriteAll();
      Bool_t WriteIntParms();
      Bool_t WriteRealParms();
      Bool_t WriteLx();
      Bool_t WriteProgram();
      Bool_t WriteRealVarImp();
      Bool_t WriteRfOut();
      Bool_t WriteRfStatus();
      Bool_t WriteRuleFitMod();
      Bool_t WriteRuleFitSum();
      Bool_t WriteTrain();
      Bool_t WriteVarNames();
      Bool_t WriteVarImp();
      Bool_t WriteYhat();
      Bool_t WriteTest();

      // read rf_go.exe i/o files
      Bool_t ReadYhat();
      Bool_t ReadIntParms();
      Bool_t ReadRealParms();
      Bool_t ReadLx();
      Bool_t ReadProgram();
      Bool_t ReadRealVarImp();
      Bool_t ReadRfOut();
      Bool_t ReadRfStatus();
      Bool_t ReadRuleFitMod();
      Bool_t ReadRuleFitSum();
      Bool_t ReadTrainX();
      Bool_t ReadTrainY();
      Bool_t ReadTrainW();
      Bool_t ReadVarNames();
      Bool_t ReadVarImp();

   private:
      // prevent empty constructor from being used
      RuleFitAPI();
      const MethodRuleFit *fMethodRuleFit; // parent method - set in constructor
      RuleFit             *fRuleFit;       // non const ptr to RuleFit class in MethodRuleFit
      //
      std::vector<Float_t> fRFYhat;      // score results from test sample
      std::vector<Float_t> fRFVarImp;    // variable importances
      std::vector<Int_t>   fRFVarImpInd; // variable index
      TString              fRFWorkDir;   // working directory
      IntParms             fRFIntParms;  // integer parameters
      RealParms            fRFRealParms; // real parameters
      std::vector<int>     fRFLx;        // variable selector
      ERFProgram           fRFProgram;   // what to run
      TString              fModelType;   // model type string

      mutable MsgLogger    fLogger;      // message logger

      ClassDef(RuleFitAPI,0);        // Friedman's RuleFit method

   };

} // namespace TMVA

//_______________________________________________________________________
void TMVA::RuleFitAPI::TrainRuleFit()
{
   // run rf_go.exe to train the model
   SetTrainParms();
   WriteAll();
   RunRuleFit();
}

//_______________________________________________________________________
void TMVA::RuleFitAPI::TestRuleFit()
{
   // run rf_go.exe with the test data
   SetTestParms();
   WriteAll();
   RunRuleFit();
   ReadYhat(); // read in the scores
}

//_______________________________________________________________________
void TMVA::RuleFitAPI::VarImp()
{
   // run rf_go.exe to get the variable importance
   SetRFVarimp();
   WriteAll();
   RunRuleFit();
   ReadVarImp(); // read in the variable importances
}

//_______________________________________________________________________
TString TMVA::RuleFitAPI::GetRFName(TString name)
{
   // get the name including the rulefit directory
   return fRFWorkDir+"/"+name;
}

//_______________________________________________________________________
Bool_t TMVA::RuleFitAPI::OpenRFile(TString name, std::ofstream & f)
{
   // open a file for writing in the rulefit directory
   TString fullName = GetRFName(name);
   f.open(fullName);
   if (!f.is_open()) {
      fLogger << kERROR << "Error opening RuleFit file for output: "
              << fullName << Endl;
      return kFALSE;
   }
   return kTRUE;
}

//_______________________________________________________________________
Bool_t TMVA::RuleFitAPI::OpenRFile(TString name, std::ifstream & f)
{
   // open a file for reading in the rulefit directory
   TString fullName = GetRFName(name);
   f.open(fullName);
   if (!f.is_open()) {
      fLogger << kERROR << "Error opening RuleFit file for input: "
              << fullName << Endl;
      return kFALSE;
   }
   return kTRUE;
}

//_______________________________________________________________________
Bool_t TMVA::RuleFitAPI::WriteInt(std::ofstream &   f, const Int_t   *v, Int_t n)
{
   // write an int
   if (!f.is_open()) return kFALSE;
   return (Bool_t)f.write(reinterpret_cast<char const *>(v), n*sizeof(Int_t));
}

//_______________________________________________________________________
Bool_t TMVA::RuleFitAPI::WriteFloat(std::ofstream & f, const Float_t *v, Int_t n)
{
   // write a float
   if (!f.is_open()) return kFALSE;
   return (Bool_t)f.write(reinterpret_cast<char const *>(v), n*sizeof(Float_t));
}

//_______________________________________________________________________
Int_t TMVA::RuleFitAPI::ReadInt(std::ifstream & f,   Int_t *v, Int_t n) const
{
   // read an int
   if (!f.is_open()) return 0;
   if (f.read(reinterpret_cast<char *>(v), n*sizeof(Int_t))) return 1;
   return 0;
}

//_______________________________________________________________________
Int_t TMVA::RuleFitAPI::ReadFloat(std::ifstream & f, Float_t *v, Int_t n) const
{
   // read a float
   if (!f.is_open()) return 0;
   if (f.read(reinterpret_cast<char *>(v), n*sizeof(Float_t))) return 1;
   return 0;
}

#endif // RuleFitAPI_H
