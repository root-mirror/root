#ifndef PlotFoams__HH
#define PlotFoams__HH

#include "tmvaglob.h"
#include "TControlBar.h"
#include "TMap.h"
#include "TVectorT.h"
#include "TString.h"
#include "TPaveText.h"
#include "TMVA/PDEFoamKernelBase.h"
#include "TMVA/PDEFoamKernelTrivial.h"

#include <cfloat>

#include "TMVA/PDEFoam.h"

namespace TMVA{

   void PlotFoams( TString fileName = "weights/TMVAClassification_PDEFoam.weights_foams.root",
                   bool useTMVAStyle = kTRUE );
   // foam plotting macro
   void Plot(TString fileName, TMVA::ECellValue cv, TString cv_long, bool useTMVAStyle = kTRUE);

   void Plot1DimFoams(TList& foam_list, TMVA::ECellValue cell_value,
                      const TString& cell_value_description,
                      TMVA::PDEFoamKernelBase* kernel);
   void PlotNDimFoams(TList& foam_list, TMVA::ECellValue cell_value,
                      const TString& cell_value_description,
                      TMVA::PDEFoamKernelBase* kernel);
   void PlotCellTree(TString fileName, TString cv_long, bool useTMVAStyle = kTRUE);

   void DrawCell( TMVA::PDEFoamCell *cell, TMVA::PDEFoam *foam,
                  Double_t x, Double_t y,
                  Double_t xscale,  Double_t yscale );
}

#endif
