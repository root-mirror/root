/* @(#)root/hist:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ enum EErrorType;

#pragma link C++ class TAxis-;
#pragma link C++ class TAxisModLab+;
#pragma link C++ class TBinomialEfficiencyFitter+;
#pragma link C++ class TFormula-;
#pragma link C++ class ROOT::v5::TFormula-;
#pragma link C++ class ROOT::v5::TFormulaPrimitive+;
#pragma link C++ class TFractionFitter+;
#pragma link C++ class TFitResult+;
#pragma link C++ class TFitResultPtr+;
#pragma link C++ class TF1NormSum+;
#pragma link C++ class TF1Convolution+;
#pragma link C++ class TF1-;
#pragma link C++ class ROOT::v5::TF1Data-;
#pragma read sourceClass="TF1" targetClass="ROOT::v5::TF1Data";
#pragma read sourceClass="TFormula" targetClass="ROOT::v5::TFormula";
#pragma link C++ class TF1Parameters+;
#pragma link C++ class TFormulaParamOrder+;
#pragma link C++ class std::map<TString,int,TFormulaParamOrder>+;
#pragma link C++ class TF12+;
#pragma link C++ class TF1AbsComposition + ;
#pragma link C++ class TF1Convolution + ;
#pragma link C++ class TF1NormSum + ;
#pragma link C++ class std::vector < std::unique_ptr < TF1 >> +;
#pragma link C++ class std::vector < std::unique_ptr < TF1AbsComposition >> +;
#pragma link C++ class TF2-;
#pragma link C++ class TF3-;
#pragma link C++ class Foption_t+;
#pragma link C++ class TGraph-;
#pragma link C++ class TGraphErrors-;
#pragma link C++ class TGraphAsymmErrors-;
#pragma link C++ class TGraphMultiErrors+;
#pragma link C++ class TGraphBentErrors+;
#pragma link C++ class TGraph2D-;
#pragma link C++ class TGraph2DErrors-;
#pragma link C++ class TGraphDelaunay+;
#pragma link C++ class TGraphDelaunay2D+;
#pragma link C++ class TGraphSmooth+;
#pragma link C++ class TGraphTime+;
#pragma link C++ class TH1-;
#pragma link C++ class TH1C+;
#pragma link C++ class TH1D+;
#pragma link C++ class TH1F+;
#pragma link C++ class TH1S+;
#pragma link C++ class TH1I+;
#pragma link C++ class TH1K+;
#pragma link C++ class TH2-;
#pragma link C++ class TH2C-;
#pragma link C++ class TH2D-;
#pragma link C++ class TH2F-;
#pragma link C++ class TH2Poly+;
#pragma link C++ class TH2PolyBin+;
#pragma link C++ class THistRange+;
#pragma link C++ class TBinIterator+;
#pragma link C++ class TProfile2Poly+;
#pragma link C++ class TProfile2PolyBin+;
#pragma link C++ class TH2S-;
#pragma link C++ class TH2I+;
#pragma link C++ class TH3-;
#pragma link C++ class TH3C-;
#pragma link C++ class TH3D-;
#pragma link C++ class TH3F-;
#pragma link C++ class TH3S-;
#pragma link C++ class TH3I+;
#pragma link C++ class THLimitsFinder+;
#pragma link C++ class THnBase+;
#pragma link C++ class THnIter+;
#pragma link C++ class TNDArray+;
#pragma link C++ class TNDArrayT<Float_t>+;
//#pragma link C++ class TNDArrayT<Float16_t>+;
#pragma link C++ class TNDArrayT<Double_t>+;
//#pragma link C++ class TNDArrayT<Double32_t>+;
#pragma link C++ class TNDArrayT<Long64_t>+;
#pragma link C++ class TNDArrayT<Long_t>+;
#pragma link C++ class TNDArrayT<Int_t>+;
#pragma link C++ class TNDArrayT<Short_t>+;
#pragma link C++ class TNDArrayT<Char_t>+;
#pragma link C++ class TNDArrayT<ULong64_t>+;
#pragma link C++ class TNDArrayT<ULong_t>+;
#pragma link C++ class TNDArrayT<UInt_t>+;
#pragma link C++ class TNDArrayT<UShort_t>+;
#pragma link C++ class TNDArrayRef<Float_t>+;
//#pragma link C++ class TNDArrayRef<Float16_t>+;
#pragma link C++ class TNDArrayRef<Double_t>+;
//#pragma link C++ class TNDArrayRef<Double32_t>+;
#pragma link C++ class TNDArrayRef<Long64_t>+;
#pragma link C++ class TNDArrayRef<Long_t>+;
#pragma link C++ class TNDArrayRef<Int_t>+;
#pragma link C++ class TNDArrayRef<Short_t>+;
#pragma link C++ class TNDArrayRef<Char_t>+;
#pragma link C++ class TNDArrayRef<ULong64_t>+;
#pragma link C++ class TNDArrayRef<ULong_t>+;
#pragma link C++ class TNDArrayRef<UInt_t>+;
#pragma link C++ class TNDArrayRef<UShort_t>+;
/*
#pragma link C++ class TNDArrayRef<const Float_t>+;
//#pragma link C++ class TNDArrayRef<const Float16_t>+;
#pragma link C++ class TNDArrayRef<const Double_t>+;
//#pragma link C++ class TNDArrayRef<const Double32_t>+;
#pragma link C++ class TNDArrayRef<const Long64_t>+;
#pragma link C++ class TNDArrayRef<const Long_t>+;
#pragma link C++ class TNDArrayRef<const Int_t>+;
#pragma link C++ class TNDArrayRef<const Short_t>+;
#pragma link C++ class TNDArrayRef<const Char_t>+;
#pragma link C++ class TNDArrayRef<const ULong64_t>+;
#pragma link C++ class TNDArrayRef<const ULong_t>+;
#pragma link C++ class TNDArrayRef<const UInt_t>+;
#pragma link C++ class TNDArrayRef<const UShort_t>+;
*/
#pragma link C++ class THn+;
#pragma link C++ class THnChain+;
#pragma link C++ class THnT<Float_t>+;
//#pragma link C++ class THnT<Float16_t>+;
#pragma link C++ class THnT<Double_t>+;
//#pragma link C++ class THnT<Double32_t>+;
#pragma link C++ class THnT<Long64_t>+;
#pragma link C++ class THnT<Long_t>+;
#pragma link C++ class THnT<Int_t>+;
#pragma link C++ class THnT<Short_t>+;
#pragma link C++ class THnT<Char_t>+;
#pragma link C++ class THnT<ULong64_t>+;
#pragma link C++ class THnT<ULong_t>+;
#pragma link C++ class THnT<UInt_t>+;
#pragma link C++ class THnT<UShort_t>+;
#pragma link C++ class THnSparse+;
#pragma link C++ class THnSparseT<TArrayD>+;
#pragma link C++ class THnSparseT<TArrayF>+;
#pragma link C++ class THnSparseT<TArrayL>+;
#pragma link C++ class THnSparseT<TArrayI>+;
#pragma link C++ class THnSparseT<TArrayS>+;
#pragma link C++ class THnSparseT<TArrayC>+;
#pragma link C++ class THnSparseArrayChunk+;
#pragma link C++ class THStack+;
#pragma link C++ class TLimit+;
#pragma link C++ class TLimitDataSource+;
#pragma link C++ class TConfidenceLevel+;
#pragma link C++ class TMultiGraph+;
#pragma link C++ class TMultiDimFit+;
#pragma link C++ class TPolyMarker-;
#pragma link C++ class TPrincipal+;
#pragma link C++ class TProfile-;
#pragma link C++ class TProfile2D-;
#pragma link C++ class TProfile3D+;
#pragma link C++ class TSpline-;
#pragma link C++ class TSpline5-;
#pragma link C++ class TSpline3-;
#pragma link C++ class TSplinePoly+;
#pragma link C++ class TSplinePoly5+;
#pragma link C++ class TSplinePoly3+;
#pragma link C++ class TVirtualHistPainter+;
#pragma link C++ class TVirtualGraphPainter+;
#pragma link C++ class TVirtualFitter+;
#pragma link C++ class TVirtualPaveStats+;
#pragma link C++ class TBackCompFitter+;
#pragma link C++ class TSVDUnfold+;
#pragma link C++ class TEfficiency+;
#pragma link C++ class TKDE+;


#pragma link C++ typedef THnSparseD;
#pragma link C++ typedef THnSparseF;
#pragma link C++ typedef THnSparseL;
#pragma link C++ typedef THnSparseI;
#pragma link C++ typedef THnSparseS;
#pragma link C++ typedef THnSparseC;

#pragma link C++ typedef THnD;
#pragma link C++ typedef THnF;
#pragma link C++ typedef THnL;
#pragma link C++ typedef THnI;
#pragma link C++ typedef THnS;
#pragma link C++ typedef THnC;


// for autoloading of typedef's (make some dummy ifdef)
// which are not recognized by the autoloading
#ifdef DO_AUTOLOAD_TYPEDEF
#pragma link C++ class THnSparseD;
#pragma link C++ class THnSparseF;
#pragma link C++ class THnSparseL;
#pragma link C++ class THnSparseI;
#pragma link C++ class THnSparseS;
#pragma link C++ class THnSparseC;

#pragma link C++ class THnD;
#pragma link C++ class THnF;
#pragma link C++ class THnL;
#pragma link C++ class THnI;
#pragma link C++ class THnS;
#pragma link C++ class THnC;
#endif


#pragma link C++ function operator*(Float_t,TH1C&);
#pragma link C++ function operator*(TH1C&, Float_t);
#pragma link C++ function operator+(TH1C&, TH1C&);
#pragma link C++ function operator-(TH1C&, TH1C&);
#pragma link C++ function operator*(TH1C&, TH1C&);
#pragma link C++ function operator/(TH1C&, TH1C&);

#pragma link C++ function operator*(Float_t,TH1S&);
#pragma link C++ function operator*(TH1S&, Float_t);
#pragma link C++ function operator+(TH1S&, TH1S&);
#pragma link C++ function operator-(TH1S&, TH1S&);
#pragma link C++ function operator*(TH1S&, TH1S&);
#pragma link C++ function operator/(TH1S&, TH1S&);

#pragma link C++ function operator*(Float_t,TH1I&);
#pragma link C++ function operator*(TH1I&, Float_t);
#pragma link C++ function operator+(TH1I&, TH1I&);
#pragma link C++ function operator-(TH1I&, TH1I&);
#pragma link C++ function operator*(TH1I&, TH1I&);
#pragma link C++ function operator/(TH1I&, TH1I&);

#pragma link C++ function operator*(Float_t,TH1F&);
#pragma link C++ function operator*(TH1F&, Float_t);
#pragma link C++ function operator+(TH1F&, TH1F&);
#pragma link C++ function operator-(TH1F&, TH1F&);
#pragma link C++ function operator*(TH1F&, TH1F&);
#pragma link C++ function operator/(TH1F&, TH1F&);

#pragma link C++ function operator*(Float_t,TH1D&);
#pragma link C++ function operator*(TH1D&, Float_t);
#pragma link C++ function operator+(TH1D&, TH1D&);
#pragma link C++ function operator-(TH1D&, TH1D&);
#pragma link C++ function operator*(TH1D&, TH1D&);
#pragma link C++ function operator/(TH1D&, TH1D&);

#pragma link C++ function operator*(Float_t,TH2C&);
#pragma link C++ function operator*(TH2C&, Float_t);
#pragma link C++ function operator+(TH2C&, TH2C&);
#pragma link C++ function operator-(TH2C&, TH2C&);
#pragma link C++ function operator*(TH2C&, TH2C&);
#pragma link C++ function operator/(TH2C&, TH2C&);

#pragma link C++ function operator*(Float_t,TH2S&);
#pragma link C++ function operator*(TH2S&, Float_t);
#pragma link C++ function operator+(TH2S&, TH2S&);
#pragma link C++ function operator-(TH2S&, TH2S&);
#pragma link C++ function operator*(TH2S&, TH2S&);
#pragma link C++ function operator/(TH2S&, TH2S&);

#pragma link C++ function operator*(Float_t,TH2I&);
#pragma link C++ function operator*(TH2I&, Float_t);
#pragma link C++ function operator+(TH2I&, TH2I&);
#pragma link C++ function operator-(TH2I&, TH2I&);
#pragma link C++ function operator*(TH2I&, TH2I&);
#pragma link C++ function operator/(TH2I&, TH2I&);

#pragma link C++ function operator*(Float_t,TH2F&);
#pragma link C++ function operator*(TH2F&, Float_t);
#pragma link C++ function operator+(TH2F&, TH2F&);
#pragma link C++ function operator-(TH2F&, TH2F&);
#pragma link C++ function operator*(TH2F&, TH2F&);
#pragma link C++ function operator/(TH2F&, TH2F&);

#pragma link C++ function operator*(Float_t,TH2D&);
#pragma link C++ function operator*(TH2D&, Float_t);
#pragma link C++ function operator+(TH2D&, TH2D&);
#pragma link C++ function operator-(TH2D&, TH2D&);
#pragma link C++ function operator*(TH2D&, TH2D&);
#pragma link C++ function operator/(TH2D&, TH2D&);

#pragma link C++ function operator*(Float_t,TH3C&);
#pragma link C++ function operator*(TH3C&, Float_t);
#pragma link C++ function operator+(TH3C&, TH3C&);
#pragma link C++ function operator-(TH3C&, TH3C&);
#pragma link C++ function operator*(TH3C&, TH3C&);
#pragma link C++ function operator/(TH3C&, TH3C&);

#pragma link C++ function operator*(Float_t,TH3S&);
#pragma link C++ function operator*(TH3S&, Float_t);
#pragma link C++ function operator+(TH3S&, TH3S&);
#pragma link C++ function operator-(TH3S&, TH3S&);
#pragma link C++ function operator*(TH3S&, TH3S&);
#pragma link C++ function operator/(TH3S&, TH3S&);

#pragma link C++ function operator*(Float_t,TH3I&);
#pragma link C++ function operator*(TH3I&, Float_t);
#pragma link C++ function operator+(TH3I&, TH3I&);
#pragma link C++ function operator-(TH3I&, TH3I&);
#pragma link C++ function operator*(TH3I&, TH3I&);
#pragma link C++ function operator/(TH3I&, TH3I&);

#pragma link C++ function operator*(Float_t,TH3F&);
#pragma link C++ function operator*(TH3F&, Float_t);
#pragma link C++ function operator+(TH3F&, TH3F&);
#pragma link C++ function operator-(TH3F&, TH3F&);
#pragma link C++ function operator*(TH3F&, TH3F&);
#pragma link C++ function operator/(TH3F&, TH3F&);

#pragma link C++ function operator*(Float_t,TH3D&);
#pragma link C++ function operator*(TH3D&, Float_t);
#pragma link C++ function operator+(TH3D&, TH3D&);
#pragma link C++ function operator-(TH3D&, TH3D&);
#pragma link C++ function operator*(TH3D&, TH3D&);
#pragma link C++ function operator/(TH3D&, TH3D&);

#pragma link C++ function R__H(Int_t);
#pragma link C++ function R__H(const char*);

#pragma link C++ class ROOT::Internal::THnBaseBrowsable;
#pragma link C++ class ROOT::Math::WrappedTF1;
#pragma link C++ class ROOT::Math::WrappedMultiTF1;

#pragma link C++ namespace ROOT::Fit;
#pragma link C++ function ROOT::Fit::FillData(ROOT::Fit::BinData &, const TH1 *, TF1 * );
#pragma link C++ function ROOT::Fit::FillData(ROOT::Fit::BinData &, const TGraph2D *, TF1 * );

#pragma link C++ namespace ROOT::Fit;
#pragma link C++ function ROOT::Fit::FillData(ROOT::Fit::BinData &, const TGraph *,  TF1 * );
#pragma link C++ function ROOT::Fit::FillData(ROOT::Fit::BinData &, const TMultiGraph *,  TF1 * );

#pragma link C++ function ROOT::Fit::FitResult::GetCovarianceMatrix<TMatrixDSym>( TMatrixDSym & );
#pragma link C++ function ROOT::Fit::FitResult::GetCorrelationMatrix<TMatrixDSym>( TMatrixDSym & );

// for having backward compatibility after new data member in TProfile
#pragma read sourceClass="TProfile" version="[1-5]" targetClass="TProfile" \
  source="" target="fBinSumw2" \
  code="{ fBinSumw2.Reset(); }"
#pragma read sourceClass="TProfile2D" version="[1-6]" targetClass="TProfile2D" \
  source="" target="fBinSumw2" \
  code="{ fBinSumw2.Reset(); }"
#pragma read sourceClass="TProfile3D" version="[1-6]" targetClass="TProfile3D" \
  source="" target="fBinSumw2" \
  code="{ fBinSumw2.Reset(); }"



#endif
