// @(#)root/roostats:$Id$
// Author: George Lewis, Kyle Cranmer
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef HISTFACTORY_SYSTEMATICS_H
#define HISTFACTORY_SYSTEMATICS_H

#include <string>
#include <fstream>
#include <iostream>

#include "TH1.h"
#include "RooStats/HistFactory/HistRef.h"

namespace RooStats{
namespace HistFactory {

  namespace Constraint {
    enum Type{ Gaussian, Poisson };
    std::string Name( Type type );
    Type GetType( const std::string& Name );
  }


  // Base class for common functions
  /*
  class Systematic {

  public:

    virtual void Print(std::ostream& = std::cout);
    virtual void writeToFile(const std::string& FileName,
			     const std::string& Directory);


  };
  */

/** \class OverallSys
 * \ingroup HistFactory
 * Configuration for a constrained overall systematic to scale sample normalisations.
 */
  class OverallSys {

  public:

    OverallSys() : fLow(0), fHigh(0) {}

    void SetName( const std::string& Name ) { fName = Name; }
    const std::string& GetName() const { return fName; }

    void SetLow( double Low )   { fLow  = Low; }
    void SetHigh( double High ) { fHigh = High; }
    double GetLow() const { return fLow; }
    double GetHigh() const { return fHigh; }

    void Print(std::ostream& = std::cout) const;
    void PrintXML(std::ostream&) const;

  protected:
    std::string fName;
    double fLow;
    double fHigh;

  };

/** \class NormFactor
 * \ingroup HistFactory
 * Configuration for an \a un- constrained overall systematic to scale sample normalisations.
 */
  class NormFactor {

  public:

    NormFactor();

    void SetName( const std::string& Name ) { fName = Name; }
    std::string GetName() const { return fName; }

    void SetVal( double Val ) { fVal = Val; }
    double GetVal() const { return fVal; }

    void SetConst( bool Const=true ) { fConst = Const; }
    bool GetConst() const { return fConst; }

    void SetLow( double Low )   { fLow  = Low; }
    void SetHigh( double High ) { fHigh = High; }
    double GetLow() const { return fLow; }
    double GetHigh() const { return fHigh; }

    void Print(std::ostream& = std::cout) const;
    void PrintXML(std::ostream&) const;

  protected:

    std::string fName;
    double fVal;
    double fLow;
    double fHigh;
    bool fConst;

  };


  /** ////////////////////////////////////////////////////////////////////////////////////////////
   * \class HistogramUncertaintyBase
   * \ingroup HistFactory
   * Base class to store the up and down variations for histogram uncertainties.
   * Use the derived classes for actual models.
   */
  class HistogramUncertaintyBase {

  public:

    HistogramUncertaintyBase() : fhLow(nullptr), fhHigh(nullptr) {}
    HistogramUncertaintyBase(const std::string& Name) : fName(Name), fhLow(nullptr), fhHigh(nullptr) {}
    HistogramUncertaintyBase(const HistogramUncertaintyBase& oth) :
    fName{oth.fName},
    fInputFileLow{oth.fInputFileLow}, fHistoNameLow{oth.fHistoNameLow}, fHistoPathLow{oth.fHistoPathLow},
    fInputFileHigh{oth.fInputFileHigh}, fHistoNameHigh{oth.fHistoNameHigh}, fHistoPathHigh{oth.fHistoPathHigh},
    fhLow{oth.fhLow ? static_cast<TH1*>(oth.fhLow->Clone()) : nullptr},
    fhHigh{oth.fhHigh ? static_cast<TH1*>(oth.fhHigh->Clone()) : nullptr} {

    }
    HistogramUncertaintyBase(HistogramUncertaintyBase&&) = default;

    virtual ~HistogramUncertaintyBase() {};


    // Need deep copies because the class owns its histograms.
    HistogramUncertaintyBase& operator=(const HistogramUncertaintyBase& oth) {
      fName = oth.fName;
      fInputFileLow = oth.fInputFileLow;
      fHistoNameLow = oth.fHistoNameLow;
      fHistoPathLow = oth.fHistoPathLow;
      fInputFileHigh = oth.fInputFileHigh;
      fHistoNameHigh = oth.fHistoNameHigh;
      fHistoPathHigh = oth.fHistoPathHigh;
      fhLow.reset(oth.fhLow ? static_cast<TH1*>(oth.fhLow->Clone()) : nullptr);
      fhHigh.reset(oth.fhHigh ? static_cast<TH1*>(oth.fhHigh->Clone()) : nullptr);

      return *this;
    }
    HistogramUncertaintyBase& operator=(HistogramUncertaintyBase&&) = default;

    virtual void Print(std::ostream& = std::cout) const;
    virtual void PrintXML(std::ostream&) const = 0;
    virtual void writeToFile( const std::string& FileName, const std::string& DirName );

    void SetHistoLow(TH1* Low ) {Low->SetDirectory(nullptr); fhLow.reset(Low);}
    void SetHistoHigh(TH1* High ) {High->SetDirectory(nullptr); fhHigh.reset(High);}

    const TH1* GetHistoLow() const {return fhLow.get();}
    const TH1* GetHistoHigh() const {return fhHigh.get();}

    void SetName( const std::string& Name ) { fName = Name; }
    const std::string& GetName() const { return fName; }

    void SetInputFileLow( const std::string& InputFileLow ) { fInputFileLow = InputFileLow; }
    void SetInputFileHigh( const std::string& InputFileHigh ) { fInputFileHigh = InputFileHigh; }

    const std::string& GetInputFileLow() const { return fInputFileLow; }
    const std::string& GetInputFileHigh() const { return fInputFileHigh; }

    void SetHistoNameLow( const std::string& HistoNameLow ) { fHistoNameLow = HistoNameLow; }
    void SetHistoNameHigh( const std::string& HistoNameHigh ) { fHistoNameHigh = HistoNameHigh; }

    const std::string& GetHistoNameLow() const { return fHistoNameLow; }
    const std::string& GetHistoNameHigh() const { return fHistoNameHigh; }

    void SetHistoPathLow( const std::string& HistoPathLow ) { fHistoPathLow = HistoPathLow; }
    void SetHistoPathHigh( const std::string& HistoPathHigh ) { fHistoPathHigh = HistoPathHigh; }

    const std::string& GetHistoPathLow() const { return fHistoPathLow; }
    const std::string& GetHistoPathHigh() const { return fHistoPathHigh; }

  protected:

    std::string fName;

    std::string fInputFileLow;
    std::string fHistoNameLow;
    std::string fHistoPathLow;

    std::string fInputFileHigh;
    std::string fHistoNameHigh;
    std::string fHistoPathHigh;

    // The Low and High Histograms
    std::unique_ptr<TH1> fhLow;
    std::unique_ptr<TH1> fhHigh;

  };

/** \class HistoSys
 * \ingroup HistFactory
 * Configuration for a constrained, coherent shape variation of affected samples.
 */
class HistoSys final : public HistogramUncertaintyBase {
public:
  virtual ~HistoSys() {}
  virtual void PrintXML(std::ostream&) const override;
};

/** \class HistoFactor
 * \ingroup HistFactory
 * Configuration for an *un*constrained, coherent shape variation of affected samples.
 */
  class HistoFactor final : public HistogramUncertaintyBase {
  public:
    virtual ~HistoFactor() {}
    void PrintXML(std::ostream&) const override;
  };

/** \class ShapeSys
 * \ingroup HistFactory
 * Constrained bin-by-bin variation of affected histogram.
 */
  class ShapeSys final : public HistogramUncertaintyBase {

  public:

    ShapeSys() :
      HistogramUncertaintyBase(),
      fConstraintType(Constraint::Gaussian) {}
    ShapeSys(const ShapeSys& other) :
      HistogramUncertaintyBase(other),
      fConstraintType(other.fConstraintType) {}
    ShapeSys& operator=(const ShapeSys& oth) {
       if (this == &oth) return *this;
       HistogramUncertaintyBase::operator=(oth);
       fConstraintType = oth.fConstraintType;
       return *this;
    }
    ShapeSys& operator=(ShapeSys&&) = default;

    void SetInputFile( const std::string& InputFile ) { fInputFileHigh = InputFile; }
    std::string GetInputFile() const { return fInputFileHigh; }

    void SetHistoName( const std::string& HistoName ) { fHistoNameHigh = HistoName; }
    std::string GetHistoName() const { return fHistoNameHigh; }

    void SetHistoPath( const std::string& HistoPath ) { fHistoPathHigh = HistoPath; }
    std::string GetHistoPath() const { return fHistoPathHigh; }

    void Print(std::ostream& = std::cout) const override;
    void PrintXML(std::ostream&) const override;
    void writeToFile( const std::string& FileName, const std::string& DirName ) override;

    const TH1* GetErrorHist() const {
      return fhHigh.get();
    }
    void SetErrorHist(TH1* hError) {
      fhHigh.reset(hError);
    }

    void SetConstraintType( Constraint::Type ConstrType ) { fConstraintType = ConstrType; }
    Constraint::Type GetConstraintType() const { return fConstraintType; }

  protected:
    Constraint::Type fConstraintType;
  };

/** \class ShapeFactor
 * \ingroup HistFactory
 * *Un*constrained bin-by-bin variation of affected histogram.
 */
  class ShapeFactor : public HistogramUncertaintyBase {

  public:

    ShapeFactor() :
      HistogramUncertaintyBase(),
      fConstant{false},
      fHasInitialShape{false} {}

    void Print(std::ostream& = std::cout) const override;
    void PrintXML(std::ostream&) const override;
    void writeToFile( const std::string& FileName, const std::string& DirName) override;

    void SetInitialShape(TH1* shape) {
      fhHigh.reset(shape);
    }
    const TH1* GetInitialShape() const { return fhHigh.get(); }

    void SetConstant(bool constant) { fConstant = constant; }
    bool IsConstant() const { return fConstant; }

    bool HasInitialShape() const { return fHasInitialShape; }

    void SetInputFile( const std::string& InputFile ) {
      fInputFileHigh = InputFile;
      fHasInitialShape=true;
    }
    const std::string& GetInputFile() const { return fInputFileHigh; }

    void SetHistoName( const std::string& HistoName ) {
      fHistoNameHigh = HistoName;
      fHasInitialShape=true;
    }
    const std::string& GetHistoName() const { return fHistoNameHigh; }

    void SetHistoPath( const std::string& HistoPath ) {
      fHistoPathHigh = HistoPath;
      fHasInitialShape=true;
    }
    const std::string& GetHistoPath() const { return fHistoPathHigh; }

  protected:

    bool fConstant;

    // A histogram representing
    // the initial shape
    bool fHasInitialShape;
  };

/** \class StatError
 * \ingroup HistFactory
 * Statistical error of Monte Carlo predictions.
 */
  class StatError : public HistogramUncertaintyBase {

  public:

    StatError() :
      HistogramUncertaintyBase(),
      fActivate(false), fUseHisto(false) {}

    void Print(std::ostream& = std::cout) const override;
    void PrintXML(std::ostream&) const override;
    void writeToFile( const std::string& FileName, const std::string& DirName ) override;

    void Activate( bool IsActive=true ) { fActivate = IsActive; }
    bool GetActivate() const { return fActivate; }

    void SetUseHisto( bool UseHisto=true ) { fUseHisto = UseHisto; }
    bool GetUseHisto() const { return fUseHisto; }

    void SetInputFile( const std::string& InputFile ) { fInputFileHigh = InputFile; }
    const std::string& GetInputFile() const { return fInputFileHigh; }

    void SetHistoName( const std::string& HistoName ) { fHistoNameHigh = HistoName; }
    const std::string& GetHistoName() const { return fHistoNameHigh; }

    void SetHistoPath( const std::string& HistoPath ) { fHistoPathHigh = HistoPath; }
    const std::string& GetHistoPath() const { return fHistoPathHigh; }


    const TH1* GetErrorHist() const {
      return fhHigh.get();
    }
    void SetErrorHist(TH1* Error) {
      fhHigh.reset(Error);
    }

  protected:

    bool fActivate;
    bool fUseHisto; // Use an external histogram for the errors
  };

/** \class StatErrorConfig
 * \ingroup HistFactory
 * Configuration to automatically assign nuisance parameters for the statistical
 * error of the Monte Carlo simulations.
 * The default is to assign a Poisson uncertainty to a bin when its statistical uncertainty
 * is larger than 5% of the bin content.
 */
  class StatErrorConfig {

  public:

    StatErrorConfig() : fRelErrorThreshold( .05 ), fConstraintType( Constraint::Poisson ) {;}
    void Print(std::ostream& = std::cout) const;
    void PrintXML(std::ostream&) const;

    void SetRelErrorThreshold( double Threshold ) { fRelErrorThreshold = Threshold; }
    double GetRelErrorThreshold() const { return fRelErrorThreshold; }

    void SetConstraintType( Constraint::Type ConstrType ) { fConstraintType = ConstrType; }
    Constraint::Type GetConstraintType() const { return fConstraintType; }

  protected:

    double fRelErrorThreshold;
    Constraint::Type fConstraintType;

  };


}
}

#endif
