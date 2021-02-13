// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DataInputHandler                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Contains all the data information                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_DataInputHandler
#define ROOT_TMVA_DataInputHandler

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DataInputHandler                                                     //
//                                                                      //
// Class that contains all the data information                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include <string>

#include "TTree.h"
#include "TCut.h"

#include "TMVA/Types.h"

namespace TMVA {

   class MsgLogger;

   class TreeInfo:public TObject {

   public:

      TreeInfo( TTree* tr, const TString& className, Double_t weight=1.0, Types::ETreeType tt = Types::kMaxTreeType, Bool_t own=kFALSE ) 
      : fTree(tr), fClassName(className), fWeight(weight), fTreeType(tt), fOwner(own) {}
      TreeInfo():fTree(0),fClassName(""),fWeight(1.0), fTreeType(Types::kMaxTreeType), fOwner(kFALSE) {}
      ~TreeInfo() { if (fOwner) delete fTree; }

      TTree*           GetTree()      const { return fTree; }
      Double_t         GetWeight()    const { return fWeight; }
      UInt_t           GetEntries()   const { if( !fTree ) return 0; else return fTree->GetEntries(); }
      Types::ETreeType GetTreeType()  const { return fTreeType; }
      const TString&   GetClassName() const { return fClassName; }

   private:

      TTree*           fTree;     // pointer to the tree
      TString          fClassName;// name of the class the tree belongs to
      Double_t         fWeight;   // weight for the tree
      Types::ETreeType fTreeType; // tree is for training/testing/both
      Bool_t           fOwner;    // true if created from file
   protected:
       ClassDef(TreeInfo,1);      
   };

   class DataInputHandler :public TObject {

   public:

      DataInputHandler();
      ~DataInputHandler();

      // setters
      void     AddSignalTree    ( TTree* tr, Double_t weight=1.0, Types::ETreeType tt = Types::kMaxTreeType );
      void     AddBackgroundTree( TTree* tr, Double_t weight=1.0, Types::ETreeType tt = Types::kMaxTreeType );
      void     AddSignalTree    ( const TString& tr, Double_t weight=1.0, Types::ETreeType tt = Types::kMaxTreeType );
      void     AddBackgroundTree( const TString& tr, Double_t weight=1.0, Types::ETreeType tt = Types::kMaxTreeType );
      void     AddInputTrees    ( TTree* inputTree, const TCut& SigCut, const TCut& BgCut);

      void     AddTree          ( TTree* tree, const TString& className, Double_t weight=1.0, 
                                  const TCut& cut = "", Types::ETreeType tt = Types::kMaxTreeType );
      void     AddTree          ( const TString& tr, const TString& className, Double_t weight=1.0, 
                                  const TCut& cut = "", Types::ETreeType tt = Types::kMaxTreeType );

      // accessors
      std::vector< TString >* GetClassList() const;

      UInt_t           GetEntries( const TString& name ) const { return GetEntries( fInputTrees[name] ); }
      UInt_t           GetNTrees ( const TString& name ) const { return fInputTrees[name].size(); }

      UInt_t           GetNSignalTrees()           const { return fInputTrees["Signal"].size(); }
      UInt_t           GetNBackgroundTrees()       const { return fInputTrees["Background"].size(); }
      UInt_t           GetSignalEntries()          const { return GetEntries(fInputTrees["Signal"]); }
      UInt_t           GetBackgroundEntries()      const { return GetEntries(fInputTrees["Background"]); }
      UInt_t           GetEntries()                const;
      const TreeInfo&  SignalTreeInfo(Int_t i)     const { return fInputTrees["Signal"][i]; }
      const TreeInfo&  BackgroundTreeInfo(Int_t i) const { return fInputTrees["Background"][i]; }

      std::vector<TreeInfo>::const_iterator begin( const TString& className ) const { return fInputTrees[className].begin(); }
      std::vector<TreeInfo>::const_iterator end( const TString& className )   const { return fInputTrees[className].end(); }
      std::vector<TreeInfo>::const_iterator Sbegin() const { return begin("Signal"); }
      std::vector<TreeInfo>::const_iterator Send()   const { return end  ("Signal"); }
      std::vector<TreeInfo>::const_iterator Bbegin() const { return begin("Background"); }
      std::vector<TreeInfo>::const_iterator Bend()   const { return end  ("Background"); }

      // reset the list of trees
      void     ClearSignalTreeList()     { ClearTreeList("Signal"); }
      void     ClearBackgroundTreeList() { ClearTreeList("Background"); }
      void     ClearTreeList( const TString& className );

   private:

      UInt_t GetEntries(const std::vector<TreeInfo>& tiV) const;

      TTree * ReadInputTree( const TString& dataFile );
      
      mutable std::map< TString, std::vector<TreeInfo> > fInputTrees;        // list of input trees per class (classname is given as first parameter in the map)
      std::map< std::string, Bool_t   >                  fExplicitTrainTest; // if set to true the user has specified training and testing data explicitly
      mutable MsgLogger*                                 fLogger;            //! message logger
      MsgLogger& Log() const { return *fLogger; }
   protected:
       ClassDef(DataInputHandler,1);
   };
}

#endif
