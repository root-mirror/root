// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DataInputHandler                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
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

/*! \class TMVA::DataInputHandler
\ingroup TMVA

Class that contains all the data information.

*/

#include "TMVA/DataInputHandler.h"

#include "TMVA/DataLoader.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTree.h"

#include "TMVA/Configurable.h"

#include <vector>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::DataInputHandler::DataInputHandler()
   : fLogger( new MsgLogger("DataInputHandler", kINFO) )
{
   fExplicitTrainTest["Signal"] = fExplicitTrainTest["Background"] = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::DataInputHandler::~DataInputHandler()
{
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// add a *className* tree to the dataset to be used as input

void TMVA::DataInputHandler::AddTree( const TString& fn,
                                      const TString& className,
                                      Double_t weight,
                                      const TCut& cut,
                                      Types::ETreeType tt  )
{
   TTree * tr = ReadInputTree(fn);
   tr->SetName( TString("Tree")+className );
   AddTree( tr, className, weight, cut, tt );
}

////////////////////////////////////////////////////////////////////////////////
/// add  tree of *className* events  for tt (Training;Testing..) type as input ..

void TMVA::DataInputHandler::AddTree( TTree* tree,
                                      const TString& className,
                                      Double_t weight,
                                      const TCut& cut,
                                      Types::ETreeType tt )
{
   if (!tree) Log() << kFATAL << "Zero pointer for tree of class " << className.Data() << Endl;
   if (tree->GetEntries()==0) Log() << kFATAL << "Encountered empty TTree or TChain of class " << className.Data() << Endl;
   if (fInputTrees[className.Data()].empty()) {
      // on the first tree (of the class) check if explicit treetype is given
      fExplicitTrainTest[className.Data()] = (tt != Types::kMaxTreeType);
   }
   else {
      // if the first tree has a specific type, all later tree's must also have one
      if (fExplicitTrainTest[className.Data()] != (tt!=Types::kMaxTreeType)) {
         if (tt==Types::kMaxTreeType)
            Log() << kFATAL << "For the tree " << tree->GetName() << " of class " << className.Data()
                  << " you did "<< (tt==Types::kMaxTreeType?"not ":"") << "specify a type,"
                  << " while you did "<< (tt==Types::kMaxTreeType?"":"not ") << "for the first tree "
                  << fInputTrees[className.Data()][0].GetTree()->GetName() << " of class " << className.Data()
                  << Endl;
      }
   }
   if (cut.GetTitle()[0] != 0) {
      fInputTrees[className.Data()].push_back(TreeInfo( tree->CopyTree(cut.GetTitle()), className, weight, tt ));
   }
   else {
      fInputTrees[className.Data()].push_back(TreeInfo( tree, className, weight, tt ));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// add a signal tree to the dataset to be used as input

void TMVA::DataInputHandler::AddSignalTree( TTree* tr, Double_t weight, Types::ETreeType tt )
{
   AddTree( tr, "Signal", weight, "", tt );
}

////////////////////////////////////////////////////////////////////////////////
/// add a background tree to the dataset to be used as input

void TMVA::DataInputHandler::AddBackgroundTree( TTree* tr, Double_t weight, Types::ETreeType tt )
{
   AddTree( tr, "Background", weight, "", tt );
}

////////////////////////////////////////////////////////////////////////////////
/// add a signal tree to the dataset to be used as input

void TMVA::DataInputHandler::AddSignalTree( const TString& fn, Double_t weight, Types::ETreeType tt )
{
   TTree * tr = ReadInputTree(fn);
   tr->SetName("TreeS");
   AddTree( tr, "Signal", weight, "", tt );
}

////////////////////////////////////////////////////////////////////////////////
/// add a background tree to the dataset to be used as input

void TMVA::DataInputHandler::AddBackgroundTree( const TString& fn, Double_t weight, Types::ETreeType tt )
{
   TTree * tr = ReadInputTree(fn);
   tr->SetName("TreeB");
   AddTree( tr, "Background", weight, "", tt );
}

////////////////////////////////////////////////////////////////////////////////
/// create trees from these ascii files

TTree* TMVA::DataInputHandler::ReadInputTree( const TString& dataFile )
{
   TTree* tr = new TTree( "tmp", dataFile );
   std::ifstream in(dataFile);
   tr->SetDirectory(0); Log() << kWARNING << "Watch out, I (Helge) made the Tree not associated to the current directory .. Hopefully that does not have unwanted consequences" << Endl;
   if (!in.good()) Log() << kFATAL << "Could not open file: " << dataFile << Endl;
   in.close();

   tr->ReadFile( dataFile );

   return tr;
}

////////////////////////////////////////////////////////////////////////////////
/// define the input trees for signal and background from single input tree,
/// containing both signal and background events distinguished by the type
/// identifiers: SigCut and BgCut

void TMVA::DataInputHandler::AddInputTrees(TTree* inputTree, const TCut& SigCut, const TCut& BgCut)
{
   if (!inputTree) Log() << kFATAL << "Zero pointer for input tree: " << inputTree << Endl;

   AddTree( inputTree, "Signal",     1.0, SigCut );
   AddTree( inputTree, "Background", 1.0, BgCut  );
}


////////////////////////////////////////////////////////////////////////////////

void TMVA::DataInputHandler::ClearTreeList( const TString& className )
{
   try {
      fInputTrees.find(className)->second.clear();
   }
   catch(int) {
      Log() << kINFO << "   Clear treelist for class " << className << " failed, since class does not exist." << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

std::vector< TString >* TMVA::DataInputHandler::GetClassList() const
{
   std::vector< TString >* ret = new std::vector< TString >();
   for ( std::map< TString, std::vector<TreeInfo> >::iterator it = fInputTrees.begin(); it != fInputTrees.end(); ++it ){
      ret->push_back( it->first );
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// return number of entries in tree

UInt_t TMVA::DataInputHandler::GetEntries(const std::vector<TreeInfo>& tiV) const
{
   UInt_t entries = 0;
   std::vector<TreeInfo>::const_iterator tiIt = tiV.begin();
   for (;tiIt != tiV.end();++tiIt) entries += tiIt->GetEntries();
   return entries;
}

////////////////////////////////////////////////////////////////////////////////
/// return number of entries in tree

UInt_t TMVA::DataInputHandler::GetEntries() const
{
   UInt_t number = 0;
   for (std::map< TString, std::vector<TreeInfo> >::iterator it = fInputTrees.begin(); it != fInputTrees.end(); ++it) {
      number += GetEntries( it->second );
   }
   return number;
}
