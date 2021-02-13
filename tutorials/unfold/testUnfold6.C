/// \file
/// \ingroup tutorial_unfold
/// \notebook
/// Test program for the class TUnfoldBinning.
///
/// read a simple binning scheme and create/test bin maps
///
/// \macro_output
/// \macro_code
///
///  **Version 17.6, in parallel to changes in TUnfold**
///
/// #### History:
///  - Version 17.5, in parallel to changes in TUnfold
///  - Version 17.4, in parallel to changes in TUnfold
///  - Version 17.3, test bin map functionality in TUnfoldBinning
///
///  This file is part of TUnfold.
///
///  TUnfold is free software: you can redistribute it and/or modify
///  it under the terms of the GNU General Public License as published by
///  the Free Software Foundation, either version 3 of the License, or
///  (at your option) any later version.
///
///  TUnfold is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
///  GNU General Public License for more details.
///
///  You should have received a copy of the GNU General Public License
///  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
///
/// \author Stefan Schmitt DESY, 14.10.2008

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"

using namespace std;

void PrintBinMap(TUnfoldBinning *binning,const char * where,
                 const Int_t *binMap);

void testUnfold6()
{
  TDOMParser parser;
  ofstream dtdFile("tunfoldbinning.dtd");
  TUnfoldBinningXML::WriteDTD(dtdFile);
  dtdFile.close();
  TString dir = gSystem->UnixPathName(gSystem->GetDirName(__FILE__));
  Int_t error=parser.ParseFile(dir+"/testUnfold6binning.xml");
  if(error) cout<<"error="<<error<<" from TDOMParser\n";
  TXMLDocument const *XMLdocument=parser.GetXMLDocument();
  TUnfoldBinningXML *binning=
     TUnfoldBinningXML::ImportXML(XMLdocument,"binning");
  if(!binning) {
     cout<<"error: can not read binning (document empty?)\n";
  } else {
     cout<<"Binning scheme:\n =================================\n";
     binning->PrintStream(cout);
     Int_t *binMap = binning->CreateEmptyBinMap();
     PrintBinMap(binning,"CreateEmptyBinMap",binMap);

     TUnfoldBinning const *branch1 = binning->FindNode("branch1");
     branch1->FillBinMap1D(binMap,"y[C]",2);
     PrintBinMap(binning,"branch1->FillBinMap1D(...,\"y[C]\",...,2)",binMap);

     delete [] binMap;
     binMap = binning->CreateEmptyBinMap();
     TUnfoldBinning const *branch2=binning->FindNode("branch2");
     branch2->FillBinMap1D(binMap,"x[C]",7);
     PrintBinMap(binning,"branch2->FillBinMap1D(...,\"x[C]\",...,7)",binMap);

     delete [] binMap;
     binMap = binning->CreateEmptyBinMap();
     binning->FillBinMap1D(binMap,"y[C]",1);
     PrintBinMap(binning,"binning->FillBinMap1D(...,\"y[C]\",...,1)",binMap);

     binning->ExportXML("testUnfold6.out.xml");

     delete [] binMap;

  }
}

void PrintBinMap(TUnfoldBinning *binning,const char * where,
                 const Int_t *binMap) {

   cout<<"\n"<<where<<"\n=======================\n";
   cout<<"global bin:";
   for(int i=0;i<binning->GetEndBin()+1;i++) {
      cout<<setw(3)<<i;
   }
   cout<<"\n";
   cout<<"mapped to: ";
   for(int i=0;i<binning->GetEndBin()+1;i++) {
      cout<<setw(3)<<binMap[i];
   }
   cout<<"\n";
   map<int,vector<int> > destBin;
   for(int i=0;i<binning->GetEndBin()+1;i++) {
      destBin[binMap[i]].push_back(i);
   }
   bool printed=false;
   for(map<int,vector<int> >::const_iterator i=destBin.begin();i!=destBin.end();i++) {
      if((*i).first>=0) {
         if(!printed) {
            cout<<"\ndest |contributing bins\n"
                <<"=====+======================================\n";
            printed=true;
         }
         for(size_t j=0;j<(*i).second.size();j++) {
            cout<<setw(4)<<(*i).first<<" |";
            cout<<setw(3)<<binning->GetBinName((*i).second[j])<<"\n";
         }
         cout<<"=====+======================================\n";
      }
   }
}
