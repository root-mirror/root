#include "TMVA/paracoor.h"

#include "TTree.h"
#include "TLeaf.h"
#include "TParallelCoord.h"
#include "TParallelCoordVar.h"
#include "TParallelCoordRange.h"

#include <vector>
using std::vector;


// plot parallel coordinates

void TMVA::paracoor(TString dataset, TString fin , Bool_t useTMVAStyle )
{
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );  
   TTree* tree = (TTree*)file->GetDirectory(dataset.Data())->Get("TestTree");
   if(!tree) {
      cout << "--- No TestTree saved in ROOT file. Parallel coordinates will not be plotted" << endl;
      return;
   }

   // first get list of leaves in tree
   TObjArray* leafList = tree->GetListOfLeaves();
   vector<TString> vars;
   vector<TString> mvas;   
   for (Int_t iar=0; iar<leafList->GetSize(); iar++) {
      TLeaf* leaf = (TLeaf*)leafList->At(iar);
      if (leaf != 0) {
         TString leafName = leaf->GetName();
         if (leafName != "type" && leafName != "weight"  && leafName != "boostweight" &&
             leafName != "class" && leafName != "className" && leafName != "classID" && 
             !leafName.Contains("prob_")) {
            // is MVA ?
            if (TMVAGlob::ExistMethodName( leafName,file->GetDirectory(dataset.Data()) )) {
               mvas.push_back( leafName );
            }
            else {
               vars.push_back( leafName );
            }
         }
      }
   }

   cout << "--- Found: " << vars.size() << " variables" << endl;
   cout << "--- Found: " << mvas.size() << " MVA(s)" << endl;
   

   TString type[2] = { "Signal", "Background" };
   for (UInt_t imva=0; imva<mvas.size(); imva++) {
      cout << "--- Plotting parallel coordinates for : " << mvas[imva] << " & input variables" << endl;

      for (Int_t itype=0; itype<2; itype++) {

         // create draw option
         TString varstr = mvas[imva] + ":";
         for (UInt_t ivar=0; ivar<vars.size(); ivar++) varstr += vars[ivar] + ":";
         varstr.Resize( varstr.Last( ':' ) );

         // create canvas
         TString mvashort = mvas[imva]; mvashort.ReplaceAll("MVA_","");
         auto c1 = new TCanvas( Form( "c1_%i_%s",itype,mvashort.Data() ),
                           Form( "Parallel coordinate representation for %s and input variables (%s events)", 
                                 mvashort.Data(), type[itype].Data() ), 
                           50*(itype), 50*(itype), 750, 500 );      
         tree->Draw( varstr.Data(), Form("classID==%i",1-itype) , "para" );
         c1->ToggleEditor();
         gStyle->SetOptTitle(0);

         TParallelCoord*    para   = (TParallelCoord*)gPad->GetListOfPrimitives()->FindObject( "ParaCoord" );
         TParallelCoordVar* mvavar = (TParallelCoordVar*)para->GetVarList()->FindObject( mvas[imva] );
         Double_t minrange = tree->GetMinimum( mvavar->GetName() );
         Double_t maxrange = tree->GetMaximum( mvavar->GetName() );
         Double_t width    = 0.2*(maxrange - minrange);
         Double_t x1 = minrange, x2 = x1 + width;
         TParallelCoordRange* parrange = new TParallelCoordRange( mvavar, x1, x2 );
         parrange->SetLineColor(4);
         mvavar->AddRange( parrange );

         para->AddSelection("-1");

         for (Int_t ivar=1; ivar<TMath::Min(Int_t(vars.size()) + 1,3); ivar++) {
            TParallelCoordVar* var = (TParallelCoordVar*)para->GetVarList()->FindObject( vars[ivar] );
            minrange = tree->GetMinimum( var->GetName() );
            maxrange = tree->GetMaximum( var->GetName() );
            width    = 0.2*(maxrange - minrange);

            switch (ivar) {
            case 0: { x1 = minrange; x2 = x1 + width; break; }
            case 1: { x1 = 0.5*(maxrange + minrange - width)*0.02; x2 = x1 + width*0.02; break; }
            case 2: { x1 = maxrange - width; x2 = x1 + width; break; }
            }

            parrange = new TParallelCoordRange( var, x1, x2 );
            parrange->SetLineColor( ivar == 0 ? 2 : ivar == 1 ? 5 : 6 );
            var->AddRange( parrange );

            para->AddSelection( Form("%i",ivar) );
         }

         c1->Update();

         TString fname = Form( "%s/plots/paracoor_c%i_%s",dataset.Data(), imva, itype == 0 ? "S" : "B" );
         TMVAGlob::imgconv( c1, fname );
      }
   }
}

