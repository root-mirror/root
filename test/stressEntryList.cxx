
// Author: Anna Kreshuk, March 2007

/////////////////////////////////////////////////////////////////
//
//___A stress test for the TEntryList class and operations with it___
//
//   The functions below test different properties of TEntryList
//   - Test1() - assembling entry lists for smaller chains from the lists
//               for bigger chains and vice versa + using them in TTree::Draw
//   - Test2() - adding and subtracting entry lists in different order
//               and using ">>+elist" in TTree::Draw
//   - Test3() - transforming TEventList objects into TEntryList objects for a TChain
//   - Test4() - same as Test3() but for a TTree
//
//   To run in batch mode, do
//     stressEntryList
//     stressEntryList 1000
//     stressEntryList 1000 10
//   Here the 1st parameter is the number of entries in each TTree,
//            2nd parameter is the number of created files
//   Default values are 10000 10
//
//   An example of output when all tests pass:
// **********************************************************************
// ***************Starting TEntryList stress test************************
// **********************************************************************
// **********Generating 10 data files, 2 trees of 10000 in each**********
// **********************************************************************
// Test1: Applying different entry lists to different chains --------- OK
// Test2: Adding and subtracting entry lists-------------------------- OK
// Test3: TEntryList and TEventList for TChain------------------------ OK
// Test4: TEntryList and TEventList for TTree------------------------- OK
// **********************************************************************
// *******************Deleting the data files****************************
// **********************************************************************

#include <map>
#include <list>
#include <array>
#include <functional>
#include <stdlib.h>
#include "TApplication.h"
#include "TEntryList.h"
#include "TEventList.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCut.h"
#include "TFile.h"
#include "TSystem.h"
#include "TError.h"
#include "snprintf.h"

Int_t stressEntryList(Int_t nentries = 10000, Int_t nfiles = 10);
void MakeTrees(Int_t nentries, Int_t nfiles);

Bool_t Test1(bool fixedCut)
{
   //Test the functionality of entry lists for chains:
   //making new entry lists out of parts of other entry lists
   //applying same entry lists to different chains, etc

   Int_t wrongentries1, wrongentries2, wrongentries3, wrongentries4, wrongentries5;
   TChain *bigchain = new TChain("bigchain", "bigchain");
   bigchain->Add("stressEntryListTrees*.root/tree1");
   bigchain->Add("stressEntryListTrees*.root/tree2");

   TChain *smallchain = new TChain("smallchain", "smallchain");
   smallchain->Add("stressEntryListTrees*.root/tree1");

   //create an entry list for the small chain
   TString fixedCutStr; fixedCutStr.Form("Entry$ >= %lld", smallchain->GetEntries());
   TCut cut = fixedCut ? fixedCutStr.Data() : "x<0 && y>0";
   smallchain->Draw(">>elist_small", cut, "entrylist");
   TEntryList *elist_small = (TEntryList*)gDirectory->Get("elist_small");

   if (fixedCut && elist_small->GetN() != 0) {
      printf("Test1: Cut \"Entry$ >= %lld\" found entries in the small chain\n", smallchain->GetEntries());
      return false;
   } else if (!fixedCut && elist_small->GetN() == 0) {
      printf("Test1: Cut \"x<0 && y>0\" found no entries in the small chain\n");
      return false;
   }

   //check if the entry list contains correct entries
   Int_t range = 100;
   TH1F *hx = new TH1F("hx", "hx", range, -range, range);
   smallchain->Draw("x >> hx", cut, "goff");
   TH1F *hcheck = new TH1F("hcheck", "hcheck", range, -range, range);

   smallchain->SetEntryList(elist_small);
   smallchain->Draw("x >> hcheck", "", "goff");
   wrongentries1 = 0;
   for (Int_t i=1; i<=range; i++){
      if (TMath::Abs(hx->GetBinContent(i)-hcheck->GetBinContent(i)) > 0.1){
         wrongentries1++;
      }
   }
   if (wrongentries1 >0)
      printf("\nsmall list and small chain: number of wrong bins=%d\n", wrongentries1);

   //set this small entry list to the big chain and check the results
   bigchain->SetEntryList(elist_small);
   bigchain->Draw("x >> hcheck_", "", "goff");
   wrongentries2 = 0;
   for (Int_t i=1; i<=range; i++){
      if (TMath::Abs(hx->GetBinContent(i)-hcheck->GetBinContent(i)) > 0.1){
         wrongentries2++;
      }
   }
   if (wrongentries2 >0)
      printf("\nsmall elist and big chain: number of wrong bins=%d\n", wrongentries2);

   smallchain->SetEntryList(0);
   bigchain->SetEntryList(0);

   //make an entry list for a big chain
   bigchain->Draw(">>elist_big", cut, "entrylist");
   TEntryList* elist_big = (TEntryList*)gDirectory->Get("elist_big");

   if (fixedCut && elist_big->GetN() != smallchain->GetEntries()) {
      printf("Test1: Cut \"Entry$ >= %lld\" did not find the right number of entries in the big chain (expected %lld got %lld\n",
            smallchain->GetEntries(), smallchain->GetEntries(), elist_big->GetN());
      return false;
   } else if (!fixedCut && elist_big->GetN() == 0) {
      printf("Test1: Cut \"x<0 && y>0\" found no entries in the big chain\n");
      return false;
   }

   //make a small entry list by extracting the lists, corresponding to the trees in
   //the small chain, from the big entry list
   TEntryList *list_extracted = new TEntryList("list_extracted", "list_extracted");
   TEntryList *elist_temp;
   TList *lists = elist_big->GetLists();
   TIter next(lists);
   while ((elist_temp = (TEntryList*)next())){
      if (!strcmp(elist_temp->GetTreeName(),"tree1"))
         list_extracted->Add(elist_temp);
   }

   //compare this extracted list to the list, generated by smallchain->Draw()
   Long64_t entry1, entry2;
   Int_t n=list_extracted->GetN();
   wrongentries3 = 0;
   for (Int_t i=0; i<n; i++){
      entry1 = list_extracted->GetEntry(i);
      entry2 = elist_small->GetEntry(i);
      if (entry1 != entry2){
         if (wrongentries3<10) printf("wrong entry: %d list2=%lld elist_small=%lld\n", i, entry1, entry2);
         wrongentries3++;
      }
   }
   if (wrongentries3 >0)
      printf("\nsmall list and extracted list: number of wrong entries = %d, n=%d\n", wrongentries3,n);

   //add another entry list to the extracted list
   elist_temp = (TEntryList*)lists->Last();
   list_extracted->Add(elist_temp);
   smallchain->SetEntryList(list_extracted);
   smallchain->Draw("x>>hcheck", "", "goff");
   wrongentries4 = 0;
   for (Int_t i=1; i<=range; i++){
      if (TMath::Abs(hx->GetBinContent(i)-hcheck->GetBinContent(i)) > 0.1){
         //printf("%d hx: %f hcheck %f\n", i, hx->GetBinContent(i), hcheck->GetBinContent(i));
         wrongentries4++;
      }
   }
   if (wrongentries4 >0)
      printf("\nextracted list with 1 wrong: number of wrong bins=%d\n", wrongentries4);

   //set the big entry list to the small chain and compare the results with
   //the entry list, generated by smallchain->Draw()
   smallchain->SetEntryList(elist_big);
   smallchain->Draw("x >> hcheck", "", "goff");
   wrongentries5 = 0;
   for (Int_t i=1; i<=range; i++){
      if (TMath::Abs(hx->GetBinContent(i)-hcheck->GetBinContent(i)) > 0.1){
         //printf("i=%d hx(i)=%f, hcheck(i)=%f\n", i, hx->GetBinContent(i), hcheck->GetBinContent(i));
         wrongentries5++;
      }
   }
   if (wrongentries5 >0)
      printf("\nbig elist and small chain: number of wrong bins = %d\n", wrongentries5);

   delete bigchain;
   delete smallchain;
   delete hx;
   delete hcheck;
   delete elist_big;
   delete elist_small;
   delete list_extracted;

   if (wrongentries1>0 || wrongentries2>0 || wrongentries3>0 || wrongentries4>0 || wrongentries5>0)
      return kFALSE;
   return kTRUE;
}

Bool_t Test1a() {
   return Test1(false);
}

Bool_t Test1b() {
   return Test1(true);
}

Bool_t Test2()
{
   //Test adding and subtracting entry lists

   Int_t wrongentries1, wrongentries2, wrongentries3, wrongentries4, wrongentries5;
   TChain *chain = new TChain("chain", "chain");
   chain->Add("stressEntryListTrees_0.root/tree1");
   chain->Add("stressEntryListTrees_0.root/tree2");
   //chain->Add("stressEntryListTrees*.root/tree1");
   //chain->Add("stressEntryListTrees*.root/tree2");
   TCut cut1("cut1", "x>0");
   TCut cut2("cut2", "y<0.1 && y>-0.1");
   TEntryList *elist1 = new TEntryList("elist1", "elist1");
   chain->Draw(">>elist1", cut1, "entrylist");
   TEntryList *elist2 = new TEntryList("elist2", "elist2");
   chain->Draw(">>elist2", cut2, "entrylist");

   //add those 2 lists (1+2)
   TEntryList *elistsum = new TEntryList("elistsum", "elistsum");
   elistsum->Add(elist1);
   elistsum->Add(elist2);

   TEntryList *elistcheck = new TEntryList("elistcheck", "elistcheck");
   chain->Draw(">>elistcheck", cut1 || cut2, "entrylist");

   Int_t n=elistcheck->GetN();
   Long64_t entry1, entry2;
   wrongentries1=0;
   for (Int_t i=0; i<n; i++){
      entry1 = elistsum->GetEntry(i);
      entry2 = elistcheck->GetEntry(i);
      if (entry1 != entry2) {
         //printf("%d, sum=%lld, check=%lld\n", i, entry1, entry2);
         wrongentries1++;
      }
   }
   if (wrongentries1>0)
      printf("\nwrong entries (1+2)=%d\n", wrongentries1);

   //add in different order
   TEntryList *elistsum2 = new TEntryList("elistsum2", "elistsum2");
   elistsum2->Add(elist2);
   elistsum2->Add(elist1);
   wrongentries2 = 0;
   for (Int_t i=0; i<n; i++){
      entry1 = elistsum2->GetEntry(i);
      entry2 = elistcheck->GetEntry(i);
      if (entry1 != entry2) {
         //printf("%d, sum=%lld, check=%lld\n", i, entry1, entry2);
         wrongentries2++;
      }
    }
    if (wrongentries2>0)
      printf("\nwrong entries (2+1)=%d\n", wrongentries2);


   //add by using "+" in TTree::Draw
   TEntryList *elistsum3 = new TEntryList("elistsum3", "elistsum3");
   chain->Draw(">>elistsum3", cut1, "entrylist");
   chain->Draw(">>+elistsum3", cut2, "entrylist");
   wrongentries3 = 0;
   for (Int_t i=0; i<n; i++){
      entry1 = elistsum3->GetEntry(i);
      entry2 = elistcheck->GetEntry(i);
      if (entry1 != entry2) {
         //printf("%d, sum=%lld, check=%lld\n", i, entry1, entry2);
         wrongentries3++;
      }
    }
   if (wrongentries3>0)
      printf("\nwrong entries with \"+\" in TChain::Draw =%d\n", wrongentries3);

   //subtract the second list
   elistsum->Subtract(elist2);
   n = elistsum->GetN();
   TEntryList *elistcheck2 = new TEntryList("elistcheck2","elistcheck2");
   chain->Draw(">>elistcheck2", cut1 && !cut2, "entrylist");

   wrongentries4 = 0;
   for (Int_t i=0; i<n; i++){
      entry1 = elistsum->GetEntry(i);
      entry2 = elistcheck2->GetEntry(i);
      if (entry1 != entry2){
          //printf("%d elist1=%lld elistsum=%lld\n", i, entry1, entry2);
         wrongentries4++;
      }
   }
   if (wrongentries4>0)
      printf("\nwrong entries after subtract 2 = %d\n", wrongentries4);

   //subtract the first list
   elistsum2->Subtract(elist1);
   elistcheck2->Reset();
   chain->Draw(">>elistcheck2", !cut1 && cut2, "entrylist");
   wrongentries5 = 0;
   n = elistcheck2->GetN();
   for (Int_t i=0; i<n; i++){
      entry1 = elistsum2->GetEntry(i);
      entry2 = elistcheck2->GetEntry(i);
      if (entry1 != entry2){
         //printf("%d elist1=%lld elistsum=%lld\n", i, entry1, entry2);
         wrongentries5++;
      }
   }
   if (wrongentries5>0)
      printf("\nwrong entries after subtract 1 = %d\n", wrongentries5);

   delete elist1;
   delete elist2;
   delete elistsum;
   delete elistsum2;
   delete elistsum3;
   delete elistcheck;
   delete elistcheck2;

   if (wrongentries1>0 || wrongentries2>0 || wrongentries3>0 || wrongentries4>0 || wrongentries5>0)
      return kFALSE;
   return kTRUE;
}

Bool_t Test3()
{
   //Test correspondence of event lists and entry lists

   TChain *chain = new TChain("chain", "chain");
   chain->Add("stressEntryListTrees*.root/tree1");
   chain->Add("stressEntryListTrees*.root/tree2");

   TCut cut = "x<0 && y>0";

   chain->Draw(">>evlist", cut, "");
   TEventList *evlist = (TEventList*)gDirectory->Get("evlist");
   chain->Draw("x>>h1", cut, "goff");
   TH1F *h1 = (TH1F*)gDirectory->Get("h1");
   chain->SetEventList(evlist);
   chain->Draw("x>>h2", "", "goff");
   TH1F *h2 = (TH1F*)gDirectory->Get("h2");

   chain->SetEventList(0);
   chain->Draw(">>enlist", cut, "entrylist");
   TEntryList *enlist = (TEntryList*)gDirectory->Get("enlist");
   chain->SetEntryList(enlist);

   chain->Draw("x>>h3", "", "goff");
   TH1F *h3 = (TH1F*)gDirectory->Get("h3");

   Int_t wrongbins = 0;
   Int_t nbins1 = h1->GetNbinsX();

   Double_t bin1,bin2,bin3;
   for (Int_t i=0; i<nbins1; i++){
      bin1 = h1->GetBinContent(i);
      bin2 = h2->GetBinContent(i);
      bin3 = h3->GetBinContent(i);
      if (TMath::Abs(bin1-bin2) > 0.1 || TMath::Abs(bin1-bin3) || TMath::Abs(bin2-bin3) > 0.1) {
         //printf("bin1=%f, bin2=%f, bin3=%f\n", bin1, bin2, bin3);
         wrongbins++;
      }
   }
   if (wrongbins>0)
      printf("wrongbins=%d\n", wrongbins);

   delete chain;
   delete h1;
   delete h2;
   delete h3;
   delete evlist;
   delete enlist;
   if (wrongbins>0)
      return kFALSE;
   return kTRUE;
}

Bool_t Test4()
{
   //Like Test3() but for trees

   TFile f("stressEntryListTrees_0.root");
   TTree *tree = (TTree*)f.Get("tree1");
   TCut cut = "x<0 && y>0";
   tree->Draw(">>evlist", cut, "");
   TEventList *evlist = (TEventList*)gDirectory->Get("evlist");
   tree->Draw("x>>h1", cut, "goff");
   TH1F *h1 = (TH1F*)gDirectory->Get("h1");
   tree->SetEventList(evlist);
   tree->Draw("x>>h2", "", "goff");
   TH1F *h2 = (TH1F*)gDirectory->Get("h2");

   tree->SetEventList(0);
   tree->Draw(">>enlist", cut, "entrylist");
   TEntryList *enlist = (TEntryList*)gDirectory->Get("enlist");
   tree->SetEntryList(enlist);
   tree->Draw("x>>h3", "", "goff");
   TH1F *h3 = (TH1F*)gDirectory->Get("h3");
   Int_t wrongbins = 0;
   Int_t nbins1 = h1->GetNbinsX();

   Double_t bin1,bin2,bin3;
   for (Int_t i=0; i<nbins1; i++){
      bin1 = h1->GetBinContent(i);
      bin2 = h2->GetBinContent(i);
      bin3 = h3->GetBinContent(i);
      if (TMath::Abs(bin1-bin2) > 0.1 || TMath::Abs(bin1-bin3) || TMath::Abs(bin2-bin3) > 0.1) {
         //printf("bin1=%f, bin2=%f, bin3=%f\n", bin1, bin2, bin3);
         wrongbins++;
      }
   }
   if (wrongbins>0)
      printf("wrongbins=%d\n", wrongbins);

   delete h1;
   delete h2;
   delete h3;
   delete evlist;
   delete enlist;
   f.Close();

   if (wrongbins>0)
      return kFALSE;
   return kTRUE;
}

Bool_t Test5And6(const std::list<const char*>& treeNamesForChain )
{
//Test entry lists with very many or very few events
//Only makes sense to check if there are > 64000 events

   TChain *chain = new TChain("chain", "chain");
   for (auto treeName : treeNamesForChain){
      chain->Add(treeName);
   }
   Int_t wrongentries1=0;
   Int_t wrongentries2=0;
   Int_t wrongentries3=0;
   Int_t wrongentries4=0;
   Int_t wrongentries5=0;
   //Let's make a full entry list
   chain->Draw(">>elfull", "", "entrylist");
   TEntryList *elfull = (TEntryList*)gDirectory->Get("elfull");

   //printf("entries in the list: %lld\n", elfull->GetN());
   //check the contents
   Long64_t cur, real;
   Int_t ntrees = chain->GetNtrees();
   Long64_t *offset = chain->GetTreeOffset();
   //This loop will check if TEntryList::Next() is correct
   for (Int_t itree=0; itree<ntrees; itree++){
      for (Int_t i=offset[itree]; i<offset[itree+1]; i++){
         real = i-offset[itree];
         cur = elfull->GetEntry(i);
         if (TMath::Abs(real-cur)>0.1){
            //printf("real=%lld, cur=%lld\n", real, cur);
            wrongentries1++;
         }
      }
   }
   //printf("wrongentries1=%d\n", wrongentries1);
//    //This loop will check if TEntryList::GetEntry() is correct
   for (Int_t itree=0; itree<ntrees; itree++){
      for (Int_t i=offset[itree]; i<offset[itree+1]; i+=2){
         real = i-offset[itree];
         cur = elfull->GetEntry(i);
         if (TMath::Abs(real-cur)>0.1){
            //printf("real=%lld, cur=%lld\n", real, cur);
            wrongentries2++;
         }
      }
   }
   //printf("wrongentries2=%d\n", wrongentries2);


   //now let's make an empty entry list
   chain->Draw(">>elempty", "x>0 && x<0", "entrylist");
   TEntryList *elempty = (TEntryList*)gDirectory->Get("elempty");
   //just a check
   Long64_t temp = elempty->GetEntry(3);
   if (TMath::Abs(temp+1)>0.1)
   wrongentries5++;

   //Merge the almost full list with the almost empty list
   //Note, how the chain pointer is passed to the Remove() and Enter()
   //functions. This is needed, because we want to remove entry #3 in
   //the chain, not entry that has #3 in the currently active sublist

   elfull->Remove(3, chain);
   elempty->Enter(3, chain);
   elfull->Add(elempty);

   //This loop will check if TEntryList::Next() is correct
   for (Int_t itree=0; itree<ntrees-1; itree++){
      for (Int_t i=offset[itree]; i<offset[itree+1]; i++){
         real = i-offset[itree];
         cur = elfull->GetEntry(i);
         if (TMath::Abs(real-cur)>0.1){
            //printf("real=%lld, cur=%lld\n", real, cur);
            wrongentries3++;
         }
      }
   }
   //printf("wrongentries3=%d\n", wrongentries3);

   //This loop will check if TEntryList::GetEntry() is correct
   for (Int_t itree=0; itree<ntrees-1; itree++){
      for (Int_t i=offset[itree]; i<offset[itree+1]; i+=2){
         real = i-offset[itree];
         cur = elfull->GetEntry(i);
         if (TMath::Abs(real-cur)>0.1){
            //printf("real=%lld, cur=%lld\n", real, cur);
            wrongentries4++;
         }
      }
   }
   // printf("wrongentries4=%d\n", wrongentries4);

   //Now let's use the full and empty lists on the chain
   chain->SetEntryList(elfull);
   chain->Draw("x>>hx", "", "goff");
   TH1F *hx = (TH1F*)gDirectory->Get("hx");
   if (TMath::Abs(hx->GetEntries()-chain->GetEntries())>0.1){
      wrongentries5++;
      //printf("entries in chain: %lld, entries in histo: %f\n", chain->GetEntries(), hx->GetEntries());
   }

   elempty->Remove(3);
   //elempty->Print("all");
   chain->SetEntryList(elempty);
   auto oldlevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError; // Draw will print a warning about the empty histogram.
   Long64_t nen = chain->Draw("x", "", "goff");
   gErrorIgnoreLevel = oldlevel;
   if (nen!=0) wrongentries5++;

   //printf("wrongentries5=%d\n", wrongentries5);

   delete elempty;
   delete elfull;
   delete hx;
   if (wrongentries1>0 || wrongentries2>0 || wrongentries3>0 || wrongentries4>0 || wrongentries5>0)
      return kFALSE;
   else
      return kTRUE;
}

Bool_t Test5()
{
   return Test5And6({"stressEntryListTrees*.root/tree1",
                     "stressEntryListTrees*.root/tree2"});
}

Bool_t Test6()
{
   return Test5And6({"stressEntryListTrees*.root/tree1",
                     "stressEntryListTrees*.root/tree2",
                     "stressEntryListTrees*.root/Dir1/tree1",
                     "stressEntryListTrees*.root/Dir1/tree2",
                     "stressEntryListTrees*.root/Dir2/tree1",
                     "stressEntryListTrees*.root/Dir2/tree2"});
}


void SetupTree(TTree* tree, Double_t &x, Double_t &y, Double_t &z)
{
   tree->Branch("x", &x, "x/D");
   tree->Branch("y", &y, "y/D");
   tree->Branch("z", &z, "z/D");
}

const char* gRootFileNameTemplate = "stressEntryListTrees_%d.root";

void MakeTrees(Int_t nentries, Int_t nfiles)
{
   //Creates nfiles files with 2 trees of nentries each

   Double_t x=0., y=0., z=0.;
   Double_t range = nentries*0.01;

   char buffer[50];
   for (Int_t ifile=0; ifile<nfiles; ifile++){
      snprintf(buffer,50, gRootFileNameTemplate, ifile);
      TFile f1(buffer, "UPDATE");
      auto dir1 = f1.mkdir("Dir1");
      auto dir2 = f1.mkdir("Dir2");
      // Init trees
      std::array<TTree*,6> trees = {{
      new TTree("tree1", "tree1"),
      new TTree("tree2", "tree2"),
      new TTree("tree1", "tree3"),
      new TTree("tree2", "tree4"),
      new TTree("tree1", "tree4"),
      new TTree("tree2", "tree5")}};

      // Set up branches
      for (auto tree : trees) {
         SetupTree(tree,x,y,z);
      }

      // Fill trees
      Int_t treeCount=0;
      for (auto tree : trees) {
         if (treeCount == 2 || treeCount == 3) dir1->cd();
         if (treeCount == 4 || treeCount == 5) dir2->cd();
         for (Int_t i=0; i<nentries; i++){
            x = gRandom->Uniform(-range, range);
            y = gRandom->Uniform(-range, range);
            z = gRandom->Uniform(-range, range);
            tree->Fill();
         }
         tree->Write();
         if (treeCount > 1)f1.cd("/");
         treeCount++;
      }

      // Close File
      f1.Close();
   }

}

void CleanUp(Int_t nfiles)
{
   char buffer[50];
   for (Int_t i=0; i<nfiles; i++){
      snprintf(buffer,50, gRootFileNameTemplate, i);
      gSystem->Unlink(buffer);
   }
}

Int_t stressEntryList(Int_t nentries, Int_t nfiles)
{
   // Make sure files are not existing already.
   CleanUp(nfiles);

   MakeTrees(nentries, nfiles);
   printf("*************************************************************************\n");
   printf("****************Starting TEntryList stress test**************************\n");
   printf("*************************************************************************\n");
   printf("***********Generating %d data files, 2 trees of %d in each************\n", nfiles, nentries);
   printf("*************************************************************************\n");

   Int_t retval = 0;
   using fcnCharPtrPair = std::pair<std::function<bool()>,const char*>;
   std::list<fcnCharPtrPair> testDescrList = {
      {Test1a, "Test1: Applying different entry lists to different chains------------ "},
      {Test1b, "Test1: Applying different entry lists to different chains (bad cuts)- "},
      {Test2,  "Test2: Adding and subtracting entry lists---------------------------- "},
      {Test3,  "Test3: TEntryList and TEventList for TChain-------------------------- "},
      {Test4,  "Test4: TEntryList and TEventList for TTree--------------------------- "},
      {Test5,  "Test5: Full and Empty TEntryList------------------------------------- "},
      {Test6,  "Test6: Full and Empty TEntryList w/ TTrees in TDirectories----------- "}
   };

   for (auto const & testDescrPair : testDescrList) {
      auto test = testDescrPair.first;
      auto descr = testDescrPair.second;
      Bool_t testRes = test();
      retval += !testRes; // increment by one upon failure
      printf("%s %s\n", descr, testRes ? "OK" : "FAILED" );
   }

   printf("*************************************************************************\n");
   printf("********************Deleting the data files******************************\n");
   printf("*************************************************************************\n");
   CleanUp(nfiles);
   return retval;
}
//_____________________________batch only_____________________
#ifndef __CINT__

int main(int argc, char *argv[])
{
   gROOT->SetBatch();
   TApplication theApp("App", &argc, argv);
   Int_t nentries = 10000;
   Int_t nfiles = 10;
   if (argc > 1) nentries = atoi(argv[1]);
   if (argc > 2) nfiles = atoi(argv[2]);
   return stressEntryList(nentries, nfiles);
}

#endif

