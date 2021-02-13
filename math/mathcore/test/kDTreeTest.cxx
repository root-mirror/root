/*
  Test macro for TKDTree

  TestBuild();       // test build function of kdTree for memory leaks
  TestSpeed();       // test the CPU consumption to build kdTree
  TestkdtreeIF();    // test functionality of the kdTree
  TestSizeIF();      // test the size of kdtree - search application - Alice TPC tracker situation
  //
*/

//#include <malloc.h>
#include "TSystem.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TKDTree.h"
#include "TApplication.h"
#include "TVirtualPad.h"
#include <iostream>


bool showGraphics = false;


void TestBuild(const Int_t npoints = 1000000, const Int_t bsize = 100);
void TestConstr(const Int_t npoints = 1000000, const Int_t bsize = 100);
void TestSpeed(Int_t npower2 = 20, Int_t bsize = 10);

//void TestkdtreeIF(Int_t npoints=1000, Int_t bsize=9, Int_t nloop=1000, Int_t mode = 2);
//void TestSizeIF(Int_t nsec=36, Int_t nrows=159, Int_t npoints=1000,  Int_t bsize=10, Int_t mode=1);



Float_t Mem()
{
  // get mem info
  ProcInfo_t procInfo;
  gSystem->GetProcInfo(&procInfo);
  return procInfo.fMemVirtual;
}

////////////////////////////////////////////////////////////////////////////////
///
///
///

void kDTreeTest()
{
  printf("\n\tTesting kDTree memory usage ...\n");
  TestBuild();
  printf("\n\tTesting kDTree speed ...\n");
  TestSpeed();
}

////////////////////////////////////////////////////////////////////////////////
///
/// Test kdTree for memory leaks
///

void TestBuild(const Int_t npoints, const Int_t bsize){
   Float_t *data0 =  new Float_t[npoints*2];
   Float_t *data[2];
   data[0] = &data0[0];
   data[1] = &data0[npoints];
   for (Int_t i=0;i<npoints;i++) {
      data[1][i]= gRandom->Rndm();
      data[0][i]= gRandom->Rndm();
   }
   Float_t before =Mem();
   TKDTreeIF *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
   kdtree->Build();
   Float_t after = Mem();
   printf("Memory usage %f KB\n",after-before);
   delete kdtree;
   Float_t end = Mem();
   printf("Memory leak %f KB\n", end-before);
   delete[] data0;
   return;
}

////////////////////////////////////////////////////////////////////////////////
///This is not really a test, it's a function that illustrates the internal
///behaviour of the kd-tree.
///
///Print out the internal kd-tree data-members, like fCrossNode, for
///better understading

void TestMembers()
{

   TKDTreeIF *kdtree = 0x0;
   Int_t npoints = 33;
   Int_t bsize = 10;
   Float_t *data0 = new Float_t[200]; //not to reallocate each time
   Float_t *data1 = new Float_t[200];
   for (Int_t i=0;i<npoints;i++) {
      data0[i]= gRandom->Rndm();
      data1[i]= gRandom->Rndm();
   }

   kdtree = new TKDTreeIF(npoints, 2, bsize);
   kdtree->SetData(0, data0);
   kdtree->SetData(1, data1);
   kdtree->Build();

   printf("fNNodes %d, fRowT0 %d, fCrossNode %d, fOffset %d\n",kdtree->GetNNodes(), kdtree->GetRowT0(), kdtree->GetCrossNode(), kdtree->GetOffset());
   delete kdtree;
   npoints = 44;
   for (Int_t i=0;i<npoints;i++) {
      data0[i]= gRandom->Rndm();
      data1[i]= gRandom->Rndm();
   }
   kdtree = new TKDTreeIF(npoints, 2, bsize);
   kdtree->SetData(0, data0);
   kdtree->SetData(1, data1);
   kdtree->Build();

   printf("fNNodes %d, fRowT0 %d, fCrossNode %d, fOffset %d\n",kdtree->GetNNodes(), kdtree->GetRowT0(), kdtree->GetCrossNode(), kdtree->GetOffset());
   delete kdtree;
   npoints = 55;
   for (Int_t i=0;i<npoints;i++) {
      data0[i]= gRandom->Rndm();
      data1[i]= gRandom->Rndm();
   }
   kdtree = new TKDTreeIF(npoints, 2, bsize);
   kdtree->SetData(0, data0);
   kdtree->SetData(1, data1);
   kdtree->Build();

   printf("fNNodes %d, fRowT0 %d, fCrossNode %d, fOffset %d\n",kdtree->GetNNodes(), kdtree->GetRowT0(), kdtree->GetCrossNode(), kdtree->GetOffset());
   delete kdtree;
   npoints = 66;
   for (Int_t i=0;i<npoints;i++) {
      data0[i]= gRandom->Rndm();
      data1[i]= gRandom->Rndm();
   }
      kdtree = new TKDTreeIF(npoints, 2, bsize);
   kdtree->SetData(0, data0);
   kdtree->SetData(1, data1);
   kdtree->Build();

   printf("fNNodes %d, fRowT0 %d, fCrossNode %d, fOffset %d\n",kdtree->GetNNodes(), kdtree->GetRowT0(), kdtree->GetCrossNode(), kdtree->GetOffset());
   delete kdtree;
   npoints = 77;
   for (Int_t i=0;i<npoints;i++) {
      data0[i]= gRandom->Rndm();
      data1[i]= gRandom->Rndm();
   }
   kdtree = new TKDTreeIF(npoints, 2, bsize);
   kdtree->SetData(0, data0);
   kdtree->SetData(1, data1);
   kdtree->Build();

   printf("fNNodes %d, fRowT0 %d, fCrossNode %d, fOffset %d\n",kdtree->GetNNodes(), kdtree->GetRowT0(), kdtree->GetCrossNode(), kdtree->GetOffset());
   delete kdtree;
   npoints = 88;
   for (Int_t i=0;i<npoints;i++) {
      data0[i]= gRandom->Rndm();
      data1[i]= gRandom->Rndm();
   }
   kdtree = new TKDTreeIF(npoints, 2, bsize);
   kdtree->SetData(0, data0);
   kdtree->SetData(1, data1);
   kdtree->Build();

   printf("fNNodes %d, fRowT0 %d, fCrossNode %d, fOffset %d\n",kdtree->GetNNodes(), kdtree->GetRowT0(), kdtree->GetCrossNode(), kdtree->GetOffset());
   delete kdtree;



   delete[] data0;
   delete[] data1;
}



////////////////////////////////////////////////////////////////////////////////
///
///compare the results of different data setting functions
///nothing printed - all works correctly

void TestConstr(const Int_t npoints, const Int_t bsize)
{
   Float_t *data0 =  new Float_t[npoints*2];
   Float_t *data[2];
   data[0] = &data0[0];
   data[1] = &data0[npoints];
   for (Int_t i=0;i<npoints;i++) {
      data[1][i]= gRandom->Rndm();
      data[0][i]= gRandom->Rndm();
   }
   Float_t before =Mem();
   TKDTreeIF *kdtree1 = new TKDTreeIF(npoints, 2, bsize, data);
   kdtree1->Build();
   TKDTreeIF *kdtree2 = new TKDTreeIF(npoints, 2, bsize);
   kdtree2->SetData(0, data[0]);
   kdtree2->SetData(1, data[1]);
   kdtree2->Build();
   Int_t nnodes = kdtree1->GetNNodes();
   if (nnodes - kdtree2->GetNNodes()>1){
      printf("different number of nodes\n");
      return;
   }
   for (Int_t inode=0; inode<nnodes; inode++){
      Float_t value1 = kdtree1->GetNodeValue(inode);
      Float_t value2 = kdtree2->GetNodeValue(inode);
      if (TMath::Abs(value1-value2 > 0.001)){
         printf("node %d value: %f %f\n", inode, kdtree1->GetNodeValue(inode), kdtree2->GetNodeValue(inode));
      }
   }
   delete kdtree1;
   delete kdtree2;
   Float_t end = Mem();
   printf("Memory leak %f KB\n", end-before);
   delete[] data0;
   return;
}


////////////////////////////////////////////////////////////////////////////////
///
/// Test of building time of kdTree
///

void TestSpeed(Int_t npower2, Int_t bsize)
{
  if(npower2 < 10){
    printf("Please specify a power of 2 greater than 10\n");
    return;
  }

  Int_t npoints = Int_t(pow(2., npower2))*bsize;
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  for (Int_t i=0;i<npoints;i++) {
    data[1][i]= gRandom->Rndm();
    data[0][i]= gRandom->Rndm();
  }

  TGraph *g = new TGraph(npower2-10);
  g->SetMarkerStyle(7);
  TStopwatch timer;
  Int_t tpoints;
  TKDTreeIF *kdtree = 0x0;
  for(int i=10; i<npower2; i++){
    tpoints = Int_t(pow(2., i))*bsize;
    timer.Start(kTRUE);
    kdtree = new TKDTreeIF(tpoints, 2, bsize, data);
    kdtree->Build();
    timer.Stop();
    g->SetPoint(i-10, i, timer.CpuTime());
    printf("npoints [%d] nodes [%d] cpu time %f [s]\n", tpoints, kdtree->GetNNodes(), timer.CpuTime());
    //timer.Print("u");
    delete kdtree;
  }
  if (showGraphics) {
     g->Draw("apl");
     gPad->Update();
  }

  delete[] data0;
  return;
}

/*
////////////////////////////////////////////////////////////////////////////////
///
/// Test size to build kdtree
///

void TestSizeIF(Int_t nsec, Int_t nrows, Int_t npoints,  Int_t bsize, Int_t mode)
{
  Float_t before =Mem();
  for (Int_t isec=0; isec<nsec;isec++)
    for (Int_t irow=0;irow<nrows;irow++){
      TestkdtreeIF(npoints,1,mode,bsize);
    }
  Float_t after = Mem();
  printf("Memory usage %f\n",after-before);
}
*/

/*
////////////////////////////////////////////////////////////////////////////////
///
/// Test speed and functionality of 2D kdtree.
/// Input parametrs:
/// npoints - number of data points
/// bsize   - bucket size
/// nloop   - number of loops
/// mode    - tasks to be performed by the kdTree
///         - 0  : time building the tree
///

void  TestkdtreeIF(Int_t npoints, Int_t bsize, Int_t nloop, Int_t mode)
{

  Float_t rangey  = 100;
  Float_t rangez  = 100;
  Float_t drangey = 0.1;
  Float_t drangez = 0.1;

  //
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  //Int_t i;
  for (Int_t i=0; i<npoints; i++){
    data[0][i]          = gRandom->Uniform(-rangey, rangey);
    data[1][i]          = gRandom->Uniform(-rangez, rangez);
  }
  TStopwatch timer;

  // check time build
  printf("building kdTree ...\n");
  timer.Start(kTRUE);
  TKDTreeIF *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
  kdtree->Build();
  timer.Stop();
  timer.Print();
  if(mode == 0) return;

  Float_t countern=0;
  Float_t counteriter  = 0;
  Float_t counterfound = 0;

  if (mode ==2){
    if (nloop) timer.Start(kTRUE);
    Int_t *res = new Int_t[npoints];
    Int_t nfound = 0;
    for (Int_t kloop = 0;kloop<nloop;kloop++){
      if (kloop==0){
   counteriter = 0;
   counterfound= 0;
   countern    = 0;
      }
      for (Int_t i=0;i<npoints;i++){
   Float_t point[2]={data[0][i],data[1][i]};
   Float_t delta[2]={drangey,drangez};
   Int_t iter  =0;
   nfound =0;
   Int_t bnode =0;
   //kdtree->FindBNode(point,delta, bnode);
   //continue;
   kdtree->FindInRangeA(point,delta,res,nfound,iter,bnode);
   if (kloop==0){
     //Bool_t isOK = kTRUE;
     Bool_t isOK = kFALSE;
     for (Int_t ipoint=0;ipoint<nfound;ipoint++)
       if (res[ipoint]==i) isOK =kTRUE;
     counteriter+=iter;
     counterfound+=nfound;
     if (isOK) {
       countern++;
     }else{
       printf("Bug\n");
     }
   }
      }
    }

    if (nloop){
      timer.Stop();
      timer.Print();
    }

    delete [] res;
  }
  delete [] data0;

  counteriter/=npoints;
  counterfound/=npoints;
  if (nloop) printf("Find nearest point:\t%f\t%f\t%f\n",countern, counteriter, counterfound);
}
*/

////////////////////////////////////////////////////////////////////////////////
///Test TKDTree::FindNearestNeighbors() function

void TestNeighbors()
{
//Generate some 3d points
   Int_t npoints = 10000;
   Int_t nn = 100;
   Int_t ntimes = 100;
   Int_t bsize = 10; //bucket size of the kd-tree

   Double_t *x = new Double_t[npoints];
   Double_t *y = new Double_t[npoints];
   Double_t *z = new Double_t[npoints];
   for (Int_t i=0; i<npoints; i++){
      x[i] = gRandom->Uniform(-100, 100);
      y[i] = gRandom->Uniform(-100, 100);
      z[i] = gRandom->Uniform(-100, 100);
   }

   Int_t diff1=0;

//for the distances brute-force:
   Double_t *dist = new Double_t[npoints];
   Int_t *index = new Int_t[npoints];

//Build the tree
   TKDTreeID *kdtree = new TKDTreeID(npoints, 3, bsize);
   kdtree->SetData(0, x);
   kdtree->SetData(1, y);
   kdtree->SetData(2, z);
   kdtree->Build();
   Int_t *index2 = new Int_t[nn];
   Double_t *dist2 = new Double_t[nn];
   Double_t point[3];
//Select a random point
   for (Int_t itime=0; itime<ntimes; itime++){
      Int_t ipoint = Int_t(gRandom->Uniform(0, npoints));

      for (Int_t i=0; i<npoints; i++){
         dist[i]=0;

         dist[i]+=(x[i]-x[ipoint])*(x[i]-x[ipoint]);
         dist[i]+=(y[i]-y[ipoint])*(y[i]-y[ipoint]);
         dist[i]+=(z[i]-z[ipoint])*(z[i]-z[ipoint]);
         dist[i]=TMath::Sqrt(dist[i]);

      }
      TMath::Sort(npoints, dist, index, kFALSE);

      point[0]=x[ipoint];
      point[1]=y[ipoint];
      point[2]=z[ipoint];

      kdtree->FindNearestNeighbors(point, nn, index2, dist2);


      for (Int_t inn=0; inn<nn; inn++){
         if (TMath::Abs(dist2[inn]-dist[index[inn]])>1E-8) {
            diff1++;
            // printf("dist1=%f, dist2=%f, in1=%lld, in2=%d\n", dist[index[inn]], dist2[inn], index[inn], index2[inn]);
         }


      }
   }

   printf("Nearest neighbors found for %d random points\n", ntimes);
   printf("%d neighbors are wrong compared to \"brute force\" method\n", diff1);
//   printf("Old: %d neighbors are wrong compared to brute-force method\n", diff2);

//    printf("\n");
//    for (Int_t i=0; i<nn; i++){
//       printf("ind[%d]=%d, dist[%d]=%f\n", i, index2[i], i, dist2[i]);
//    }




   delete [] x;
   delete [] y;
   delete [] z;
   delete [] index;
   delete [] dist;
   delete [] index2;
   delete [] dist2;
}

////////////////////////////////////////////////////////////////////////////////

void TestRange()
{
//Test TKDTree::FindInRange() function

   Int_t npoints = Int_t(gRandom->Uniform(0, 100000));
   Double_t range = gRandom->Uniform(20, 100);

   printf("%d points, range=%f\n", npoints, range);
   Int_t ntimes = 10;
   Double_t *x = new Double_t[npoints];
   Double_t *y = new Double_t[npoints];
   Double_t *z = new Double_t[npoints];
   for (Int_t i=0; i<npoints; i++){
      x[i] = gRandom->Uniform(-100, 100);
      y[i] = gRandom->Uniform(-100, 100);
      z[i] = gRandom->Uniform(-100, 100);
   }

   Int_t *results1 = new Int_t[npoints];
   std::vector<Int_t> results2;
   Int_t np1;

//Compute with the kd-tree
   Int_t bsize = 10;
   TKDTreeID *kdtree = new TKDTreeID(npoints, 3, bsize);
   kdtree->SetData(0, x);
   kdtree->SetData(1, y);
   kdtree->SetData(2, z);
   kdtree->Build();
   Double_t *dist = new Double_t[npoints];
   Int_t *index = new Int_t[npoints];
   Int_t ndiff = 0;
   for (Int_t itime=0; itime<ntimes; itime++){
      Double_t point[3];
      point[0]=gRandom->Uniform(-90, 90);
      point[1]=gRandom->Uniform(-90, 90);
      point[2]=gRandom->Uniform(-90, 90);

      //printf("point: (%f, %f, %f)\n\n", point[0], point[1], point[2]);
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
         dist[ipoint]=0;
         dist[ipoint]+=(x[ipoint]-point[0])*(x[ipoint]-point[0]);
         dist[ipoint]+=(y[ipoint]-point[1])*(y[ipoint]-point[1]);
         dist[ipoint]+=(z[ipoint]-point[2])*(z[ipoint]-point[2]);
         dist[ipoint]=TMath::Sqrt(dist[ipoint]);
         index[ipoint]=ipoint;
      }
      TMath::Sort(npoints, dist, index, kFALSE);
      np1=0;
      while (np1<npoints && dist[index[np1]]<=range){
         results1[np1]=index[np1];
         np1++;
      }
      results2.clear();
      kdtree->FindInRange(point, range, results2);

      if (TMath::Abs(np1 - Int_t(results2.size()))>0.1) {
         ndiff++;
         printf("different numbers of points found, %d %d\n", np1, Int_t(results2.size()));
         continue;
      }

      //have to sort the results, as they are in different order
      TMath::Sort(np1, results1, index, kFALSE);
      std::sort(results2.begin(), results2.end());

      for (Int_t i=0; i<np1; i++){
         if (TMath::Abs(results1[index[i]]-results2[i])>1E-8) ndiff++;
      }
   }
   printf("%d  differences found between \"brute force\" method and kd-tree\n", ndiff);

   delete [] x;
   delete [] y;
   delete [] z;
   delete [] index;
   delete [] dist;
   delete [] results1;
   delete kdtree;
}




////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // Parse command line arguments
  for (Int_t i=1 ;  i<argc ; i++) {
     std::string arg = argv[i] ;
     if (arg == "-g") {
      showGraphics = true;
     }
     // if (arg == "-v") {
     //  showGraphics = true;
     //  verbose = true;
     // }
     if (arg == "-h") {
        std::cerr << "Usage: " << argv[0] << " [-g] [-v]\n";
        std::cerr << "  where:\n";
        std::cerr << "     -g : graphics mode\n";
        //std::cerr << "     -v : verbose  mode";
        std::cerr << std::endl;
        return -1;
     }
   }

   TApplication* theApp = 0;
   if ( showGraphics )
      theApp = new TApplication("App",&argc,argv);

   kDTreeTest();

   if ( showGraphics )
   {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return 0;
}
