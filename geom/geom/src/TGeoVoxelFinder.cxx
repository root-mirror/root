// @(#)root/geom:$Id$
// Author: Andrei Gheata   04/02/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoVoxelFinder
\ingroup Geometry_classes

Finder class handling voxels.

Full description with examples and pictures

\image html geom_t_finder.png
\image html geom_t_voxelfind.png
\image html geom_t_voxtree.png
*/

#include "TGeoVoxelFinder.h"

#include "TBuffer.h"
#include "TMath.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoStateInfo.h"

ClassImp(TGeoVoxelFinder);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVoxelFinder::TGeoVoxelFinder()
{
   fVolume  = 0;
   fNboxes  = 0;
   fIbx     = 0;
   fIby     = 0;
   fIbz     = 0;
   fNox     = 0;
   fNoy     = 0;
   fNoz     = 0;
   fNex     = 0;
   fNey     = 0;
   fNez     = 0;
   fNx      = 0;
   fNy      = 0;
   fNz      = 0;
   fBoxes   = 0;
   fXb      = 0;
   fYb      = 0;
   fZb      = 0;
   fOBx     = 0;
   fOBy     = 0;
   fOBz     = 0;
   fOEx     = 0;
   fOEy     = 0;
   fOEz     = 0;
   fIndcX   = 0;
   fIndcY   = 0;
   fIndcZ   = 0;
   fExtraX  = 0;
   fExtraY  = 0;
   fExtraZ  = 0;
   fNsliceX = 0;
   fNsliceY = 0;
   fNsliceZ = 0;
   memset(fPriority, 0, 3*sizeof(Int_t));
   SetInvalid(kFALSE);
}
////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVoxelFinder::TGeoVoxelFinder(TGeoVolume *vol)
{
   if (!vol) {
      Fatal("TGeoVoxelFinder", "Null pointer for volume");
      return; // To make code checkers happy
   }
   fVolume  = vol;
   fVolume->SetCylVoxels(kFALSE);
   fNboxes   = 0;
   fIbx      = 0;
   fIby      = 0;
   fIbz      = 0;
   fNox      = 0;
   fNoy      = 0;
   fNoz      = 0;
   fNex     = 0;
   fNey     = 0;
   fNez     = 0;
   fNx       = 0;
   fNy       = 0;
   fNz       = 0;
   fBoxes    = 0;
   fXb       = 0;
   fYb       = 0;
   fZb       = 0;
   fOBx      = 0;
   fOBy      = 0;
   fOBz      = 0;
   fOEx     = 0;
   fOEy     = 0;
   fOEz     = 0;
   fIndcX   = 0;
   fIndcY   = 0;
   fIndcZ   = 0;
   fExtraX  = 0;
   fExtraY  = 0;
   fExtraZ  = 0;
   fNsliceX = 0;
   fNsliceY = 0;
   fNsliceZ = 0;
   memset(fPriority, 0, 3*sizeof(Int_t));
   SetNeedRebuild();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoVoxelFinder::~TGeoVoxelFinder()
{
   if (fOBx) delete [] fOBx;
   if (fOBy) delete [] fOBy;
   if (fOBz) delete [] fOBz;
   if (fOEx) delete [] fOEx;
   if (fOEy) delete [] fOEy;
   if (fOEz) delete [] fOEz;
//   printf("OBx OBy OBz...\n");
   if (fBoxes) delete [] fBoxes;
//   printf("boxes...\n");
   if (fXb) delete [] fXb;
   if (fYb) delete [] fYb;
   if (fZb) delete [] fZb;
//   printf("Xb Yb Zb...\n");
   if (fNsliceX) delete [] fNsliceX;
   if (fNsliceY) delete [] fNsliceY;
   if (fNsliceZ) delete [] fNsliceZ;
   if (fIndcX) delete [] fIndcX;
   if (fIndcY) delete [] fIndcY;
   if (fIndcZ) delete [] fIndcZ;
   if (fExtraX) delete [] fExtraX;
   if (fExtraY) delete [] fExtraY;
   if (fExtraZ) delete [] fExtraZ;
//   printf("IndX IndY IndZ...\n");
}

////////////////////////////////////////////////////////////////////////////////

Int_t TGeoVoxelFinder::GetNcandidates(TGeoStateInfo &td) const
{
   return td.fVoxNcandidates;
}

////////////////////////////////////////////////////////////////////////////////

Int_t* TGeoVoxelFinder::GetCheckList(Int_t &nelem, TGeoStateInfo &td) const
{
   nelem = td.fVoxNcandidates;
   return td.fVoxCheckList;
}

////////////////////////////////////////////////////////////////////////////////
/// build the array of bounding boxes of the nodes inside

void TGeoVoxelFinder::BuildVoxelLimits()
{
   Int_t nd = fVolume->GetNdaughters();
   if (!nd) return;
   Int_t id;
   TGeoNode *node;
   if (fBoxes) delete [] fBoxes;
   fNboxes = 6*nd;
   fBoxes = new Double_t[fNboxes];
   Double_t vert[24] = {0};
   Double_t pt[3] = {0};
   Double_t xyz[6] = {0};
//   printf("boundaries for %s :\n", GetName());
   TGeoBBox *box = 0;
   for (id=0; id<nd; id++) {
      node = fVolume->GetNode(id);
//      if (!strcmp(node->ClassName(), "TGeoNodeOffset") continue;
      box = (TGeoBBox*)node->GetVolume()->GetShape();
      box->SetBoxPoints(&vert[0]);
      for (Int_t point=0; point<8; point++) {
         DaughterToMother(id, &vert[3*point], &pt[0]);
         if (!point) {
            xyz[0] = xyz[1] = pt[0];
            xyz[2] = xyz[3] = pt[1];
            xyz[4] = xyz[5] = pt[2];
            continue;
         }
         for (Int_t j=0; j<3; j++) {
            if (pt[j] < xyz[2*j]) xyz[2*j]=pt[j];
            if (pt[j] > xyz[2*j+1]) xyz[2*j+1]=pt[j];
         }
      }
      fBoxes[6*id]   = 0.5*(xyz[1]-xyz[0]); // dX
      fBoxes[6*id+1] = 0.5*(xyz[3]-xyz[2]); // dY
      fBoxes[6*id+2] = 0.5*(xyz[5]-xyz[4]); // dZ
      fBoxes[6*id+3] = 0.5*(xyz[0]+xyz[1]); // Ox
      fBoxes[6*id+4] = 0.5*(xyz[2]+xyz[3]); // Oy
      fBoxes[6*id+5] = 0.5*(xyz[4]+xyz[5]); // Oz
   }
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point from the local reference system of node id to reference
/// system of mother volume

void TGeoVoxelFinder::DaughterToMother(Int_t id, const Double_t *local, Double_t *master) const
{
   TGeoMatrix *mat = fVolume->GetNode(id)->GetMatrix();
   if (!mat) memcpy(master,local,3*sizeof(Double_t));
   else      mat->LocalToMaster(local, master);
}
////////////////////////////////////////////////////////////////////////////////
/// Computes squared distance from POINT to the voxel(s) containing node INODE. Returns 0
/// if POINT inside voxel(s).

Bool_t TGeoVoxelFinder::IsSafeVoxel(const Double_t *point, Int_t inode, Double_t minsafe) const
{
   if (NeedRebuild()) {
      TGeoVoxelFinder *vox = (TGeoVoxelFinder*)this;
      vox->Voxelize();
      fVolume->FindOverlaps();
   }
   Double_t dxyz, minsafe2=minsafe*minsafe;
   Int_t ist = 6*inode;
   Int_t i;
   Double_t rsq = 0;
   for (i=0; i<3; i++) {
      dxyz = TMath::Abs(point[i]-fBoxes[ist+i+3])-fBoxes[ist+i];
      if (dxyz>-1E-6) rsq+=dxyz*dxyz;
      if (rsq > minsafe2*(1.+TGeoShape::Tolerance())) return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute voxelization efficiency.

Double_t TGeoVoxelFinder::Efficiency()
{
   printf("Voxelization efficiency for %s\n", fVolume->GetName());
   if (NeedRebuild()) {
      Voxelize();
      fVolume->FindOverlaps();
   }
   Double_t nd = Double_t(fVolume->GetNdaughters());
   Double_t eff = 0;
   Double_t effslice = 0;
   Int_t id;
   if (fPriority[0]) {
      for (id=0; id<fIbx-1; id++) {  // loop on boundaries
         effslice += fNsliceX[id];
      }
      if (!TGeoShape::IsSameWithinTolerance(effslice,0)) effslice = nd/effslice;
      else printf("Woops : slice X\n");
   }
   printf("X efficiency : %g\n", effslice);
   eff += effslice;
   effslice = 0;
   if (fPriority[1]) {
      for (id=0; id<fIby-1; id++) {  // loop on boundaries
         effslice += fNsliceY[id];
      }
      if (!TGeoShape::IsSameWithinTolerance(effslice,0)) effslice = nd/effslice;
      else printf("Woops : slice X\n");
   }
   printf("Y efficiency : %g\n", effslice);
   eff += effslice;
   effslice = 0;
   if (fPriority[2]) {
      for (id=0; id<fIbz-1; id++) {  // loop on boundaries
         effslice += fNsliceZ[id];
      }
      if (!TGeoShape::IsSameWithinTolerance(effslice,0)) effslice = nd/effslice;
      else printf("Woops : slice X\n");
   }
   printf("Z efficiency : %g\n", effslice);
   eff += effslice;
   eff /= 3.;
   printf("Total efficiency : %g\n", eff);
   return eff;
}
////////////////////////////////////////////////////////////////////////////////
/// create the list of nodes for which the bboxes overlap with inode's bbox

void TGeoVoxelFinder::FindOverlaps(Int_t inode) const
{
   if (!fBoxes) return;
   Double_t xmin, xmax, ymin, ymax, zmin, zmax;
   Double_t xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
   Double_t ddx1, ddx2;
   Int_t nd = fVolume->GetNdaughters();
   Int_t *ovlps = 0;
   Int_t *otmp = new Int_t[nd-1];
   Int_t novlp = 0;
   TGeoNode *node = fVolume->GetNode(inode);
   xmin = fBoxes[6*inode+3] - fBoxes[6*inode];
   xmax = fBoxes[6*inode+3] + fBoxes[6*inode];
   ymin = fBoxes[6*inode+4] - fBoxes[6*inode+1];
   ymax = fBoxes[6*inode+4] + fBoxes[6*inode+1];
   zmin = fBoxes[6*inode+5] - fBoxes[6*inode+2];
   zmax = fBoxes[6*inode+5] + fBoxes[6*inode+2];
   // loop on brothers
   for (Int_t ib=0; ib<nd; ib++) {
      if (ib == inode) continue; // everyone overlaps with itself
      xmin1 = fBoxes[6*ib+3] - fBoxes[6*ib];
      xmax1 = fBoxes[6*ib+3] + fBoxes[6*ib];
      ymin1 = fBoxes[6*ib+4] - fBoxes[6*ib+1];
      ymax1 = fBoxes[6*ib+4] + fBoxes[6*ib+1];
      zmin1 = fBoxes[6*ib+5] - fBoxes[6*ib+2];
      zmax1 = fBoxes[6*ib+5] + fBoxes[6*ib+2];

      ddx1 = xmax-xmin1;
      ddx2 = xmax1-xmin;
      if (ddx1*ddx2 <= 0.) continue;
      ddx1 = ymax-ymin1;
      ddx2 = ymax1-ymin;
      if (ddx1*ddx2 <= 0.) continue;
      ddx1 = zmax-zmin1;
      ddx2 = zmax1-zmin;
      if (ddx1*ddx2 <= 0.) continue;
      otmp[novlp++] = ib;
   }
   if (!novlp) {
      delete [] otmp;
      node->SetOverlaps(ovlps, 0);
      return;
   }
   ovlps = new Int_t[novlp];
   memcpy(ovlps, otmp, novlp*sizeof(Int_t));
   delete [] otmp;
   node->SetOverlaps(ovlps, novlp);
}

////////////////////////////////////////////////////////////////////////////////
/// Get indices for current slices on x, y, z

Bool_t TGeoVoxelFinder::GetIndices(const Double_t *point, TGeoStateInfo &td)
{
   td.fVoxSlices[0] = -2; // -2 means 'all daughters in slice'
   td.fVoxSlices[1] = -2;
   td.fVoxSlices[2] = -2;
   Bool_t flag=kTRUE;
   if (fPriority[0]) {
      td.fVoxSlices[0] = TMath::BinarySearch(fIbx, fXb, point[0]);
      if ((td.fVoxSlices[0]<0) || (td.fVoxSlices[0]>=fIbx-1)) {
         // outside slices
         flag=kFALSE;
      } else {
         if (fPriority[0]==2) {
            // nothing in current slice
            if (!fNsliceX[td.fVoxSlices[0]]) flag = kFALSE;
         }
      }
   }
   if (fPriority[1]) {
      td.fVoxSlices[1] = TMath::BinarySearch(fIby, fYb, point[1]);
      if ((td.fVoxSlices[1]<0) || (td.fVoxSlices[1]>=fIby-1)) {
         // outside slices
         flag=kFALSE;
      } else {
         if (fPriority[1]==2) {
            // nothing in current slice
            if (!fNsliceY[td.fVoxSlices[1]]) flag = kFALSE;
         }
      }
   }
   if (fPriority[2]) {
      td.fVoxSlices[2] = TMath::BinarySearch(fIbz, fZb, point[2]);
      if ((td.fVoxSlices[2]<0) || (td.fVoxSlices[2]>=fIbz-1)) return kFALSE;
      if (fPriority[2]==2) {
         // nothing in current slice
         if (!fNsliceZ[td.fVoxSlices[2]]) return kFALSE;
      }
   }
   return flag;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of extra candidates in a given X slice compared to
/// another (left or right)

Int_t *TGeoVoxelFinder::GetExtraX(Int_t islice, Bool_t left, Int_t &nextra) const
{
   Int_t *list = 0;
   nextra = 0;
   if (fPriority[0]!=2) return list;
   if (left) {
      nextra = fExtraX[fOEx[islice]];
      list = &fExtraX[fOEx[islice]+2];
   } else {
      nextra = fExtraX[fOEx[islice]+1];
      list = &fExtraX[fOEx[islice]+2+fExtraX[fOEx[islice]]];
   }
   return list;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of extra candidates in a given Y slice compared to
/// another (left or right)

Int_t *TGeoVoxelFinder::GetExtraY(Int_t islice, Bool_t left, Int_t &nextra) const
{
   Int_t *list = 0;
   nextra = 0;
   if (fPriority[1]!=2) return list;
   if (left) {
      nextra = fExtraY[fOEy[islice]];
      list = &fExtraY[fOEy[islice]+2];
   } else {
      nextra = fExtraY[fOEy[islice]+1];
      list = &fExtraY[fOEy[islice]+2+fExtraY[fOEy[islice]]];
   }
   return list;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of extra candidates in a given Z slice compared to
/// another (left or right)

Int_t *TGeoVoxelFinder::GetExtraZ(Int_t islice, Bool_t left, Int_t &nextra) const
{
   Int_t *list = 0;
   nextra = 0;
   if (fPriority[2]!=2) return list;
   if (left) {
      nextra = fExtraZ[fOEz[islice]];
      list = &fExtraZ[fOEz[islice]+2];
   } else {
      nextra = fExtraZ[fOEz[islice]+1];
      list = &fExtraZ[fOEz[islice]+2+fExtraZ[fOEz[islice]]];
   }
   return list;
}

////////////////////////////////////////////////////////////////////////////////
/// Get extra candidates that are not contained in current check list

Int_t *TGeoVoxelFinder::GetValidExtra(Int_t *list, Int_t &ncheck, TGeoStateInfo &td)
{
   td.fVoxNcandidates = 0;
   Int_t icand;
   UInt_t bitnumber, loc;
   UChar_t bit, byte;
   for (icand=0; icand<ncheck; icand++) {
      bitnumber = (UInt_t)list[icand];
      loc = bitnumber>>3;
      bit = bitnumber%8;
      byte = (~td.fVoxBits1[loc]) & (1<<bit);
      if (byte) td.fVoxCheckList[td.fVoxNcandidates++]=list[icand];
   }
   ncheck = td.fVoxNcandidates;
   return td.fVoxCheckList;
}

////////////////////////////////////////////////////////////////////////////////
/// Get extra candidates that are contained in array1 but not in current check list

Int_t *TGeoVoxelFinder::GetValidExtra(Int_t /*n1*/, UChar_t *array1, Int_t *list, Int_t &ncheck, TGeoStateInfo &td)
{
   td.fVoxNcandidates = 0;
   Int_t icand;
   UInt_t bitnumber, loc;
   UChar_t bit, byte;
   for (icand=0; icand<ncheck; icand++) {
      bitnumber = (UInt_t)list[icand];
      loc = bitnumber>>3;
      bit = bitnumber%8;
      byte = (~td.fVoxBits1[loc]) & array1[loc] & (1<<bit);
      if (byte) td.fVoxCheckList[td.fVoxNcandidates++]=list[icand];
   }
   ncheck = td.fVoxNcandidates;
   return td.fVoxCheckList;
}

////////////////////////////////////////////////////////////////////////////////
/// Get extra candidates that are contained in array1 but not in current check list

Int_t *TGeoVoxelFinder::GetValidExtra(Int_t /*n1*/, UChar_t *array1, Int_t /*n2*/, UChar_t *array2, Int_t *list, Int_t &ncheck, TGeoStateInfo &td)
{
   td.fVoxNcandidates = 0;
   Int_t icand;
   UInt_t bitnumber, loc;
   UChar_t bit, byte;
   for (icand=0; icand<ncheck; icand++) {
      bitnumber = (UInt_t)list[icand];
      loc = bitnumber>>3;
      bit = bitnumber%8;
      byte = (~td.fVoxBits1[loc]) & array1[loc] & array2[loc] & (1<<bit);
      if (byte) td.fVoxCheckList[td.fVoxNcandidates++]=list[icand];
   }
   ncheck = td.fVoxNcandidates;
   return td.fVoxCheckList;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns list of new candidates in next voxel. If NULL, nowhere to
/// go next.

Int_t *TGeoVoxelFinder::GetNextCandidates(const Double_t *point, Int_t &ncheck, TGeoStateInfo &td)
{
   if (NeedRebuild()) {
      Voxelize();
      fVolume->FindOverlaps();
   }
   ncheck = 0;
   if (td.fVoxLimits[0]<0) return 0;
   if (td.fVoxLimits[1]<0) return 0;
   if (td.fVoxLimits[2]<0) return 0;
   Int_t dind[3]; // new slices
   //---> start from old slices
   memcpy(&dind[0], &td.fVoxSlices[0], 3*sizeof(Int_t));
   Double_t dmin[3]; // distances to get to next X,Y, Z slices.
   dmin[0] = dmin[1] = dmin[2] = TGeoShape::Big();
   //---> max. possible step to be considered
   Double_t maxstep = TMath::Min(gGeoManager->GetStep(), td.fVoxLimits[TMath::LocMin(3, td.fVoxLimits)]);
//   printf("1- maxstep=%g\n", maxstep);
   Bool_t isXlimit=kFALSE, isYlimit=kFALSE, isZlimit=kFALSE;
   Bool_t isForcedX=kFALSE, isForcedY=kFALSE, isForcedZ=kFALSE;
   Double_t dforced[3];
   dforced[0] = dforced[1] = dforced[2] = TGeoShape::Big();
   Int_t iforced = 0;
   //
   //---> work on X
   if (fPriority[0] && td.fVoxInc[0]) {
      //---> increment/decrement slice
      dind[0] += td.fVoxInc[0];
      if (td.fVoxInc[0]==1) {
         if (dind[0]<0 || dind[0]>fIbx-1) return 0; // outside range
         dmin[0] = (fXb[dind[0]]-point[0])*td.fVoxInvdir[0];
      } else {
         if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>fIbx-1) return 0; // outside range
         dmin[0] = (fXb[td.fVoxSlices[0]]-point[0])*td.fVoxInvdir[0];
      }
      isXlimit = (dmin[0]>maxstep)?kTRUE:kFALSE;
//      printf("---> X : priority=%i, slice=%i/%i inc=%i\n",
//             fPriority[0], td.fVoxSlices[0], fIbx-2, td.fVoxInc[0]);
//      printf("2- step to next X (%i) = %g\n", (Int_t)isXlimit, dmin[0]);
      //---> check if propagation to next slice on this axis is forced
      if ((td.fVoxSlices[0]==-1) || (td.fVoxSlices[0]==fIbx-1)) {
         isForcedX = kTRUE;
         dforced[0] = dmin[0];
         iforced++;
//         printf("   FORCED 1\n");
         if (isXlimit) return 0;
      } else {
         if (fPriority[0]==2) {
            // if no candidates in current slice, force next slice
            if (fNsliceX[td.fVoxSlices[0]]==0) {
               isForcedX = kTRUE;
               dforced[0] = dmin[0];
               iforced++;
//               printf("   FORCED 2\n");
               if (isXlimit) return 0;
            }
         }
      }
   } else {
      // no slices on this axis -> bounding box limit
//      printf("   No slice on X\n");
      dmin[0] = td.fVoxLimits[0];
      isXlimit = kTRUE;
   }
   //---> work on Y
   if (fPriority[1] && td.fVoxInc[1]) {
      //---> increment/decrement slice
      dind[1] += td.fVoxInc[1];
      if (td.fVoxInc[1]==1) {
         if (dind[1]<0 || dind[1]>fIby-1) return 0; // outside range
         dmin[1] = (fYb[dind[1]]-point[1])*td.fVoxInvdir[1];
      } else {
         if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>fIby-1) return 0; // outside range
         dmin[1] = (fYb[td.fVoxSlices[1]]-point[1])*td.fVoxInvdir[1];
      }
      isYlimit = (dmin[1]>maxstep)?kTRUE:kFALSE;
//      printf("---> Y : priority=%i, slice=%i/%i inc=%i\n",
//             fPriority[1], td.fVoxSlices[1], fIby-2, td.fVoxInc[1]);
//      printf("3- step to next Y (%i) = %g\n", (Int_t)isYlimit, dmin[1]);

      //---> check if propagation to next slice on this axis is forced
      if ((td.fVoxSlices[1]==-1) || (td.fVoxSlices[1]==fIby-1)) {
         isForcedY = kTRUE;
         dforced[1] = dmin[1];
         iforced++;
//         printf("   FORCED 1\n");
         if (isYlimit) return 0;
      } else {
         if (fPriority[1]==2) {
            // if no candidates in current slice, force next slice
            if (fNsliceY[td.fVoxSlices[1]]==0) {
               isForcedY = kTRUE;
               dforced[1] = dmin[1];
               iforced++;
//               printf("   FORCED 2\n");
               if (isYlimit) return 0;
            }
         }
      }
   } else {
      // no slices on this axis -> bounding box limit
//      printf("   No slice on Y\n");
      dmin[1] = td.fVoxLimits[1];
      isYlimit = kTRUE;
   }
   //---> work on Z
   if (fPriority[2] && td.fVoxInc[2]) {
      //---> increment/decrement slice
      dind[2] += td.fVoxInc[2];
      if (td.fVoxInc[2]==1) {
         if (dind[2]<0 || dind[2]>fIbz-1) return 0; // outside range
         dmin[2] = (fZb[dind[2]]-point[2])*td.fVoxInvdir[2];
      } else {
         if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>fIbz-1) return 0; // outside range
         dmin[2] = (fZb[td.fVoxSlices[2]]-point[2])*td.fVoxInvdir[2];
      }
      isZlimit = (dmin[2]>maxstep)?kTRUE:kFALSE;
//      printf("---> Z : priority=%i, slice=%i/%i inc=%i\n",
//             fPriority[2], td.fVoxSlices[2], fIbz-2, td.fVoxInc[2]);
//      printf("4- step to next Z (%i) = %g\n", (Int_t)isZlimit, dmin[2]);

      //---> check if propagation to next slice on this axis is forced
      if ((td.fVoxSlices[2]==-1) || (td.fVoxSlices[2]==fIbz-1)) {
         isForcedZ = kTRUE;
         dforced[2] = dmin[2];
         iforced++;
//         printf("   FORCED 1\n");
         if (isZlimit) return 0;
      } else {
         if (fPriority[2]==2) {
            // if no candidates in current slice, force next slice
            if (fNsliceZ[td.fVoxSlices[2]]==0) {
               isForcedZ = kTRUE;
               dforced[2] = dmin[2];
               iforced++;
//               printf("   FORCED 2\n");
               if (isZlimit) return 0;
            }
         }
      }
   } else {
      // no slices on this axis -> bounding box limit
//      printf("   No slice on Z\n");
      dmin[2] = td.fVoxLimits[2];
      isZlimit = kTRUE;
   }
   //---> We are done with checking. See which is the closest slice.
   // First check if some slice is forced

   Double_t dslice = 0;
   Int_t islice = 0;
   if (iforced) {
   // some slice is forced
      if (isForcedX) {
      // X forced
         dslice = dforced[0];
         islice = 0;
         if (isForcedY) {
         // X+Y forced
            if (dforced[1]>dslice) {
               dslice = dforced[1];
               islice = 1;
            }
            if (isForcedZ) {
            // X+Y+Z forced
               if (dforced[2]>dslice) {
                  dslice = dforced[2];
                  islice = 2;
               }
            }
         } else {
         // X forced
            if (isForcedZ) {
            // X+Z forced
               if (dforced[2]>dslice) {
                  dslice = dforced[2];
                  islice = 2;
               }
            }
         }
      } else {
         if (isForcedY) {
         // Y forced
            dslice = dforced[1];
            islice = 1;
            if (isForcedZ) {
            // Y+Z forced
               if (dforced[2]>dslice) {
                  dslice = dforced[2];
                  islice = 2;
               }
            }
         } else {
         // Z forced
            dslice = dforced[2];
            islice = 2;
         }
      }
   } else {
   // Nothing forced -> get minimum distance
      islice = TMath::LocMin(3, dmin);
      dslice = dmin[islice];
      if (dslice>=maxstep) {
//         printf("DSLICE > MAXSTEP -> EXIT\n");
         return 0;
      }
   }
//   printf("5- islicenext=%i  DSLICE=%g\n", islice, dslice);
   Double_t xptnew;
   Int_t *new_list; // list of new candidates
   UChar_t *slice1 = 0;
   UChar_t *slice2 = 0;
   Int_t ndd[2] = {0,0};
   Int_t islices = 0;
   Bool_t left;
   switch (islice) {
      case 0:
         if (isXlimit) return 0;
         // increment/decrement X slice
         td.fVoxSlices[0]=dind[0];
         if (iforced) {
            // we have to recompute Y and Z slices
            if (dslice>td.fVoxLimits[1]) return 0;
            if (dslice>td.fVoxLimits[2]) return 0;
            if ((dslice>dmin[1]) && td.fVoxInc[1]) {
               xptnew = point[1]+dslice/td.fVoxInvdir[1];
               //               printf("   recomputing Y slice, pos=%g\n", xptnew);
               while (1) {
                  td.fVoxSlices[1] += td.fVoxInc[1];
                  if (td.fVoxInc[1]==1) {
                     if (td.fVoxSlices[1]<-1 || td.fVoxSlices[1]>fIby-2) break; // outside range
                     if (fYb[td.fVoxSlices[1]+1]>=xptnew) break;
                  } else {
                     if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>fIby-1) break; // outside range
                     if (fYb[td.fVoxSlices[1]]<= xptnew) break;
                  }
               }
               //               printf("   %i/%i\n", td.fVoxSlices[1], fIby-2);
            }
            if ((dslice>dmin[2]) && td.fVoxInc[2]) {
               xptnew = point[2]+dslice/td.fVoxInvdir[2];
               //               printf("   recomputing Z slice, pos=%g\n", xptnew);
               while (1) {
                  td.fVoxSlices[2] += td.fVoxInc[2];
                  if (td.fVoxInc[2]==1) {
                     if (td.fVoxSlices[2]<-1 || td.fVoxSlices[2]>fIbz-2) break; // outside range
                     if (fZb[td.fVoxSlices[2]+1]>=xptnew) break;
                  } else {
                     if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>fIbz-1) break; // outside range
                     if (fZb[td.fVoxSlices[2]]<= xptnew) break;
                  }
               }
               //               printf("   %i/%i\n", td.fVoxSlices[2], fIbz-2);
            }
         }
         // new indices are set -> Get new candidates
         if (fPriority[0]==1) {
            // we are entering the unique slice on this axis
            //---> intersect and store Y and Z
            if (fPriority[1]==2) {
               if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>=fIby-1) return td.fVoxCheckList; // outside range
               ndd[0] = fNsliceY[td.fVoxSlices[1]];
               if (!ndd[0]) return td.fVoxCheckList;
               slice1 = &fIndcY[fOBy[td.fVoxSlices[1]]];
               islices++;
            }
            if (fPriority[2]==2) {
               if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>=fIbz-1) return td.fVoxCheckList; // outside range
               ndd[1] = fNsliceZ[td.fVoxSlices[2]];
               if (!ndd[1]) return td.fVoxCheckList;
               islices++;
               if (slice1) {
                  slice2 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
               } else {
                  slice1 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
                  ndd[0] = ndd[1];
               }
            }
            if (islices<=1) {
               IntersectAndStore(ndd[0], slice1, td);
            } else {
               IntersectAndStore(ndd[0], slice1, ndd[1], slice2, td);
            }
            ncheck = td.fVoxNcandidates;
            return td.fVoxCheckList;
         }
         // We got into a new slice -> Get only new candidates
         left = (td.fVoxInc[0]>0)?kTRUE:kFALSE;
         new_list = GetExtraX(td.fVoxSlices[0], left, ncheck);
//         printf("   New list on X : %i new candidates\n", ncheck);
         if (!ncheck) return td.fVoxCheckList;
         if (fPriority[1]==2) {
            if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>=fIby-1) {
               ncheck = 0;
               return td.fVoxCheckList; // outside range
            }
            ndd[0] = fNsliceY[td.fVoxSlices[1]];
            if (!ndd[0]) {
               ncheck = 0;
               return td.fVoxCheckList;
            }
            slice1 = &fIndcY[fOBy[td.fVoxSlices[1]]];
            islices++;
         }
         if (fPriority[2]==2) {
            if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>=fIbz-1) {
               ncheck = 0;
               return td.fVoxCheckList; // outside range
            }
            ndd[1] = fNsliceZ[td.fVoxSlices[2]];
            if (!ndd[1]) {
               ncheck = 0;
               return td.fVoxCheckList;
            }
            islices++;
            if (slice1) {
               slice2 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
            } else {
               slice1 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
               ndd[0] = ndd[1];
            }
         }
         if (!islices) return GetValidExtra(new_list, ncheck, td);
         if (islices==1) {
            return GetValidExtra(ndd[0], slice1, new_list, ncheck,td);
         } else {
            return GetValidExtra(ndd[0], slice1, ndd[1], slice2, new_list, ncheck, td);
         }
      case 1:
         if (isYlimit) return 0;
         // increment/decrement Y slice
         td.fVoxSlices[1]=dind[1];
         if (iforced) {
         // we have to recompute X and Z slices
            if (dslice>td.fVoxLimits[0]) return 0;
            if (dslice>td.fVoxLimits[2]) return 0;
            if ((dslice>dmin[0]) && td.fVoxInc[0]) {
               xptnew = point[0]+dslice/td.fVoxInvdir[0];
//               printf("   recomputing X slice, pos=%g\n", xptnew);
               while (1) {
                  td.fVoxSlices[0] += td.fVoxInc[0];
                  if (td.fVoxInc[0]==1) {
                     if (td.fVoxSlices[0]<-1 || td.fVoxSlices[0]>fIbx-2) break; // outside range
                     if (fXb[td.fVoxSlices[0]+1]>=xptnew) break;
                  } else {
                     if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>fIbx-1) break; // outside range
                     if (fXb[td.fVoxSlices[0]]<= xptnew) break;
                  }
               }
//               printf("   %i/%i\n", td.fVoxSlices[0], fIbx-2);
            }
            if ((dslice>dmin[2]) && td.fVoxInc[2]) {
               xptnew = point[2]+dslice/td.fVoxInvdir[2];
//               printf("   recomputing Z slice, pos=%g\n", xptnew);
               while (1) {
                  td.fVoxSlices[2] += td.fVoxInc[2];
                  if (td.fVoxInc[2]==1) {
                     if (td.fVoxSlices[2]<-1 || td.fVoxSlices[2]>fIbz-2) break; // outside range
                     if (fZb[td.fVoxSlices[2]+1]>=xptnew) break;
                  } else {
                     if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>fIbz-1) break; // outside range
                     if (fZb[td.fVoxSlices[2]]<= xptnew) break;
                  }
               }
//               printf("   %i/%i\n", td.fVoxSlices[2], fIbz-2);
            }
         }
         // new indices are set -> Get new candidates
         if (fPriority[1]==1) {
         // we are entering the unique slice on this axis
         //---> intersect and store X and Z
            if (fPriority[0]==2) {
               if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>=fIbx-1) return td.fVoxCheckList; // outside range
               ndd[0] = fNsliceX[td.fVoxSlices[0]];
               if (!ndd[0]) return td.fVoxCheckList;
               slice1 = &fIndcX[fOBx[td.fVoxSlices[0]]];
               islices++;
            }
            if (fPriority[2]==2) {
               if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>=fIbz-1) return td.fVoxCheckList; // outside range
               ndd[1] = fNsliceZ[td.fVoxSlices[2]];
               if (!ndd[1]) return td.fVoxCheckList;
               islices++;
               if (slice1) {
                  slice2 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
               } else {
                  slice1 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
                  ndd[0] = ndd[1];
               }
            }
            if (islices<=1) {
               IntersectAndStore(ndd[0], slice1, td);
            } else {
               IntersectAndStore(ndd[0], slice1, ndd[1], slice2, td);
            }
            ncheck = td.fVoxNcandidates;
            return td.fVoxCheckList;
         }
         // We got into a new slice -> Get only new candidates
         left = (td.fVoxInc[1]>0)?kTRUE:kFALSE;
         new_list = GetExtraY(td.fVoxSlices[1], left, ncheck);
//         printf("   New list on Y : %i new candidates\n", ncheck);
         if (!ncheck) return td.fVoxCheckList;
         if (fPriority[0]==2) {
            if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>=fIbx-1) {
               ncheck = 0;
               return td.fVoxCheckList; // outside range
            }
            ndd[0] = fNsliceX[td.fVoxSlices[0]];
            if (!ndd[0]) {
               ncheck = 0;
               return td.fVoxCheckList;
            }
            slice1 = &fIndcX[fOBx[td.fVoxSlices[0]]];
            islices++;
         }
         if (fPriority[2]==2) {
            if (td.fVoxSlices[2]<0 || td.fVoxSlices[2]>=fIbz-1) {
               ncheck = 0;
               return td.fVoxCheckList; // outside range
            }
            ndd[1] = fNsliceZ[td.fVoxSlices[2]];
            if (!ndd[1]) {
               ncheck = 0;
               return td.fVoxCheckList;
            }
            islices++;
            if (slice1) {
               slice2 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
            } else {
               slice1 = &fIndcZ[fOBz[td.fVoxSlices[2]]];
               ndd[0] = ndd[1];
            }
         }
         if (!islices) return GetValidExtra(new_list, ncheck, td);
         if (islices==1) {
            return GetValidExtra(ndd[0], slice1, new_list, ncheck, td);
         } else {
            return GetValidExtra(ndd[0], slice1, ndd[1], slice2, new_list, ncheck, td);
         }
      case 2:
         if (isZlimit) return 0;
         // increment/decrement Z slice
         td.fVoxSlices[2]=dind[2];
         if (iforced) {
            // we have to recompute Y and X slices
            if (dslice>td.fVoxLimits[1]) return 0;
            if (dslice>td.fVoxLimits[0]) return 0;
            if ((dslice>dmin[1]) && td.fVoxInc[1]) {
               xptnew = point[1]+dslice/td.fVoxInvdir[1];
               //               printf("   recomputing Y slice, pos=%g\n", xptnew);
               while (1) {
                  td.fVoxSlices[1] += td.fVoxInc[1];
                  if (td.fVoxInc[1]==1) {
                     if (td.fVoxSlices[1]<-1 || td.fVoxSlices[1]>fIby-2) break; // outside range
                     if (fYb[td.fVoxSlices[1]+1]>=xptnew) break;
                  } else {
                     if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>fIby-1) break; // outside range
                     if (fYb[td.fVoxSlices[1]]<= xptnew) break;
                  }
               }
               //               printf("   %i/%i\n", td.fVoxSlices[1], fIby-2);
            }
            if ((dslice>dmin[0]) && td.fVoxInc[0]) {
               xptnew = point[0]+dslice/td.fVoxInvdir[0];
               //               printf("   recomputing X slice, pos=%g\n", xptnew);
               while (1) {
                  td.fVoxSlices[0] += td.fVoxInc[0];
                  if (td.fVoxInc[0]==1) {
                     if (td.fVoxSlices[0]<-1 || td.fVoxSlices[0]>fIbx-2) break; // outside range
                     if (fXb[td.fVoxSlices[0]+1]>=xptnew) break;
                  } else {
                     if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>fIbx-1) break; // outside range
                     if (fXb[td.fVoxSlices[0]]<= xptnew) break;
                  }
               }
//               printf("   %i/%i\n", td.fVoxSlices[0], fIbx-2);
            }
         }
         // new indices are set -> Get new candidates
         if (fPriority[2]==1) {
         // we are entering the unique slice on this axis
         //---> intersect and store Y and X
            if (fPriority[1]==2) {
               if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>=fIby-1) return td.fVoxCheckList; // outside range
               ndd[0] = fNsliceY[td.fVoxSlices[1]];
               if (!ndd[0]) return td.fVoxCheckList;
               slice1 = &fIndcY[fOBy[td.fVoxSlices[1]]];
               islices++;
            }
            if (fPriority[0]==2) {
               if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>=fIbx-1) return td.fVoxCheckList; // outside range
               ndd[1] = fNsliceX[td.fVoxSlices[0]];
               if (!ndd[1]) return td.fVoxCheckList;
               islices++;
               if (slice1) {
                  slice2 = &fIndcX[fOBx[td.fVoxSlices[0]]];
               } else {
                  slice1 = &fIndcX[fOBx[td.fVoxSlices[0]]];
                  ndd[0] = ndd[1];
               }
            }
            if (islices<=1) {
               IntersectAndStore(ndd[0], slice1, td);
            } else {
               IntersectAndStore(ndd[0], slice1, ndd[1], slice2, td);
            }
            ncheck = td.fVoxNcandidates;
            return td.fVoxCheckList;
         }
         // We got into a new slice -> Get only new candidates
         left = (td.fVoxInc[2]>0)?kTRUE:kFALSE;
         new_list = GetExtraZ(td.fVoxSlices[2], left, ncheck);
//         printf("   New list on Z : %i new candidates\n", ncheck);
         if (!ncheck) return td.fVoxCheckList;
         if (fPriority[1]==2) {
            if (td.fVoxSlices[1]<0 || td.fVoxSlices[1]>=fIby-1) {
               ncheck = 0;
               return td.fVoxCheckList; // outside range
            }
            ndd[0] = fNsliceY[td.fVoxSlices[1]];
            if (!ndd[0]) {
               ncheck = 0;
               return td.fVoxCheckList;
            }
            slice1 = &fIndcY[fOBy[td.fVoxSlices[1]]];
            islices++;
         }
         if (fPriority[0]==2) {
            if (td.fVoxSlices[0]<0 || td.fVoxSlices[0]>=fIbx-1) {
               ncheck = 0;
               return td.fVoxCheckList; // outside range
            }
            ndd[1] = fNsliceX[td.fVoxSlices[0]];
            if (!ndd[1]) {
               ncheck = 0;
               return td.fVoxCheckList;
            }
            islices++;
            if (slice1) {
               slice2 = &fIndcX[fOBx[td.fVoxSlices[0]]];
            } else {
               slice1 = &fIndcX[fOBx[td.fVoxSlices[0]]];
               ndd[0] = ndd[1];
            }
         }
         if (!islices) return GetValidExtra(new_list, ncheck, td);
         if (islices==1) {
            return GetValidExtra(ndd[0], slice1, new_list, ncheck, td);
         } else {
            return GetValidExtra(ndd[0], slice1, ndd[1], slice2, new_list, ncheck, td);
         }
      default:
         Error("GetNextCandidates", "Invalid islice=%i inside %s", islice, fVolume->GetName());
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// get the list in the next voxel crossed by a ray

void TGeoVoxelFinder::SortCrossedVoxels(const Double_t *point, const Double_t *dir, TGeoStateInfo &td)
{
   if (NeedRebuild()) {
      TGeoVoxelFinder *vox = (TGeoVoxelFinder*)this;
      vox->Voxelize();
      fVolume->FindOverlaps();
   }
   td.fVoxCurrent = 0;
//   printf("###Sort crossed voxels for %s\n", fVolume->GetName());
   td.fVoxNcandidates = 0;
   Int_t  loc = 1+((fVolume->GetNdaughters()-1)>>3);
//   printf("   LOC=%i\n", loc*sizeof(UChar_t));
//   UChar_t *bits = gGeoManager->GetBits();
   memset(td.fVoxBits1, 0, loc);
   memset(td.fVoxInc, 0, 3*sizeof(Int_t));
   for (Int_t i=0; i<3; i++) {
      td.fVoxInvdir[i] = TGeoShape::Big();
      if (TMath::Abs(dir[i])<1E-10) continue;
      td.fVoxInc[i] = (dir[i]>0)?1:-1;
      td.fVoxInvdir[i] = 1./dir[i];
   }
   Bool_t flag = GetIndices(point, td);
   TGeoBBox *box = (TGeoBBox*)(fVolume->GetShape());
   const Double_t *box_orig = box->GetOrigin();
   if (td.fVoxInc[0]==0) {
      td.fVoxLimits[0] = TGeoShape::Big();
   } else {
      if (td.fVoxSlices[0]==-2) {
         // no slice on this axis -> get limit to bounding box limit
         td.fVoxLimits[0] = (box_orig[0]-point[0]+td.fVoxInc[0]*box->GetDX())*td.fVoxInvdir[0];
      } else {
         if (td.fVoxInc[0]==1) {
            td.fVoxLimits[0] = (fXb[fIbx-1]-point[0])*td.fVoxInvdir[0];
         } else {
            td.fVoxLimits[0] = (fXb[0]-point[0])*td.fVoxInvdir[0];
         }
      }
   }
   if (td.fVoxInc[1]==0) {
      td.fVoxLimits[1] = TGeoShape::Big();
   } else {
      if (td.fVoxSlices[1]==-2) {
         // no slice on this axis -> get limit to bounding box limit
         td.fVoxLimits[1] = (box_orig[1]-point[1]+td.fVoxInc[1]*box->GetDY())*td.fVoxInvdir[1];
      } else {
         if (td.fVoxInc[1]==1) {
            td.fVoxLimits[1] = (fYb[fIby-1]-point[1])*td.fVoxInvdir[1];
         } else {
            td.fVoxLimits[1] = (fYb[0]-point[1])*td.fVoxInvdir[1];
         }
      }
   }
   if (td.fVoxInc[2]==0) {
      td.fVoxLimits[2] = TGeoShape::Big();
   } else {
      if (td.fVoxSlices[2]==-2) {
         // no slice on this axis -> get limit to bounding box limit
         td.fVoxLimits[2] = (box_orig[2]-point[2]+td.fVoxInc[2]*box->GetDZ())*td.fVoxInvdir[2];
      } else {
         if (td.fVoxInc[2]==1) {
            td.fVoxLimits[2] = (fZb[fIbz-1]-point[2])*td.fVoxInvdir[2];
         } else {
            td.fVoxLimits[2] = (fZb[0]-point[2])*td.fVoxInvdir[2];
         }
      }
   }

   if (!flag) {
//      printf("   NO candidates in first voxel\n");
//      printf("   bits[0]=%i\n", bits[0]);
      return;
   }
//   printf("   current slices : %i   %i  %i\n", td.fVoxSlices[0], td.fVoxSlices[1], td.fVoxSlices[2]);
   Int_t nd[3];
   Int_t islices = 0;
   memset(&nd[0], 0, 3*sizeof(Int_t));
   UChar_t *slicex = 0;
   if (fPriority[0]==2) {
      nd[0] = fNsliceX[td.fVoxSlices[0]];
      slicex=&fIndcX[fOBx[td.fVoxSlices[0]]];
      islices++;
   }
   UChar_t *slicey = 0;
   if (fPriority[1]==2) {
      nd[1] = fNsliceY[td.fVoxSlices[1]];
      islices++;
      if (slicex) {
         slicey=&fIndcY[fOBy[td.fVoxSlices[1]]];
      } else {
         slicex=&fIndcY[fOBy[td.fVoxSlices[1]]];
         nd[0] = nd[1];
      }
   }
   UChar_t *slicez = 0;
   if (fPriority[2]==2) {
      nd[2] = fNsliceZ[td.fVoxSlices[2]];
      islices++;
      if (slicex && slicey) {
         slicez=&fIndcZ[fOBz[td.fVoxSlices[2]]];
      } else {
         if (slicex) {
            slicey=&fIndcZ[fOBz[td.fVoxSlices[2]]];
            nd[1] = nd[2];
         } else {
            slicex=&fIndcZ[fOBz[td.fVoxSlices[2]]];
            nd[0] = nd[2];
         }
      }
   }
//   printf("Ndaughters in first voxel : %i %i %i\n", nd[0], nd[1], nd[2]);
   switch (islices) {
      case 0:
         Error("SortCrossedVoxels", "no slices for %s", fVolume->GetName());
//         printf("Slices :(%i,%i,%i) Priority:(%i,%i,%i)\n", td.fVoxSlices[0], td.fVoxSlices[1], td.fVoxSlices[2], fPriority[0], fPriority[1], fPriority[2]);
         return;
      case 1:
         IntersectAndStore(nd[0], slicex, td);
         break;
      case 2:
         IntersectAndStore(nd[0], slicex, nd[1], slicey, td);
         break;
      default:
         IntersectAndStore(nd[0], slicex, nd[1], slicey, nd[2], slicez, td);
   }
//   printf("   bits[0]=%i  END\n", bits[0]);
//   if (td.fVoxNcandidates) {
//      printf("   candidates for first voxel :\n");
//      for (Int_t i=0; i<td.fVoxNcandidates; i++) printf("    %i\n", td.fVoxCheckList[i]);
//   }
}

////////////////////////////////////////////////////////////////////////////////
/// get the list of daughter indices for which point is inside their bbox

Int_t *TGeoVoxelFinder::GetCheckList(const Double_t *point, Int_t &nelem, TGeoStateInfo &td)
{
   if (NeedRebuild()) {
      Voxelize();
      fVolume->FindOverlaps();
   }
   if (fVolume->GetNdaughters() == 1) {
      if (fXb) {
         if (point[0]<fXb[0] || point[0]>fXb[1]) return 0;
      }
      if (fYb) {
         if (point[1]<fYb[0] || point[1]>fYb[1]) return 0;
      }

      if (fZb) {
         if (point[2]<fZb[0] || point[2]>fZb[1]) return 0;
      }
      td.fVoxCheckList[0] = 0;
      nelem = 1;
      return td.fVoxCheckList;
   }
   Int_t nslices = 0;
   UChar_t *slice1 = 0;
   UChar_t *slice2 = 0;
   UChar_t *slice3 = 0;
   Int_t nd[3] = {0,0,0};
   Int_t im;
   if (fPriority[0]) {
      im = TMath::BinarySearch(fIbx, fXb, point[0]);
      if ((im==-1) || (im==fIbx-1)) return 0;
      if (fPriority[0]==2) {
         nd[0] = fNsliceX[im];
         if (!nd[0]) return 0;
         nslices++;
         slice1 = &fIndcX[fOBx[im]];
      }
   }

   if (fPriority[1]) {
      im = TMath::BinarySearch(fIby, fYb, point[1]);
      if ((im==-1) || (im==fIby-1)) return 0;
      if (fPriority[1]==2) {
         nd[1] = fNsliceY[im];
         if (!nd[1]) return 0;
         nslices++;
         if (slice1) {
            slice2 = &fIndcY[fOBy[im]];
         } else {
            slice1 = &fIndcY[fOBy[im]];
            nd[0] = nd[1];
         }
      }
   }

   if (fPriority[2]) {
      im = TMath::BinarySearch(fIbz, fZb, point[2]);
      if ((im==-1) || (im==fIbz-1)) return 0;
      if (fPriority[2]==2) {
         nd[2] = fNsliceZ[im];
         if (!nd[2]) return 0;
         nslices++;
         if (slice1 && slice2) {
            slice3 = &fIndcZ[fOBz[im]];
         } else {
            if (slice1) {
               slice2 = &fIndcZ[fOBz[im]];
               nd[1] = nd[2];
            } else {
               slice1 = &fIndcZ[fOBz[im]];
               nd[0] = nd[2];
            }
         }
      }
   }
   nelem = 0;
//   Int_t i = 0;
   Bool_t intersect = kFALSE;
   switch (nslices) {
      case 0:
         Error("GetCheckList", "No slices for %s", fVolume->GetName());
         return 0;
      case 1:
         intersect = Intersect(nd[0], slice1, nelem, td.fVoxCheckList);
         break;
      case 2:
         intersect = Intersect(nd[0], slice1, nd[1], slice2, nelem, td.fVoxCheckList);
         break;
      default:
         intersect = Intersect(nd[0], slice1, nd[1], slice2, nd[2], slice3, nelem, td.fVoxCheckList);
   }
   if (intersect) return td.fVoxCheckList;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// get the list of candidates in voxel (i,j,k) - no check

Int_t *TGeoVoxelFinder::GetVoxelCandidates(Int_t i, Int_t j, Int_t k, Int_t &ncheck, TGeoStateInfo &td)
{
   UChar_t *slice1 = 0;
   UChar_t *slice2 = 0;
   UChar_t *slice3 = 0;
   Int_t nd[3] = {0,0,0};
   Int_t nslices = 0;
   if (fPriority[0]==2) {
      nd[0] = fNsliceX[i];
      if (!nd[0]) return 0;
      nslices++;
      slice1 = &fIndcX[fOBx[i]];
   }

   if (fPriority[1]==2) {
      nd[1] = fNsliceY[j];
      if (!nd[1]) return 0;
      nslices++;
      if (slice1) {
         slice2 = &fIndcY[fOBy[j]];
      } else {
         slice1 = &fIndcY[fOBy[j]];
         nd[0] = nd[1];
      }
   }

   if (fPriority[2]==2) {
      nd[2] = fNsliceZ[k];
      if (!nd[2]) return 0;
      nslices++;
      if (slice1 && slice2) {
         slice3 = &fIndcZ[fOBz[k]];
      } else {
         if (slice1) {
            slice2 = &fIndcZ[fOBz[k]];
            nd[1] = nd[2];
         } else {
            slice1 = &fIndcZ[fOBz[k]];
            nd[0] = nd[2];
         }
      }
   }
   Bool_t intersect = kFALSE;
   switch (nslices) {
      case 0:
         Error("GetCheckList", "No slices for %s", fVolume->GetName());
         return 0;
      case 1:
         intersect = Intersect(nd[0], slice1, ncheck, td.fVoxCheckList);
         break;
      case 2:
         intersect = Intersect(nd[0], slice1, nd[1], slice2, ncheck, td.fVoxCheckList);
         break;
      default:
         intersect = Intersect(nd[0], slice1, nd[1], slice2, nd[2], slice3, ncheck, td.fVoxCheckList);
   }
   if (intersect) return td.fVoxCheckList;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// get the list of new candidates for the next voxel crossed by current ray
///   printf("### GetNextVoxel\n");

Int_t *TGeoVoxelFinder::GetNextVoxel(const Double_t *point, const Double_t * /*dir*/, Int_t &ncheck, TGeoStateInfo &td)
{
   if (NeedRebuild()) {
      Voxelize();
      fVolume->FindOverlaps();
   }
   if (td.fVoxCurrent==0) {
//      printf(">>> first voxel, %i candidates\n", ncheck);
//      printf("   bits[0]=%i\n", gGeoManager->GetBits()[0]);
      td.fVoxCurrent++;
      ncheck = td.fVoxNcandidates;
      return td.fVoxCheckList;
   }
   td.fVoxCurrent++;
//   printf(">>> voxel %i\n", td.fCurrentVoxel);
   // Get slices for next voxel
//   printf("before - td.fSlices : %i %i %i\n", td.fSlices[0], td.fSlices[1], td.fSlices[2]);
   return GetNextCandidates(point, ncheck, td);
}

////////////////////////////////////////////////////////////////////////////////
/// return the list of nodes corresponding to one array of bits

Bool_t TGeoVoxelFinder::Intersect(Int_t n1, UChar_t *array1, Int_t &nf, Int_t *result)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
   nf = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   Bool_t ibreak = kFALSE;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = array1[current_byte];
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            result[nf++] = (current_byte<<3)+current_bit;
            if (nf==n1) {
               ibreak = kTRUE;
               break;
            }
         }
      }
      if (ibreak) return kTRUE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// return the list of nodes corresponding to one array of bits

Bool_t TGeoVoxelFinder::IntersectAndStore(Int_t n1, UChar_t *array1, TGeoStateInfo &td)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
//   UChar_t *bits = gGeoManager->GetBits();
   td.fVoxNcandidates = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   if (!array1) {
      memset(td.fVoxBits1, 0xFF, nbytes*sizeof(UChar_t));
      while (td.fVoxNcandidates<nd) {
         td.fVoxCheckList[td.fVoxNcandidates] = td.fVoxNcandidates;
         ++td.fVoxNcandidates;
      }
      return kTRUE;
   }
   memcpy(td.fVoxBits1, array1, nbytes*sizeof(UChar_t));
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   Bool_t ibreak = kFALSE;
   Int_t icand;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = array1[current_byte];
      icand = current_byte<<3;
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            td.fVoxCheckList[td.fVoxNcandidates++] = icand+current_bit;
            if (td.fVoxNcandidates==n1) {
               ibreak = kTRUE;
               break;
            }
         }
      }
      if (ibreak) return kTRUE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// make union of older bits with new array
///   printf("Union - one slice\n");

Bool_t TGeoVoxelFinder::Union(Int_t n1, UChar_t *array1, TGeoStateInfo &td)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
//   UChar_t *bits = gGeoManager->GetBits();
   td.fVoxNcandidates = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   Bool_t ibreak = kFALSE;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
//      printf("   byte %i : bits=%i array=%i\n", current_byte, bits[current_byte], array1[current_byte]);
      byte = (~td.fVoxBits1[current_byte]) & array1[current_byte];
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            td.fVoxCheckList[td.fVoxNcandidates++] = (current_byte<<3)+current_bit;
            if (td.fVoxNcandidates==n1) {
               ibreak = kTRUE;
               break;
            }
         }
      }
      td.fVoxBits1[current_byte] |= byte;
      if (ibreak) return kTRUE;
   }
   return (td.fVoxNcandidates>0);
}

////////////////////////////////////////////////////////////////////////////////
/// make union of older bits with new array
///   printf("Union - two slices\n");

Bool_t TGeoVoxelFinder::Union(Int_t /*n1*/, UChar_t *array1, Int_t /*n2*/, UChar_t *array2, TGeoStateInfo &td)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
//   UChar_t *bits = gGeoManager->GetBits();
   td.fVoxNcandidates = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = (~td.fVoxBits1[current_byte]) & (array1[current_byte] & array2[current_byte]);
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            td.fVoxCheckList[td.fVoxNcandidates++] = (current_byte<<3)+current_bit;
         }
      }
      td.fVoxBits1[current_byte] |= byte;
   }
   return (td.fVoxNcandidates>0);
}

////////////////////////////////////////////////////////////////////////////////
/// make union of older bits with new array
///   printf("Union - three slices\n");
///   printf("n1=%i n2=%i n3=%i\n", n1,n2,n3);

Bool_t TGeoVoxelFinder::Union(Int_t /*n1*/, UChar_t *array1, Int_t /*n2*/, UChar_t *array2, Int_t /*n3*/, UChar_t *array3, TGeoStateInfo &td)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
//   UChar_t *bits = gGeoManager->GetBits();
   td.fVoxNcandidates = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = (~td.fVoxBits1[current_byte]) & (array1[current_byte] & array2[current_byte] & array3[current_byte]);
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            td.fVoxCheckList[td.fVoxNcandidates++] = (current_byte<<3)+current_bit;
         }
      }
      td.fVoxBits1[current_byte] |= byte;
   }
   return (td.fVoxNcandidates>0);
}

////////////////////////////////////////////////////////////////////////////////
/// return the list of nodes corresponding to the intersection of two arrays of bits

Bool_t TGeoVoxelFinder::Intersect(Int_t n1, UChar_t *array1, Int_t n2, UChar_t *array2, Int_t &nf, Int_t *result)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
   nf = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   Bool_t ibreak = kFALSE;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = array1[current_byte] & array2[current_byte];
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            result[nf++] = (current_byte<<3)+current_bit;
            if ((nf==n1) || (nf==n2)) {
               ibreak = kTRUE;
               break;
            }
         }
      }
      if (ibreak) return kTRUE;
   }
   return (nf>0);
}

////////////////////////////////////////////////////////////////////////////////
/// return the list of nodes corresponding to the intersection of two arrays of bits

Bool_t TGeoVoxelFinder::IntersectAndStore(Int_t /*n1*/, UChar_t *array1, Int_t /*n2*/, UChar_t *array2, TGeoStateInfo &td)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
//   UChar_t *bits = gGeoManager->GetBits();
   td.fVoxNcandidates = 0;
   Int_t nbytes = 1+((nd-1)>>3);
//   memset(bits, 0, nbytes*sizeof(UChar_t));
   Int_t current_byte;
   Int_t current_bit;
   Int_t icand;
   UChar_t byte;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = array1[current_byte] & array2[current_byte];
      icand = current_byte<<3;
      td.fVoxBits1[current_byte] = byte;
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            td.fVoxCheckList[td.fVoxNcandidates++] = icand+current_bit;
         }
      }
   }
   return (td.fVoxNcandidates>0);
}

////////////////////////////////////////////////////////////////////////////////
/// return the list of nodes corresponding to the intersection of three arrays of bits

Bool_t TGeoVoxelFinder::Intersect(Int_t n1, UChar_t *array1, Int_t n2, UChar_t *array2, Int_t n3, UChar_t *array3, Int_t &nf, Int_t *result)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
   nf = 0;
   Int_t nbytes = 1+((nd-1)>>3);
   Int_t current_byte;
   Int_t current_bit;
   UChar_t byte;
   Bool_t ibreak = kFALSE;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = array1[current_byte] & array2[current_byte] & array3[current_byte];
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            result[nf++] = (current_byte<<3)+current_bit;
            if ((nf==n1) || (nf==n2) || (nf==n3)) {
               ibreak = kTRUE;
               break;
            }
         }
      }
      if (ibreak) return kTRUE;
   }
   return (nf>0);
}

////////////////////////////////////////////////////////////////////////////////
/// return the list of nodes corresponding to the intersection of three arrays of bits

Bool_t TGeoVoxelFinder::IntersectAndStore(Int_t /*n1*/, UChar_t *array1, Int_t /*n2*/, UChar_t *array2, Int_t /*n3*/, UChar_t *array3, TGeoStateInfo &td)
{
   Int_t nd = fVolume->GetNdaughters(); // also number of bits to scan
//   UChar_t *bits = gGeoManager->GetBits();
   td.fVoxNcandidates = 0;
   Int_t nbytes = 1+((nd-1)>>3);
//   memset(bits, 0, nbytes*sizeof(UChar_t));
   Int_t current_byte;
   Int_t current_bit;
   Int_t icand;
   UChar_t byte;
   for (current_byte=0; current_byte<nbytes; current_byte++) {
      byte = array1[current_byte] & array2[current_byte] & array3[current_byte];
      icand = current_byte<<3;
      td.fVoxBits1[current_byte] = byte;
      if (!byte) continue;
      for (current_bit=0; current_bit<8; current_bit++) {
         if (byte & (1<<current_bit)) {
            td.fVoxCheckList[td.fVoxNcandidates++] = icand+current_bit;
         }
      }
   }
   return (td.fVoxNcandidates>0);
}
////////////////////////////////////////////////////////////////////////////////
/// order bounding boxes along x, y, z

void TGeoVoxelFinder::SortAll(Option_t *)
{
   Int_t nd = fVolume->GetNdaughters();
   Int_t nperslice  = 1+(nd-1)/(8*sizeof(UChar_t)); /*Nbytes per slice*/
   Int_t nmaxslices = 2*nd+1; // max number of slices on each axis
   Double_t xmin, xmax, ymin, ymax, zmin, zmax;
   TGeoBBox *box = (TGeoBBox*)fVolume->GetShape(); // bounding box for volume
   // compute range on X, Y, Z according to volume bounding box
   xmin = (box->GetOrigin())[0] - box->GetDX();
   xmax = (box->GetOrigin())[0] + box->GetDX();
   ymin = (box->GetOrigin())[1] - box->GetDY();
   ymax = (box->GetOrigin())[1] + box->GetDY();
   zmin = (box->GetOrigin())[2] - box->GetDZ();
   zmax = (box->GetOrigin())[2] + box->GetDZ();
   if ((xmin>=xmax) || (ymin>=ymax) || (zmin>=zmax)) {
      Error("SortAll", "Wrong bounding box for volume %s", fVolume->GetName());
      return;
   }
   Int_t id;
   // compute boundaries coordinates on X,Y,Z
   Double_t *boundaries = new Double_t[6*nd]; // list of different boundaries
   for (id=0; id<nd; id++) {
      // x boundaries
      boundaries[2*id] = fBoxes[6*id+3]-fBoxes[6*id];
      boundaries[2*id+1] = fBoxes[6*id+3]+fBoxes[6*id];
      // y boundaries
      boundaries[2*id+2*nd] = fBoxes[6*id+4]-fBoxes[6*id+1];
      boundaries[2*id+2*nd+1] = fBoxes[6*id+4]+fBoxes[6*id+1];
      // z boundaries
      boundaries[2*id+4*nd] = fBoxes[6*id+5]-fBoxes[6*id+2];
      boundaries[2*id+4*nd+1] = fBoxes[6*id+5]+fBoxes[6*id+2];
   }
   Int_t *index = new Int_t[2*nd]; // indexes for sorted boundaries on one axis
   UChar_t *ind = new UChar_t[nmaxslices*nperslice]; // ind[fOBx[i]] = ndghts in slice fInd[i]--fInd[i+1]
   Int_t *extra = new Int_t[nmaxslices*4];
   // extra[fOEx[i]]   = nextra_to_left (i/i-1)
   // extra[fOEx[i]+1] = nextra_to_right (i/i+1)
   // Int_t *extra_to_left  = extra[fOEx[i]+2]
   // Int_t *extra_to_right = extra[fOEx[i]+2+nextra_to_left]
   Double_t *temp = new Double_t[2*nd]; // temporary array to store sorted boundary positions
   Int_t current  = 0;
   Int_t indextra = 0;
   Int_t nleft, nright;
   Int_t *extra_left  = new Int_t[nd];
   Int_t *extra_right = new Int_t[nd];
   Double_t xxmin, xxmax, xbmin, xbmax, ddx1, ddx2;
   UChar_t *bits;
   UInt_t loc, bitnumber;
   UChar_t bit;

   // sort x boundaries
   Int_t ib = 0; // total number of DIFFERENT boundaries
   TMath::Sort(2*nd, &boundaries[0], &index[0], kFALSE);
   // compact common boundaries
   for (id=0; id<2*nd; id++) {
      if (!ib) {temp[ib++] = boundaries[index[id]]; continue;}
      if (TMath::Abs(temp[ib-1]-boundaries[index[id]])>1E-10)
         temp[ib++] = boundaries[index[id]];
   }
   // now find priority
   if (ib < 2) {
      Error("SortAll", "Cannot voxelize %s :less than 2 boundaries on X", fVolume->GetName());
      delete [] boundaries;
      delete [] index;
      delete [] ind;
      delete [] temp;
      delete [] extra;
      delete [] extra_left;
      delete [] extra_right;
      SetInvalid();
      return;
   }
   if (ib == 2) {
   // check range
      if (((temp[0]-xmin)<1E-10) && ((temp[1]-xmax)>-1E-10)) {
      // ordering on this axis makes no sense. Clear all arrays.
         fPriority[0] = 0;  // always skip this axis
         if (fIndcX) delete [] fIndcX;
         fIndcX = 0;
         fNx = 0;
         if (fXb) delete [] fXb;
         fXb = 0;
         fIbx = 0;
         if (fOBx) delete [] fOBx;
         fOBx = 0;
         fNox = 0;
      } else {
         fPriority[0] = 1; // all in one slice
      }
   } else {
      fPriority[0] = 2;    // check all
   }
   // store compacted boundaries
   if (fPriority[0]) {
      if (fXb) delete [] fXb;
      fXb = new Double_t[ib];
      memcpy(fXb, &temp[0], ib*sizeof(Double_t));
      fIbx = ib;
   }

   //--- now build the lists of nodes in each slice of X axis
   if (fPriority[0]==2) {
      memset(ind, 0, (nmaxslices*nperslice)*sizeof(UChar_t));
      if (fOBx) delete [] fOBx;
      fNox = fIbx-1; // number of different slices
      fOBx = new Int_t[fNox]; // offsets in ind
      if (fOEx) delete [] fOEx;
      fOEx = new Int_t[fNox]; // offsets in extra
      if (fNsliceX) delete [] fNsliceX;
      fNsliceX = new Int_t[fNox];
      current  = 0;
      indextra = 0;
      //--- now loop all slices
      for (id=0; id<fNox; id++) {
         fOBx[id] = current; // offset in dght list for this slice
         fOEx[id] = indextra; // offset in exta list for this slice
         fNsliceX[id] = 0; // ndght in this slice
         extra[indextra] = extra[indextra+1] = 0; // no extra left/right
         nleft = nright = 0;
         bits = &ind[current]; // adress of bits for this slice
         xxmin = fXb[id];
         xxmax = fXb[id+1];
         for (Int_t ic=0; ic<nd; ic++) {
            xbmin = fBoxes[6*ic+3]-fBoxes[6*ic];
            xbmax = fBoxes[6*ic+3]+fBoxes[6*ic];
            ddx1 = xbmin-xxmax;
            if (ddx1>-1E-10) continue;
            ddx2 = xbmax-xxmin;
            if (ddx2<1E-10) continue;
            // daughter ic in interval
            //---> set the bit
            fNsliceX[id]++;
            bitnumber = (UInt_t)ic;
            loc = bitnumber/8;
            bit = bitnumber%8;
            bits[loc] |= 1<<bit;
            //---> chech if it is extra to left/right
            //--- left
            ddx1 = xbmin-xxmin;
            ddx2 = xbmax-xxmax;
            if ((id==0) || (ddx1>-1E-10)) {
               extra_left[nleft++] = ic;
            }
            //---right
            if ((id==(fNoz-1)) || (ddx2<1E-10)) {
               extra_right[nright++] = ic;
            }
         }
         //--- compute offset of next slice
         if (fNsliceX[id]>0) current += nperslice;
         //--- copy extra candidates
         extra[indextra] = nleft;
         extra[indextra+1] = nright;
         if (nleft)  memcpy(&extra[indextra+2], extra_left, nleft*sizeof(Int_t));
         if (nright) memcpy(&extra[indextra+2+nleft], extra_right, nright*sizeof(Int_t));
         indextra += 2+nleft+nright;
      }
      if (fIndcX) delete [] fIndcX;
      fNx = current;
      fIndcX = new UChar_t[current];
      memcpy(fIndcX, ind, current*sizeof(UChar_t));
      if (fExtraX) delete [] fExtraX;
      fNex = indextra;
      if (indextra>nmaxslices*4) printf("Woops!!!\n");
      fExtraX = new Int_t[indextra];
      memcpy(fExtraX, extra, indextra*sizeof(Int_t));
   }

   // sort y boundaries
   ib = 0;
   TMath::Sort(2*nd, &boundaries[2*nd], &index[0], kFALSE);
   // compact common boundaries
   for (id=0; id<2*nd; id++) {
      if (!ib) {temp[ib++] = boundaries[2*nd+index[id]]; continue;}
      if (TMath::Abs(temp[ib-1]-boundaries[2*nd+index[id]])>1E-10)
         temp[ib++]=boundaries[2*nd+index[id]];
   }
   // now find priority on Y
   if (ib < 2) {
      Error("SortAll", "Cannot voxelize %s :less than 2 boundaries on Y", fVolume->GetName());
      delete [] boundaries;
      delete [] index;
      delete [] ind;
      delete [] temp;
      delete [] extra;
      delete [] extra_left;
      delete [] extra_right;
      SetInvalid();
      return;
   }
   if (ib == 2) {
   // check range
      if (((temp[0]-ymin)<1E-10) && ((temp[1]-ymax)>-1E-10)) {
      // ordering on this axis makes no sense. Clear all arrays.
         fPriority[1] = 0; // always skip this axis
         if (fIndcY) delete [] fIndcY;
         fIndcY = 0;
         fNy = 0;
         if (fYb) delete [] fYb;
         fYb = 0;
         fIby = 0;
         if (fOBy) delete [] fOBy;
         fOBy = 0;
         fNoy = 0;
      } else {
         fPriority[1] = 1; // all in one slice
      }
   } else {
      fPriority[1] = 2;    // check all
   }

   // store compacted boundaries
   if (fPriority[1]) {
      if (fYb) delete [] fYb;
      fYb = new Double_t[ib];
      memcpy(fYb, &temp[0], ib*sizeof(Double_t));
      fIby = ib;
   }


   if (fPriority[1]==2) {
   //--- now build the lists of nodes in each slice of Y axis
      memset(ind, 0, (nmaxslices*nperslice)*sizeof(UChar_t));
      if (fOBy) delete [] fOBy;
      fNoy = fIby-1; // number of different slices
      fOBy = new Int_t[fNoy]; // offsets in ind
      if (fOEy) delete [] fOEy;
      fOEy = new Int_t[fNoy]; // offsets in extra
      if (fNsliceY) delete [] fNsliceY;
      fNsliceY = new Int_t[fNoy];
      current = 0;
      indextra = 0;
      //--- now loop all slices
      for (id=0; id<fNoy; id++) {
         fOBy[id] = current; // offset of dght list
         fOEy[id] = indextra; // offset in exta list for this slice
         fNsliceY[id] = 0; // ndght in this slice
         extra[indextra] = extra[indextra+1] = 0; // no extra left/right
         nleft = nright = 0;
         bits = &ind[current]; // adress of bits for this slice
         xxmin = fYb[id];
         xxmax = fYb[id+1];
         for (Int_t ic=0; ic<nd; ic++) {
            xbmin = fBoxes[6*ic+4]-fBoxes[6*ic+1];
            xbmax = fBoxes[6*ic+4]+fBoxes[6*ic+1];
            ddx1 = xbmin-xxmax;
            if (ddx1>-1E-10) continue;
            ddx2 = xbmax-xxmin;
            if (ddx2<1E-10) continue;
            // daughter ic in interval
            //---> set the bit
            fNsliceY[id]++;
            bitnumber = (UInt_t)ic;
            loc = bitnumber/8;
            bit = bitnumber%8;
            bits[loc] |= 1<<bit;
            //---> check if it is extra to left/right
            //--- left
            ddx1 = xbmin-xxmin;
            ddx2 = xbmax-xxmax;
            if ((id==0) || (ddx1>-1E-10)) {
               extra_left[nleft++] = ic;
            }
            //---right
            if ((id==(fNoz-1)) || (ddx2<1E-10)) {
               extra_right[nright++] = ic;
            }
         }
         //--- compute offset of next slice
         if (fNsliceY[id]>0) current += nperslice;
         //--- copy extra candidates
         extra[indextra] = nleft;
         extra[indextra+1] = nright;
         if (nleft)  memcpy(&extra[indextra+2], extra_left, nleft*sizeof(Int_t));
         if (nright) memcpy(&extra[indextra+2+nleft], extra_right, nright*sizeof(Int_t));
         indextra += 2+nleft+nright;
      }
      if (fIndcY) delete [] fIndcY;
      fNy = current;
      fIndcY = new UChar_t[current];
      memcpy(fIndcY, &ind[0], current*sizeof(UChar_t));
      if (fExtraY) delete [] fExtraY;
      fNey = indextra;
      if (indextra>nmaxslices*4) printf("Woops!!!\n");
      fExtraY = new Int_t[indextra];
      memcpy(fExtraY, extra, indextra*sizeof(Int_t));
   }

   // sort z boundaries
   ib = 0;
   TMath::Sort(2*nd, &boundaries[4*nd], &index[0], kFALSE);
   // compact common boundaries
   for (id=0; id<2*nd; id++) {
      if (!ib) {temp[ib++] = boundaries[4*nd+index[id]]; continue;}
      if ((TMath::Abs(temp[ib-1]-boundaries[4*nd+index[id]]))>1E-10)
         temp[ib++]=boundaries[4*nd+index[id]];
   }
   // now find priority on Z
   if (ib < 2) {
      Error("SortAll", "Cannot voxelize %s :less than 2 boundaries on Z", fVolume->GetName());
      delete [] boundaries;
      delete [] index;
      delete [] ind;
      delete [] temp;
      delete [] extra;
      delete [] extra_left;
      delete [] extra_right;
      SetInvalid();
      return;
   }
   if (ib == 2) {
   // check range
      if (((temp[0]-zmin)<1E-10) && ((temp[1]-zmax)>-1E-10)) {
      // ordering on this axis makes no sense. Clear all arrays.
         fPriority[2] = 0;
         if (fIndcZ) delete [] fIndcZ;
         fIndcZ = 0;
         fNz = 0;
         if (fZb) delete [] fZb;
         fZb = 0;
         fIbz = 0;
         if (fOBz) delete [] fOBz;
         fOBz = 0;
         fNoz = 0;
      } else {
         fPriority[2] = 1; // all in one slice
      }
   } else {
      fPriority[2] = 2;    // check all
   }

   // store compacted boundaries
   if (fPriority[2]) {
      if (fZb) delete [] fZb;
      fZb = new Double_t[ib];
      memcpy(fZb, &temp[0], ib*sizeof(Double_t));
      fIbz = ib;
   }


   if (fPriority[2]==2) {
   //--- now build the lists of nodes in each slice of Y axis
      memset(ind, 0, (nmaxslices*nperslice)*sizeof(UChar_t));
      if (fOBz) delete [] fOBz;
      fNoz = fIbz-1; // number of different slices
      fOBz = new Int_t[fNoz]; // offsets in ind
      if (fOEz) delete [] fOEz;
      fOEz = new Int_t[fNoz]; // offsets in extra
      if (fNsliceZ) delete [] fNsliceZ;
      fNsliceZ = new Int_t[fNoz];
      current = 0;
      indextra = 0;
      //--- now loop all slices
      for (id=0; id<fNoz; id++) {
         fOBz[id] = current; // offset of dght list
         fOEz[id] = indextra; // offset in exta list for this slice
         fNsliceZ[id] = 0; // ndght in this slice
         extra[indextra] = extra[indextra+1] = 0; // no extra left/right
         nleft = nright = 0;
         bits = &ind[current]; // adress of bits for this slice
         xxmin = fZb[id];
         xxmax = fZb[id+1];
         for (Int_t ic=0; ic<nd; ic++) {
            xbmin = fBoxes[6*ic+5]-fBoxes[6*ic+2];
            xbmax = fBoxes[6*ic+5]+fBoxes[6*ic+2];
            ddx1 = xbmin-xxmax;
            if (ddx1>-1E-10) continue;
            ddx2 = xbmax-xxmin;
            if (ddx2<1E-10) continue;
            // daughter ic in interval
            //---> set the bit
            fNsliceZ[id]++;
            bitnumber = (UInt_t)ic;
            loc = bitnumber/8;
            bit = bitnumber%8;
            bits[loc] |= 1<<bit;
            //---> check if it is extra to left/right
            //--- left
            ddx1 = xbmin-xxmin;
            ddx2 = xbmax-xxmax;
            if ((id==0) || (ddx1>-1E-10)) {
               extra_left[nleft++] = ic;
            }
            //---right
            if ((id==(fNoz-1)) || (ddx2<1E-10)) {
               extra_right[nright++] = ic;
            }
         }
         //--- compute offset of next slice
         if (fNsliceZ[id]>0) current += nperslice;
         //--- copy extra candidates
         extra[indextra] = nleft;
         extra[indextra+1] = nright;
         if (nleft)  memcpy(&extra[indextra+2], extra_left, nleft*sizeof(Int_t));
         if (nright) memcpy(&extra[indextra+2+nleft], extra_right, nright*sizeof(Int_t));
         indextra += 2+nleft+nright;
      }
      if (fIndcZ) delete [] fIndcZ;
      fNz = current;
      fIndcZ = new UChar_t[current];
      memcpy(fIndcZ, &ind[0], current*sizeof(UChar_t));
      if (fExtraZ) delete [] fExtraZ;
      fNez = indextra;
      if (indextra>nmaxslices*4) printf("Woops!!!\n");
      fExtraZ = new Int_t[indextra];
      memcpy(fExtraZ, extra, indextra*sizeof(Int_t));
   }
   delete [] boundaries;
   delete [] index;
   delete [] temp;
   delete [] ind;
   delete [] extra;
   delete [] extra_left;
   delete [] extra_right;

   if ((!fPriority[0]) && (!fPriority[1]) && (!fPriority[2])) {
      SetInvalid();
      if (nd>1) Error("SortAll", "Volume %s: Cannot make slices on any axis", fVolume->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print the voxels.

void TGeoVoxelFinder::Print(Option_t *) const
{
   if (NeedRebuild()) {
      TGeoVoxelFinder *vox = (TGeoVoxelFinder*)this;
      vox->Voxelize();
      fVolume->FindOverlaps();
   }
   Int_t id, i;
   Int_t nd = fVolume->GetNdaughters();
   printf("Voxels for volume %s (nd=%i)\n", fVolume->GetName(), fVolume->GetNdaughters());
   printf("priority : x=%i y=%i z=%i\n", fPriority[0], fPriority[1], fPriority[2]);
//   return;
   Int_t nextra;
   Int_t nbytes = 1+((fVolume->GetNdaughters()-1)>>3);
   UChar_t byte, bit;
   UChar_t *slice;
   printf("XXX\n");
   if (fPriority[0]==2) {
      for (id=0; id<fIbx; id++) {
         printf("%15.10f\n",fXb[id]);
         if (id == (fIbx-1)) break;
         printf("slice %i : %i\n", id, fNsliceX[id]);
         if (fNsliceX[id]) {
            slice = &fIndcX[fOBx[id]];
            for (i=0; i<nbytes; i++) {
               byte = slice[i];
               for (bit=0; bit<8; bit++) {
                  if (byte & (1<<bit)) printf(" %i ", 8*i+bit);
               }
            }
            printf("\n");
         }
         GetExtraX(id,kTRUE,nextra);
         printf("   extra_about_left  = %i\n", nextra);
         GetExtraX(id,kFALSE,nextra);
         printf("   extra_about_right = %i\n", nextra);
      }
   } else if (fPriority[0]==1) {
      printf("%15.10f\n",fXb[0]);
      for (id=0; id<nd; id++) printf(" %i ",id);
      printf("\n");
      printf("%15.10f\n",fXb[1]);
   }
   printf("YYY\n");
   if (fPriority[1]==2) {
      for (id=0; id<fIby; id++) {
         printf("%15.10f\n", fYb[id]);
         if (id == (fIby-1)) break;
         printf("slice %i : %i\n", id, fNsliceY[id]);
         if (fNsliceY[id]) {
            slice = &fIndcY[fOBy[id]];
            for (i=0; i<nbytes; i++) {
               byte = slice[i];
               for (bit=0; bit<8; bit++) {
                  if (byte & (1<<bit)) printf(" %i ", 8*i+bit);
               }
            }
         }
         GetExtraY(id,kTRUE,nextra);
         printf("   extra_about_left  = %i\n", nextra);
         GetExtraY(id,kFALSE,nextra);
         printf("   extra_about_right = %i\n", nextra);
      }
   } else if (fPriority[1]==1) {
      printf("%15.10f\n",fYb[0]);
      for (id=0; id<nd; id++) printf(" %i ",id);
      printf("\n");
      printf("%15.10f\n",fYb[1]);
   }

   printf("ZZZ\n");
   if (fPriority[2]==2) {
      for (id=0; id<fIbz; id++) {
         printf("%15.10f\n", fZb[id]);
         if (id == (fIbz-1)) break;
         printf("slice %i : %i\n", id, fNsliceZ[id]);
         if (fNsliceZ[id]) {
            slice = &fIndcZ[fOBz[id]];
            for (i=0; i<nbytes; i++) {
               byte = slice[i];
               for (bit=0; bit<8; bit++) {
                  if (byte & (1<<bit)) printf(" %i ", 8*i+bit);
               }
            }
            printf("\n");
         }
         GetExtraZ(id,kTRUE,nextra);
         printf("   extra_about_left  = %i\n", nextra);
         GetExtraZ(id,kFALSE,nextra);
         printf("   extra_about_right = %i\n", nextra);
      }
   } else if (fPriority[2]==1) {
      printf("%15.10f\n",fZb[0]);
      for (id=0; id<nd; id++) printf(" %i ",id);
      printf("\n");
      printf("%15.10f\n",fZb[1]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print the voxel containing point

void TGeoVoxelFinder::PrintVoxelLimits(const Double_t *point) const
{
   if (NeedRebuild()) {
      TGeoVoxelFinder *vox = (TGeoVoxelFinder*)this;
      vox->Voxelize();
      fVolume->FindOverlaps();
   }
   Int_t im=0;
   if (fPriority[0]) {
      im = TMath::BinarySearch(fIbx, fXb, point[0]);
      if ((im==-1) || (im==fIbx-1)) {
         printf("Voxel X limits: OUT\n");
      } else {
         printf("Voxel X limits: %g  %g\n", fXb[im], fXb[im+1]);
      }
   }
   if (fPriority[1]) {
      im = TMath::BinarySearch(fIby, fYb, point[1]);
      if ((im==-1) || (im==fIby-1)) {
         printf("Voxel Y limits: OUT\n");
      } else {
         printf("Voxel Y limits: %g  %g\n", fYb[im], fYb[im+1]);
      }
   }
   if (fPriority[2]) {
      im = TMath::BinarySearch(fIbz, fZb, point[2]);
      if ((im==-1) || (im==fIbz-1)) {
         printf("Voxel Z limits: OUT\n");
      } else {
         printf("Voxel Z limits: %g  %g\n", fZb[im], fZb[im+1]);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Voxelize attached volume according to option
/// If the volume is an assembly, make sure the bbox is computed.

void TGeoVoxelFinder::Voxelize(Option_t * /*option*/)
{
   if (fVolume->IsAssembly()) fVolume->GetShape()->ComputeBBox();
   Int_t nd = fVolume->GetNdaughters();
   TGeoVolume *vd;
   for (Int_t i=0; i<nd; i++) {
      vd = fVolume->GetNode(i)->GetVolume();
      if (vd->IsAssembly()) vd->GetShape()->ComputeBBox();
   }
   BuildVoxelLimits();
   SortAll();
   SetNeedRebuild(kFALSE);
}
////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGeoVoxelFinder.

void TGeoVoxelFinder::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TGeoVoxelFinder::Class(), this, R__v, R__s, R__c);
         return;
      }
      // Process old versions of the voxel finder. Just read the data
      // from the buffer in a temp variable then mark voxels as garbage.
      UChar_t *dummy = new UChar_t[R__c-2];
      R__b.ReadFastArray(dummy, R__c-2);
      delete [] dummy;
      SetInvalid(kTRUE);
   } else {
      R__b.WriteClassBuffer(TGeoVoxelFinder::Class(), this);
   }
}
