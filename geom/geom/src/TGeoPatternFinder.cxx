// @(#)root/geom:$Id$
// Author: Andrei Gheata   30/10/01

/** \class TGeoPatternFinder
\ingroup Geometry_classes

Base finder class for patterns.

  A pattern is specifying a division type which applies only to a given
shape type. The implemented patterns are for the moment equidistant slices
on different axis. Implemented patterns are:

  - TGeoPatternX - a X axis divison pattern
  - TGeoPatternY - a Y axis divison pattern
  - TGeoPatternZ - a Z axis divison pattern
  - TGeoPatternParaX - a X axis divison pattern for PARA shape
  - TGeoPatternParaY - a Y axis divison pattern for PARA shape
  - TGeoPatternParaZ - a Z axis divison pattern for PARA shape
  - TGeoPatternTrapZ - a Z axis divison pattern for TRAP or GTRA shapes
  - TGeoPatternCylR - a cylindrical R divison pattern
  - TGeoPatternCylPhi - a cylindrical phi divison pattern
  - TGeoPatternSphR - a spherical R divison pattern
  - TGeoPatternSphTheta - a spherical theta divison pattern
  - TGeoPatternSphPhi - a spherical phi divison pattern
  - TGeoPatternHoneycomb - a divison pattern specialized for honeycombs
*/

#include "TGeoPatternFinder.h"

#include "TBuffer.h"
#include "TObject.h"
#include "TGeoMatrix.h"
#include "TGeoPara.h"
#include "TGeoArb8.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TMath.h"

ClassImp(TGeoPatternFinder);
ClassImp(TGeoPatternX);
ClassImp(TGeoPatternY);
ClassImp(TGeoPatternZ);
ClassImp(TGeoPatternParaX);
ClassImp(TGeoPatternParaY);
ClassImp(TGeoPatternParaZ);
ClassImp(TGeoPatternTrapZ);
ClassImp(TGeoPatternCylR);
ClassImp(TGeoPatternCylPhi);
ClassImp(TGeoPatternSphR);
ClassImp(TGeoPatternSphTheta);
ClassImp(TGeoPatternSphPhi);
ClassImp(TGeoPatternHoneycomb);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoPatternFinder::ThreadData_t::ThreadData_t() :
   fMatrix(0), fCurrent(-1), fNextIndex(-1)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoPatternFinder::ThreadData_t::~ThreadData_t()
{
//   if (fMatrix != gGeoIdentity) delete fMatrix;
}

////////////////////////////////////////////////////////////////////////////////

TGeoPatternFinder::ThreadData_t& TGeoPatternFinder::GetThreadData() const
{
   Int_t tid = TGeoManager::ThreadId();
   return *fThreadData[tid];
}

////////////////////////////////////////////////////////////////////////////////

void TGeoPatternFinder::ClearThreadData() const
{
   std::lock_guard<std::mutex> guard(fMutex);
   std::vector<ThreadData_t*>::iterator i = fThreadData.begin();
   while (i != fThreadData.end())
   {
      delete *i;
      ++i;
   }
   fThreadData.clear();
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create thread data for n threads max.

void TGeoPatternFinder::CreateThreadData(Int_t nthreads)
{
   std::lock_guard<std::mutex> guard(fMutex);
   fThreadData.resize(nthreads);
   fThreadSize = nthreads;
   for (Int_t tid=0; tid<nthreads; tid++) {
      if (fThreadData[tid] == 0) {
         fThreadData[tid] = new ThreadData_t;
         fThreadData[tid]->fMatrix = CreateMatrix();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternFinder::TGeoPatternFinder()
{
   fNdivisions = 0;
   fDivIndex   = 0;
   fStep       = 0;
   fStart      = 0;
   fEnd        = 0;
   fVolume     = 0;
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternFinder::TGeoPatternFinder(TGeoVolume *vol, Int_t ndiv)
{
   fVolume     = vol;
   fNdivisions = ndiv;
   fDivIndex   = 0;
   fStep       = 0;
   fStart      = 0;
   fEnd        = 0;
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternFinder::TGeoPatternFinder(const TGeoPatternFinder& pf) :
  TObject(pf),
  fStep(pf.fStep),
  fStart(pf.fStart),
  fEnd(pf.fEnd),
  fNdivisions(pf.fNdivisions),
  fDivIndex(pf.fDivIndex),
  fVolume(pf.fVolume)
{
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternFinder& TGeoPatternFinder::operator=(const TGeoPatternFinder& pf)
{
   if(this!=&pf) {
      TObject::operator=(pf);
      fStep=pf.fStep;
      fStart=pf.fStart;
      fEnd=pf.fEnd;
      fNdivisions=pf.fNdivisions;
      fDivIndex=pf.fDivIndex;
      fVolume=pf.fVolume;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternFinder::~TGeoPatternFinder()
{
   ClearThreadData();
}

////////////////////////////////////////////////////////////////////////////////
/// Return current index.

Int_t TGeoPatternFinder::GetCurrent()
{
   return GetThreadData().fCurrent;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current matrix.

TGeoMatrix* TGeoPatternFinder::GetMatrix()
{
   return GetThreadData().fMatrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Get index of next division.

Int_t TGeoPatternFinder::GetNext() const
{
   return GetThreadData().fNextIndex;
}

////////////////////////////////////////////////////////////////////////////////
/// Set index of next division.

void TGeoPatternFinder::SetNext(Int_t index)
{
   GetThreadData().fNextIndex = index;
}

////////////////////////////////////////////////////////////////////////////////
/// Make next node (if any) current.

TGeoNode *TGeoPatternFinder::CdNext()
{
   ThreadData_t& td = GetThreadData();
   if (td.fNextIndex < 0) return NULL;
   cd(td.fNextIndex);
   return GetNodeOffset(td.fCurrent);
}

////////////////////////////////////////////////////////////////////////////////
/// Set division range. Use this method only when dividing an assembly.

void TGeoPatternFinder::SetRange(Double_t start, Double_t step, Int_t ndivisions)
{
   fStart = start;
   fEnd = fStart + ndivisions*step;
   fStep = step;
   fNdivisions = ndivisions;
}

//______________________________________________________________________________
// TGeoPatternX - a X axis divison pattern
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternX::TGeoPatternX()
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternX::TGeoPatternX(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dx = ((TGeoBBox*)vol->GetShape())->GetDX();
   fStart = -dx;
   fEnd = dx;
   fStep = 2*dx/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternX::TGeoPatternX(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dx = ((TGeoBBox*)vol->GetShape())->GetDX();
   fStart = -dx;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternX::TGeoPatternX(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternX::TGeoPatternX(const TGeoPatternX& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternX& TGeoPatternX::operator=(const TGeoPatternX& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternX::~TGeoPatternX()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternX::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
   td.fMatrix->SetDx(fStart+idiv*fStep+0.5*fStep);
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternX::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternX::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   matrix.SetDx(fStart+idiv*fStep+0.5*fStep);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternX::IsOnBoundary(const Double_t *point) const
{
   Double_t seg = (point[0]-fStart)/fStep;
   Double_t diff = seg - Long64_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the cell corresponding to point and next cell along dir (if asked)

TGeoNode *TGeoPatternX::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Int_t ind = (Int_t)(1.+(point[0]-fStart)/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      if (dir[0]>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance to next division layer returning the index of next section.
/// Point is in the frame of the divided volume.

Double_t TGeoPatternX::FindNextBoundary(Double_t *point, Double_t *dir, Int_t &indnext)
{
   ThreadData_t& td = GetThreadData();
   indnext = -1;
   Double_t dist = TGeoShape::Big();
   if (TMath::Abs(dir[0])<TGeoShape::Tolerance()) return dist;
   if (td.fCurrent<0) {
      Error("FindNextBoundary", "Must call FindNode first");
      return dist;
   }
   Int_t inc = (dir[0]>0)?1:0;
   dist = (fStep*(td.fCurrent+inc)-point[0])/dir[0];
   if (dist<0.) Error("FindNextBoundary", "Negative distance d=%g",dist);
   if (!inc) inc = -1;
   indnext = td.fCurrent+inc;
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternX::MakeCopy(Bool_t reflect)
{
   TGeoPatternX *finder = new TGeoPatternX(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternX::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 1;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

//______________________________________________________________________________
// TGeoPatternY - a Y axis divison pattern
//______________________________________________________________________________


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternY::TGeoPatternY()
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternY::TGeoPatternY(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dy = ((TGeoBBox*)vol->GetShape())->GetDY();
   fStart = -dy;
   fEnd = dy;
   fStep = 2*dy/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternY::TGeoPatternY(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dy = ((TGeoBBox*)vol->GetShape())->GetDY();
   fStart = -dy;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternY::TGeoPatternY(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternY::TGeoPatternY(const TGeoPatternY& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternY& TGeoPatternY::operator=(const TGeoPatternY& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternY::~TGeoPatternY()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternY::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
   td.fMatrix->SetDy(fStart+idiv*fStep+0.5*fStep);
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternY::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternY::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   matrix.SetDy(fStart+idiv*fStep+0.5*fStep);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternY::IsOnBoundary(const Double_t *point) const
{
   Double_t seg = (point[1]-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the cell corresponding to point and next cell along dir (if asked)

TGeoNode *TGeoPatternY::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Int_t ind = (Int_t)(1.+(point[1]-fStart)/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      if (dir[1]>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance to next division layer returning the index of next section.
/// Point is in the frame of the divided volume.

Double_t TGeoPatternY::FindNextBoundary(Double_t *point, Double_t *dir, Int_t &indnext)
{
   ThreadData_t& td = GetThreadData();
   indnext = -1;
   Double_t dist = TGeoShape::Big();
   if (TMath::Abs(dir[1])<TGeoShape::Tolerance()) return dist;
   if (td.fCurrent<0) {
      Error("FindNextBoundary", "Must call FindNode first");
      return dist;
   }
   Int_t inc = (dir[1]>0)?1:0;
   dist = (fStep*(td.fCurrent+inc)-point[1])/dir[1];
   if (dist<0.) Error("FindNextBoundary", "Negative distance d=%g",dist);
   if (!inc) inc = -1;
   indnext = td.fCurrent+inc;
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternY::MakeCopy(Bool_t reflect)
{
   TGeoPatternY *finder = new TGeoPatternY(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternY::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 2;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

//______________________________________________________________________________
// TGeoPatternZ - a Z axis divison pattern
//______________________________________________________________________________


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternZ::TGeoPatternZ()
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternZ::TGeoPatternZ(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dz = ((TGeoBBox*)vol->GetShape())->GetDZ();
   fStart = -dz;
   fEnd = dz;
   fStep = 2*dz/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternZ::TGeoPatternZ(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dz = ((TGeoBBox*)vol->GetShape())->GetDZ();
   fStart = -dz;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternZ::TGeoPatternZ(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternZ::TGeoPatternZ(const TGeoPatternZ& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternZ& TGeoPatternZ::operator=(const TGeoPatternZ& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternZ::~TGeoPatternZ()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternZ::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
   td.fMatrix->SetDz(((IsReflected())?-1.:1.)*(fStart+idiv*fStep+0.5*fStep));
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternZ::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternZ::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   matrix.SetDz(((IsReflected())?-1.:1.)*(fStart+idiv*fStep+0.5*fStep));
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternZ::IsOnBoundary(const Double_t *point) const
{
   Double_t seg = (point[2]-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the cell corresponding to point and next cell along dir (if asked)

TGeoNode *TGeoPatternZ::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Int_t ind = (Int_t)(1.+(point[2]-fStart)/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      if (dir[2]>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance to next division layer returning the index of next section.
/// Point is in the frame of the divided volume.

Double_t TGeoPatternZ::FindNextBoundary(Double_t *point, Double_t *dir, Int_t &indnext)
{
   indnext = -1;
   ThreadData_t& td = GetThreadData();
   Double_t dist = TGeoShape::Big();
   if (TMath::Abs(dir[2])<TGeoShape::Tolerance()) return dist;
   if (td.fCurrent<0) {
      Error("FindNextBoundary", "Must call FindNode first");
      return dist;
   }
   Int_t inc = (dir[2]>0)?1:0;
   dist = (fStep*(td.fCurrent+inc)-point[2])/dir[2];
   if (dist<0.) Error("FindNextBoundary", "Negative distance d=%g",dist);
   if (!inc) inc = -1;
   indnext = td.fCurrent+inc;
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternZ::MakeCopy(Bool_t reflect)
{
   TGeoPatternZ *finder = new TGeoPatternZ(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternZ::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 3;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

//______________________________________________________________________________
// TGeoPatternParaX - a X axis divison pattern for PARA shape
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternParaX::TGeoPatternParaX()
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaX::TGeoPatternParaX(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dx = ((TGeoPara*)vol->GetShape())->GetX();
   fStart = -dx;
   fEnd = dx;
   fStep = 2*dx/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaX::TGeoPatternParaX(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t dx = ((TGeoPara*)vol->GetShape())->GetX();
   fStart = -dx;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaX::TGeoPatternParaX(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternParaX::TGeoPatternParaX(const TGeoPatternParaX& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternParaX& TGeoPatternParaX::operator=(const TGeoPatternParaX& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternParaX::~TGeoPatternParaX()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternParaX::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
   td.fMatrix->SetDx(fStart+idiv*fStep+0.5*fStep);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternParaX::IsOnBoundary(const Double_t *point) const
{
   Double_t txy = ((TGeoPara*)fVolume->GetShape())->GetTxy();
   Double_t txz = ((TGeoPara*)fVolume->GetShape())->GetTxz();
   Double_t tyz = ((TGeoPara*)fVolume->GetShape())->GetTyz();
   Double_t xt = point[0]-txz*point[2]-txy*(point[1]-tyz*point[2]);
   Double_t seg = (xt-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// get the node division containing the query point

TGeoNode *TGeoPatternParaX::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Double_t txy = ((TGeoPara*)fVolume->GetShape())->GetTxy();
   Double_t txz = ((TGeoPara*)fVolume->GetShape())->GetTxz();
   Double_t tyz = ((TGeoPara*)fVolume->GetShape())->GetTyz();
   Double_t xt = point[0]-txz*point[2]-txy*(point[1]-tyz*point[2]);
   Int_t ind = (Int_t)(1.+(xt-fStart)/fStep)-1;
   if (dir) {
      Double_t ttsq = txy*txy + (txz-txy*tyz)*(txz-txy*tyz);
      Double_t divdirx = 1./TMath::Sqrt(1.+ttsq);
      Double_t divdiry = -txy*divdirx;
      Double_t divdirz = -(txz-txy*tyz)*divdirx;
      Double_t dot = dir[0]*divdirx + dir[1]*divdiry + dir[2]*divdirz;
      td.fNextIndex = ind;
      if (dot>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternParaX::MakeCopy(Bool_t reflect)
{
   TGeoPatternParaX *finder = new TGeoPatternParaX(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternParaX::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 1;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternParaX::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternParaX::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   matrix.SetDx(fStart+idiv*fStep+0.5*fStep);
}

//______________________________________________________________________________
// TGeoPatternParaY - a Y axis divison pattern for PARA shape
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternParaY::TGeoPatternParaY()
{
   fTxy = 0;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaY::TGeoPatternParaY(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   fTxy = ((TGeoPara*)vol->GetShape())->GetTxy();
   Double_t dy = ((TGeoPara*)vol->GetShape())->GetY();
   fStart = -dy;
   fEnd = dy;
   fStep = 2*dy/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaY::TGeoPatternParaY(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   fTxy = ((TGeoPara*)vol->GetShape())->GetTxy();
   Double_t dy = ((TGeoPara*)vol->GetShape())->GetY();
   fStart = -dy;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaY::TGeoPatternParaY(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   fTxy = ((TGeoPara*)vol->GetShape())->GetTxy();
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternParaY::TGeoPatternParaY(const TGeoPatternParaY& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternParaY& TGeoPatternParaY::operator=(const TGeoPatternParaY& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternParaY::~TGeoPatternParaY()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternParaY::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent = idiv;
   Double_t dy = fStart+idiv*fStep+0.5*fStep;
   td.fMatrix->SetDx(fTxy*dy);
   td.fMatrix->SetDy(dy);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternParaY::IsOnBoundary(const Double_t *point) const
{
   Double_t tyz = ((TGeoPara*)fVolume->GetShape())->GetTyz();
   Double_t yt = point[1]-tyz*point[2];
   Double_t seg = (yt-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// get the node division containing the query point

TGeoNode *TGeoPatternParaY::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Double_t tyz = ((TGeoPara*)fVolume->GetShape())->GetTyz();
   Double_t yt = point[1]-tyz*point[2];
   Int_t ind = (Int_t)(1.+(yt-fStart)/fStep) - 1;
   if (dir) {
      Double_t divdiry = 1./TMath::Sqrt(1.+tyz*tyz);
      Double_t divdirz = -tyz*divdiry;
      Double_t dot = dir[1]*divdiry + dir[2]*divdirz;
      td.fNextIndex = ind;
      if (dot>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternParaY::MakeCopy(Bool_t reflect)
{
   TGeoPatternParaY *finder = new TGeoPatternParaY(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternParaY::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 2;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternParaY::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternParaY::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   Double_t dy = fStart+idiv*fStep+0.5*fStep;
   matrix.SetDx(fTxy*dy);
   matrix.SetDy(dy);
}

//______________________________________________________________________________
// TGeoPatternParaZ - a Z axis divison pattern for PARA shape
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternParaZ::TGeoPatternParaZ()
{
   fTxz = 0;
   fTyz = 0;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaZ::TGeoPatternParaZ(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   fTxz = ((TGeoPara*)vol->GetShape())->GetTxz();
   fTyz = ((TGeoPara*)vol->GetShape())->GetTyz();
   Double_t dz = ((TGeoPara*)vol->GetShape())->GetZ();
   fStart = -dz;
   fEnd = dz;
   fStep = 2*dz/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaZ::TGeoPatternParaZ(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   fTxz = ((TGeoPara*)vol->GetShape())->GetTxz();
   fTyz = ((TGeoPara*)vol->GetShape())->GetTyz();
   Double_t dz = ((TGeoPara*)vol->GetShape())->GetZ();
   fStart = -dz;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternParaZ::TGeoPatternParaZ(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   fTxz = ((TGeoPara*)vol->GetShape())->GetTxz();
   fTyz = ((TGeoPara*)vol->GetShape())->GetTyz();
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternParaZ::TGeoPatternParaZ(const TGeoPatternParaZ& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternParaZ& TGeoPatternParaZ::operator=(const TGeoPatternParaZ& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternParaZ::~TGeoPatternParaZ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternParaZ::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent = idiv;
   Double_t dz = fStart+idiv*fStep+0.5*fStep;
   td.fMatrix->SetDx(fTxz*dz);
   td.fMatrix->SetDy(fTyz*dz);
   td.fMatrix->SetDz((IsReflected())?-dz:dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternParaZ::IsOnBoundary(const Double_t *point) const
{
   Double_t seg = (point[2]-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// get the node division containing the query point

TGeoNode *TGeoPatternParaZ::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Double_t zt = point[2];
   Int_t ind = (Int_t)(1.+(zt-fStart)/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      if (dir[2]>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternParaZ::MakeCopy(Bool_t reflect)
{
   TGeoPatternParaZ *finder = new TGeoPatternParaZ(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternParaZ::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 3;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternParaZ::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternParaZ::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   Double_t dz = fStart+idiv*fStep+0.5*fStep;
   matrix.SetDx(fTxz*dz);
   matrix.SetDy(fTyz*dz);
   matrix.SetDz((IsReflected())?-dz:dz);
}

//______________________________________________________________________________
// TGeoPatternTrapZ - a Z axis divison pattern for TRAP or GTRA shapes
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternTrapZ::TGeoPatternTrapZ()
{
   fTxz = 0;
   fTyz = 0;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternTrapZ::TGeoPatternTrapZ(TGeoVolume *vol, Int_t ndivisions)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t theta = ((TGeoTrap*)vol->GetShape())->GetTheta();
   Double_t phi   = ((TGeoTrap*)vol->GetShape())->GetPhi();
   fTxz = TMath::Tan(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad());
   fTyz = TMath::Tan(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad());
   Double_t dz = ((TGeoArb8*)vol->GetShape())->GetDz();
   fStart = -dz;
   fEnd = dz;
   fStep = 2*dz/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternTrapZ::TGeoPatternTrapZ(TGeoVolume *vol, Int_t ndivisions, Double_t step)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t theta = ((TGeoTrap*)vol->GetShape())->GetTheta();
   Double_t phi   = ((TGeoTrap*)vol->GetShape())->GetPhi();
   fTxz = TMath::Tan(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad());
   fTyz = TMath::Tan(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad());
   Double_t dz = ((TGeoArb8*)vol->GetShape())->GetDz();
   fStart = -dz;
   fEnd = fStart + ndivisions*step;
   fStep       = step;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternTrapZ::TGeoPatternTrapZ(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
             :TGeoPatternFinder(vol, ndivisions)
{
   Double_t theta = ((TGeoTrap*)vol->GetShape())->GetTheta();
   Double_t phi   = ((TGeoTrap*)vol->GetShape())->GetPhi();
   fTxz = TMath::Tan(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad());
   fTyz = TMath::Tan(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad());
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternTrapZ::TGeoPatternTrapZ(const TGeoPatternTrapZ& pf) :
  TGeoPatternFinder(pf),
  fTxz(pf.fTxz),
  fTyz(pf.fTyz)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternTrapZ& TGeoPatternTrapZ::operator=(const TGeoPatternTrapZ& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      fTxz = pf.fTxz;
      fTyz = pf.fTyz;
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternTrapZ::~TGeoPatternTrapZ()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternTrapZ::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent = idiv;
   Double_t dz = fStart+idiv*fStep+0.5*fStep;
   td.fMatrix->SetDx(fTxz*dz);
   td.fMatrix->SetDy(fTyz*dz);
   td.fMatrix->SetDz((IsReflected())?-dz:dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternTrapZ::IsOnBoundary(const Double_t *point) const
{
   Double_t seg = (point[2]-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// get the node division containing the query point

TGeoNode *TGeoPatternTrapZ::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Double_t zt = point[2];
   Int_t ind = (Int_t)(1. + (zt-fStart)/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      if (dir[2]>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternTrapZ::MakeCopy(Bool_t reflect)
{
   TGeoPatternTrapZ *finder = new TGeoPatternTrapZ(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternTrapZ::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 3;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternTrapZ::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoMatrix *matrix = new TGeoTranslation(0.,0.,0.);
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoCombiTrans *combi = new TGeoCombiTrans();
   combi->RegisterYourself();
   combi->ReflectZ(kTRUE);
   combi->ReflectZ(kFALSE);
   return combi;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternTrapZ::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   Double_t dz = fStart+idiv*fStep+0.5*fStep;
   matrix.SetDx(fTxz*dz);
   matrix.SetDy(fTyz*dz);
   matrix.SetDz((IsReflected())?-dz:dz);
}

//______________________________________________________________________________
// TGeoPatternCylR - a cylindrical R divison pattern
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternCylR::TGeoPatternCylR()
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternCylR::TGeoPatternCylR(TGeoVolume *vol, Int_t ndivisions)
                :TGeoPatternFinder(vol, ndivisions)
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternCylR::TGeoPatternCylR(TGeoVolume *vol, Int_t ndivisions, Double_t step)
                :TGeoPatternFinder(vol, ndivisions)
{
   fStep       = step;
   CreateThreadData(1);
// compute start, end
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternCylR::TGeoPatternCylR(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
                :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternCylR::TGeoPatternCylR(const TGeoPatternCylR& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternCylR& TGeoPatternCylR::operator=(const TGeoPatternCylR& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternCylR::~TGeoPatternCylR()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternCylR::IsOnBoundary(const Double_t *point) const
{
   Double_t r = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
   Double_t seg = (r-fStart)/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternCylR::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
}

////////////////////////////////////////////////////////////////////////////////
/// find the node containing the query point

TGeoNode *TGeoPatternCylR::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   if (!td.fMatrix) td.fMatrix = gGeoIdentity;
   TGeoNode *node = 0;
   Double_t r = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
   Int_t ind = (Int_t)(1. + (r-fStart)/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      Double_t dot = point[0]*dir[0] + point[1]*dir[1];
      if (dot>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternCylR::MakeCopy(Bool_t reflect)
{
   TGeoPatternCylR *finder = new TGeoPatternCylR(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternCylR::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 1;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternCylR::CreateMatrix() const
{
   return gGeoIdentity;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternCylR::UpdateMatrix(Int_t, TGeoHMatrix &matrix) const
{
   matrix.Clear();
}

//______________________________________________________________________________
// TGeoPatternCylPhi - a cylindrical phi divison pattern
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternCylPhi::TGeoPatternCylPhi()
{
   fSinCos = 0;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor
/// compute step, start, end

TGeoPatternCylPhi::TGeoPatternCylPhi(TGeoVolume *vol, Int_t ndivisions)
                  :TGeoPatternFinder(vol, ndivisions)
{
   fStart = 0;
   fEnd = 0;
   fStep = 0;
   fSinCos     = new Double_t[2*fNdivisions];
   for (Int_t i = 0; i<fNdivisions; i++) {
      fSinCos[2*i] = TMath::Sin(TMath::DegToRad()*(fStart+0.5*fStep+i*fStep));
      fSinCos[2*i+1] = TMath::Cos(TMath::DegToRad()*(fStart+0.5*fStep+i*fStep));
   }
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternCylPhi::TGeoPatternCylPhi(TGeoVolume *vol, Int_t ndivisions, Double_t step)
                  :TGeoPatternFinder(vol, ndivisions)
{
   fStep       = step;
   fSinCos     = new Double_t[2*ndivisions];
   for (Int_t i = 0; i<fNdivisions; i++) {
      fSinCos[2*i] = TMath::Sin(TMath::DegToRad()*(fStart+0.5*fStep+i*fStep));
      fSinCos[2*i+1] = TMath::Cos(TMath::DegToRad()*(fStart+0.5*fStep+i*fStep));
   }
   CreateThreadData(1);
// compute start, end
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternCylPhi::TGeoPatternCylPhi(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
                  :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   if (fStart<0) fStart+=360;
   fEnd        = end;
   if (fEnd<0) fEnd+=360;
   if ((end-start)<0)
      fStep       = (end-start+360)/ndivisions;
   else
      fStep       = (end-start)/ndivisions;
   fSinCos     = new Double_t[2*ndivisions];
   for (Int_t idiv = 0; idiv<ndivisions; idiv++) {
      fSinCos[2*idiv] = TMath::Sin(TMath::DegToRad()*(start+0.5*fStep+idiv*fStep));
      fSinCos[2*idiv+1] = TMath::Cos(TMath::DegToRad()*(start+0.5*fStep+idiv*fStep));
   }
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternCylPhi::~TGeoPatternCylPhi()
{
   if (fSinCos) delete [] fSinCos;
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternCylPhi::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent = idiv;
   ((TGeoRotation*)td.fMatrix)->FastRotZ(&fSinCos[2*idiv]);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternCylPhi::IsOnBoundary(const Double_t *point) const
{
   Double_t phi = TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
   if (phi<0) phi += 360;
   Double_t ddp = phi - fStart;
   if (ddp<0) ddp+=360;
   Double_t seg = ddp/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// find the node containing the query point

TGeoNode *TGeoPatternCylPhi::FindNode(Double_t *point, const Double_t *dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Double_t phi = TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
   if (phi<0) phi += 360;
//   Double_t dphi = fStep*fNdivisions;
   Double_t ddp = phi - fStart;
   if (ddp<0) ddp+=360;
//   if (ddp>360) ddp-=360;
   Int_t ind = (Int_t)(1. + ddp/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      Double_t dot = point[0]*dir[1]-point[1]*dir[0];
      if (dot>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternCylPhi::MakeCopy(Bool_t reflect)
{
   TGeoPatternCylPhi *finder = new TGeoPatternCylPhi(*this);
   if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternCylPhi::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 2;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGeoVolume.

void TGeoPatternCylPhi::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TGeoPatternCylPhi::Class(), this);
      if (fNdivisions) {
         fSinCos     = new Double_t[2*fNdivisions];
         for (Int_t idiv = 0; idiv<fNdivisions; idiv++) {
            fSinCos[2*idiv] = TMath::Sin(TMath::DegToRad()*(fStart+0.5*fStep+idiv*fStep));
            fSinCos[2*idiv+1] = TMath::Cos(TMath::DegToRad()*(fStart+0.5*fStep+idiv*fStep));
         }
      }
   } else {
      R__b.WriteClassBuffer(TGeoPatternCylPhi::Class(), this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternCylPhi::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoRotation *matrix = new TGeoRotation();
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoRotation *rot = new TGeoRotation();
   rot->RegisterYourself();
   rot->ReflectZ(kTRUE);
   rot->ReflectZ(kFALSE);
   return rot;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternCylPhi::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   matrix.Clear();
   matrix.FastRotZ(&fSinCos[2*idiv]);
}

//______________________________________________________________________________
// TGeoPatternSphR - a spherical R divison pattern
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternSphR::TGeoPatternSphR()
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor
/// compute step, start, end

TGeoPatternSphR::TGeoPatternSphR(TGeoVolume *vol, Int_t ndivisions)
                :TGeoPatternFinder(vol, ndivisions)
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternSphR::TGeoPatternSphR(TGeoVolume *vol, Int_t ndivisions, Double_t step)
                :TGeoPatternFinder(vol, ndivisions)
{
   fStep       = step;
   CreateThreadData(1);
// compute start, end
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternSphR::TGeoPatternSphR(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
                :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternSphR::TGeoPatternSphR(const TGeoPatternSphR& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternSphR& TGeoPatternSphR::operator=(const TGeoPatternSphR& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternSphR::~TGeoPatternSphR()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternSphR::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent = idiv;
}
////////////////////////////////////////////////////////////////////////////////
/// find the node containing the query point

TGeoNode *TGeoPatternSphR::FindNode(Double_t * /*point*/, const Double_t * /*dir*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternSphR::MakeCopy(Bool_t)
{
   TGeoPatternSphR *finder = new TGeoPatternSphR(*this);
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternSphR::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 1;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternSphR::CreateMatrix() const
{
   return gGeoIdentity;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternSphR::UpdateMatrix(Int_t, TGeoHMatrix &matrix) const
{
   matrix.Clear();
}

//______________________________________________________________________________
// TGeoPatternSphTheta - a spherical theta divison pattern
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternSphTheta::TGeoPatternSphTheta()
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor
/// compute step, start, end

TGeoPatternSphTheta::TGeoPatternSphTheta(TGeoVolume *vol, Int_t ndivisions)
                    :TGeoPatternFinder(vol, ndivisions)
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternSphTheta::TGeoPatternSphTheta(TGeoVolume *vol, Int_t ndivisions, Double_t step)
                    :TGeoPatternFinder(vol, ndivisions)
{
   fStep       = step;
   CreateThreadData(1);
// compute start, end
}
////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoPatternSphTheta::TGeoPatternSphTheta(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
                    :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   fEnd        = end;
   fStep       = (end - start)/ndivisions;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternSphTheta::TGeoPatternSphTheta(const TGeoPatternSphTheta& pf) :
  TGeoPatternFinder(pf)
{
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternSphTheta& TGeoPatternSphTheta::operator=(const TGeoPatternSphTheta& pf)
{
   if(this!=&pf) {
      TGeoPatternFinder::operator=(pf);
      CreateThreadData(1);
   }
   return *this;
}
////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternSphTheta::~TGeoPatternSphTheta()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternSphTheta::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
}
////////////////////////////////////////////////////////////////////////////////
/// find the node containing the query point

TGeoNode *TGeoPatternSphTheta::FindNode(Double_t * /*point*/, const Double_t * /*dir*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternSphTheta::MakeCopy(Bool_t)
{
   TGeoPatternSphTheta *finder = new TGeoPatternSphTheta(*this);
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternSphTheta::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 2;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternSphTheta::CreateMatrix() const
{
   return gGeoIdentity;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternSphTheta::UpdateMatrix(Int_t, TGeoHMatrix &matrix) const
{
   matrix.Clear();
}

//______________________________________________________________________________
// TGeoPatternSphPhi - a spherical phi divison pattern
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternSphPhi::TGeoPatternSphPhi()
{
   fSinCos = 0;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor
/// compute step, start, end

TGeoPatternSphPhi::TGeoPatternSphPhi(TGeoVolume *vol, Int_t ndivisions)
                  :TGeoPatternFinder(vol, ndivisions)
{
   fStart = 0;
   fEnd = 360.;
   fStep = 360./ndivisions;
   CreateSinCos();
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor
/// compute start, end

TGeoPatternSphPhi::TGeoPatternSphPhi(TGeoVolume *vol, Int_t ndivisions, Double_t step)
                  :TGeoPatternFinder(vol, ndivisions)
{
   fStep       = step;
   CreateSinCos();
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// constructor
/// compute step

TGeoPatternSphPhi::TGeoPatternSphPhi(TGeoVolume *vol, Int_t ndivisions, Double_t start, Double_t end)
                  :TGeoPatternFinder(vol, ndivisions)
{
   fStart      = start;
   if (fStart<0) fStart+=360;
   fEnd        = end;
   if (fEnd<0) fEnd+=360;
   if ((end-start)<0)
      fStep       = (end-start+360)/ndivisions;
   else
      fStep       = (end-start)/ndivisions;
   CreateSinCos();
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPatternSphPhi::~TGeoPatternSphPhi()
{
   delete [] fSinCos;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the sincos table if it does not exist

Double_t *TGeoPatternSphPhi::CreateSinCos()
{
   fSinCos     = new Double_t[2*fNdivisions];
   for (Int_t idiv = 0; idiv<fNdivisions; idiv++) {
      fSinCos[2*idiv] = TMath::Sin(TMath::DegToRad()*(fStart+0.5*fStep+idiv*fStep));
      fSinCos[2*idiv+1] = TMath::Cos(TMath::DegToRad()*(fStart+0.5*fStep+idiv*fStep));
   }
   return fSinCos;
}

////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternSphPhi::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent = idiv;
   if (!fSinCos) CreateSinCos();
   ((TGeoRotation*)td.fMatrix)->FastRotZ(&fSinCos[2*idiv]);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the current point is on division boundary

Bool_t TGeoPatternSphPhi::IsOnBoundary(const Double_t *point) const
{
   Double_t phi = TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
   if (phi<0) phi += 360;
   Double_t ddp = phi - fStart;
   if (ddp<0) ddp+=360;
   Double_t seg = ddp/fStep;
   Double_t diff = seg - Int_t(seg);
   if (diff>0.5) diff = 1.-diff;
   if (diff<1e-8) return kTRUE;
   return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// find the node containing the query point

TGeoNode *TGeoPatternSphPhi::FindNode(Double_t * point, const Double_t * dir)
{
   ThreadData_t& td = GetThreadData();
   TGeoNode *node = 0;
   Double_t phi = TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
   if (phi<0) phi += 360;
//   Double_t dphi = fStep*fNdivisions;
   Double_t ddp = phi - fStart;
   if (ddp<0) ddp+=360;
//   if (ddp>360) ddp-=360;
   Int_t ind = (Int_t)(1. + ddp/fStep) - 1;
   if (dir) {
      td.fNextIndex = ind;
      Double_t dot = point[0]*dir[1]-point[1]*dir[0];
      if (dot>0) td.fNextIndex++;
      else td.fNextIndex--;
      if ((td.fNextIndex<0) || (td.fNextIndex>=fNdivisions)) td.fNextIndex = -1;
   }
   if ((ind<0) || (ind>=fNdivisions)) return node;
   node = GetNodeOffset(ind);
   cd(ind);
   return node;
}
////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this finder. Reflect by Z if required.

TGeoPatternFinder *TGeoPatternSphPhi::MakeCopy(Bool_t reflect)
{
   TGeoPatternSphPhi *finder = new TGeoPatternSphPhi(fVolume, fNdivisions, fStart, fEnd);
      if (!reflect) return finder;
   finder->Reflect();
   return finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPatternSphPhi::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   Int_t iaxis = 2;
   out << iaxis << ", " << fNdivisions << ", " << fStart << ", " << fStep;
}
////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternSphPhi::CreateMatrix() const
{
   if (!IsReflected()) {
      TGeoRotation *matrix = new TGeoRotation();
      matrix->RegisterYourself();
      return matrix;
   }
   TGeoRotation *rot = new TGeoRotation();
   rot->RegisterYourself();
   rot->ReflectZ(kTRUE);
   rot->ReflectZ(kFALSE);
   return rot;
}
////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternSphPhi::UpdateMatrix(Int_t idiv, TGeoHMatrix &matrix) const
{
   if (!fSinCos) ((TGeoPatternSphPhi*)this)->CreateSinCos();
   matrix.Clear();
   matrix.FastRotZ(&fSinCos[2*idiv]);
}

//______________________________________________________________________________
// TGeoPatternHoneycomb - a divison pattern specialized for honeycombs
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternHoneycomb::TGeoPatternHoneycomb()
{
   fNrows       = 0;
   fAxisOnRows  = 0;
   fNdivisions  = 0;
   fStart       = 0;
   CreateThreadData(1);
}
////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPatternHoneycomb::TGeoPatternHoneycomb(TGeoVolume *vol, Int_t nrows)
                     :TGeoPatternFinder(vol, nrows)
{
   fNrows = nrows;
   fAxisOnRows  = 0;
   fNdivisions  = 0;
   fStart       = 0;
   CreateThreadData(1);
// compute everything else
}
////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoPatternHoneycomb::TGeoPatternHoneycomb(const TGeoPatternHoneycomb& pfh) :
  TGeoPatternFinder(pfh),
  fNrows(pfh.fNrows),
  fAxisOnRows(pfh.fAxisOnRows),
  fNdivisions(pfh.fNdivisions),
  fStart(pfh.fStart)
{
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoPatternHoneycomb& TGeoPatternHoneycomb::operator=(const TGeoPatternHoneycomb& pfh)
{
   if(this!=&pfh) {
      TGeoPatternFinder::operator=(pfh);
      fNrows=pfh.fNrows;
      fAxisOnRows=pfh.fAxisOnRows;
      fNdivisions=pfh.fNdivisions;
      fStart=pfh.fStart;
      CreateThreadData(1);
   }
   return *this;
}
////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoPatternHoneycomb::~TGeoPatternHoneycomb()
{
}
////////////////////////////////////////////////////////////////////////////////
/// Update current division index and global matrix to point to a given slice.

void TGeoPatternHoneycomb::cd(Int_t idiv)
{
   ThreadData_t& td = GetThreadData();
   td.fCurrent=idiv;
}
////////////////////////////////////////////////////////////////////////////////
/// find the node containing the query point

TGeoNode *TGeoPatternHoneycomb::FindNode(Double_t * /*point*/, const Double_t * /*dir*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return new matrix of type used by  this finder.

TGeoMatrix* TGeoPatternHoneycomb::CreateMatrix() const
{
   return gGeoIdentity;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills external matrix with the local one corresponding to the given division
/// index.

void TGeoPatternHoneycomb::UpdateMatrix(Int_t, TGeoHMatrix &matrix) const
{
   matrix.Clear();
}
