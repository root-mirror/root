// @(#):$Id$
// Author: Andrei Gheata   30/05/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoBoolNode
#define ROOT_TGeoBoolNode

#include "TObject.h"

#include <mutex>
#include <vector>

// forward declarations
class TGeoShape;
class TGeoMatrix;
class TGeoHMatrix;

class TGeoBoolNode : public TObject
{
public:
enum EGeoBoolType {
   kGeoUnion,
   kGeoIntersection,
   kGeoSubtraction
};
   struct ThreadData_t
   {
      Int_t          fSelected;       // ! selected branch

      ThreadData_t();
      ~ThreadData_t();
   };
   ThreadData_t&     GetThreadData()   const;
   void              ClearThreadData() const;
   void              CreateThreadData(Int_t nthreads);
private:
   TGeoBoolNode(const TGeoBoolNode&) = delete;
   TGeoBoolNode& operator=(const TGeoBoolNode&) = delete;

protected:
   TGeoShape        *fLeft{nullptr};         // shape on the left branch
   TGeoShape        *fRight{nullptr};        // shape on the right branch
   TGeoMatrix       *fLeftMat{nullptr};      // transformation that applies to the left branch
   TGeoMatrix       *fRightMat{nullptr};     // transformation that applies to the right branch
   Int_t             fNpoints{0};            //! number of points on the mesh
   Double_t         *fPoints{nullptr};      //! array of mesh points

   mutable std::vector<ThreadData_t*> fThreadData;    //! Navigation data per thread
   mutable Int_t                      fThreadSize{0}; //! Size for the navigation data array
   mutable std::mutex                 fMutex;         //! Mutex for thread data access
// methods
   Bool_t            MakeBranch(const char *expr, Bool_t left);
   void              AssignPoints(Int_t npoints, Double_t *points);

public:
   // constructors
   TGeoBoolNode();
   TGeoBoolNode(const char *expr1, const char *expr2);
   TGeoBoolNode(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat = nullptr, TGeoMatrix *rmat = nullptr);

   // destructor
   virtual ~TGeoBoolNode();
   // methods
   virtual void      ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin) = 0;
   virtual void      ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) = 0;
   virtual Bool_t    Contains(const Double_t *point) const = 0;
   virtual Int_t     DistanceToPrimitive(Int_t px, Int_t py) = 0;
   virtual Double_t  DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                               Double_t step=0, Double_t *safe=nullptr) const = 0;
   virtual Double_t  DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                               Double_t step=0, Double_t *safe=nullptr) const = 0;
   virtual EGeoBoolType GetBooleanOperator() const = 0;
   virtual Int_t     GetNpoints() = 0;
   TGeoMatrix       *GetLeftMatrix() const {return fLeftMat;}
   TGeoMatrix       *GetRightMatrix() const {return fRightMat;}
   TGeoShape        *GetLeftShape() const {return fLeft;}
   TGeoShape        *GetRightShape() const {return fRight;}
   virtual TGeoBoolNode *MakeClone() const = 0;
           void      Paint(Option_t *option) override;
   void              RegisterMatrices();
   Bool_t            ReplaceMatrix(TGeoMatrix *mat, TGeoMatrix *newmat);
   virtual Double_t  Safety(const Double_t *point, Bool_t in=kTRUE) const = 0;
           void      SavePrimitive(std::ostream &out, Option_t *option = "")  override;
   virtual void      SetPoints(Double_t *points) const;
   virtual void      SetPoints(Float_t *points)  const;
   void              SetSelected(Int_t sel);
   virtual void      Sizeof3D() const;

   ClassDefOverride(TGeoBoolNode, 1)              // a boolean node
};

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TGeoUnion - Boolean node representing a union between two components.    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TGeoUnion : public TGeoBoolNode
{
public:
   // constructors
   TGeoUnion();
   TGeoUnion(const char *expr1, const char *expr2);
   TGeoUnion(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat = nullptr, TGeoMatrix *rmat = nullptr);

   // destructor
   virtual ~TGeoUnion();
   // methods
   void      ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin) override;
   void      ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) override;
   Bool_t    Contains(const Double_t *point) const override;
   Int_t     DistanceToPrimitive(Int_t px, Int_t py) override;
   Double_t  DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                            Double_t step = 0, Double_t *safe = nullptr) const override;
   Double_t  DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                             Double_t step = 0, Double_t *safe=nullptr) const  override;
   EGeoBoolType GetBooleanOperator() const  override {return kGeoUnion;}
   Int_t     GetNpoints() override;
   Double_t  Safety(const Double_t *point, Bool_t in=kTRUE) const  override;
   void      SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void      Sizeof3D() const override;

   //CS specific
   TGeoBoolNode *MakeClone() const override;
   void      Paint(Option_t *option) override;

   ClassDefOverride(TGeoUnion, 1)              // union node
};

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TGeoIntersection - Boolean node representing an intersection between two //
// components.                                                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TGeoIntersection : public TGeoBoolNode
{
public:
   // constructors
   TGeoIntersection();
   TGeoIntersection(const char *expr1, const char *expr2);
   TGeoIntersection(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat = nullptr, TGeoMatrix *rmat = nullptr);

   // destructor
   virtual ~TGeoIntersection();
   // methods
   void      ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin) override;
   void      ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) override;
   Bool_t    Contains(const Double_t *point) const override;
   Int_t     DistanceToPrimitive(Int_t px, Int_t py) override;
   Double_t  DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                            Double_t step = 0, Double_t *safe = nullptr) const override;
   Double_t  DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                             Double_t step = 0, Double_t *safe = nullptr) const override;
   EGeoBoolType GetBooleanOperator() const override { return kGeoIntersection; }
   Int_t     GetNpoints() override;
   Double_t  Safety(const Double_t *point, Bool_t in=kTRUE) const override;
   void      SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void      Sizeof3D() const override;

   //CS specific
   TGeoBoolNode *MakeClone() const override;
   void      Paint(Option_t *option) override;

   ClassDefOverride(TGeoIntersection, 1)              // intersection node
};

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TGeoSubtraction - Boolean node representing a subtraction.               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TGeoSubtraction : public TGeoBoolNode
{
public:
   // constructors
   TGeoSubtraction();
   TGeoSubtraction(const char *expr1, const char *expr2);
   TGeoSubtraction(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat = nullptr, TGeoMatrix *rmat = nullptr);

   // destructor
   virtual ~TGeoSubtraction();
   // methods
   void      ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin) override;
   void      ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) override;
   Bool_t    Contains(const Double_t *point) const override;
   Int_t     DistanceToPrimitive(Int_t px, Int_t py) override;
   Double_t  DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                            Double_t step = 0, Double_t *safe = nullptr) const override;
   Double_t  DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                             Double_t step = 0, Double_t *safe = nullptr) const override;
   EGeoBoolType GetBooleanOperator() const override { return kGeoSubtraction; }
   Int_t     GetNpoints() override;
   Double_t  Safety(const Double_t *point, Bool_t in=kTRUE) const override;
   void      SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void      Sizeof3D() const override;

   //CS specific
   TGeoBoolNode *MakeClone() const override;
   void      Paint(Option_t *option) override;

   ClassDefOverride(TGeoSubtraction, 1)              // subtraction node
};

#endif
