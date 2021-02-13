// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveElement_hxx
#define ROOT7_REveElement_hxx

#include <ROOT/REveTypes.hxx>
#include <ROOT/REveVector.hxx>
#include <ROOT/REveProjectionBases.hxx>

#include <memory>

#include <nlohmann/json.hpp>

class TGeoMatrix;

namespace ROOT {
namespace Experimental {

class REveAunt;
class REveScene;
class REveCompound;
class REveTrans;
class REveRenderData;

//==============================================================================
// REveElement
// Base class for ROOT Event Visualization Environment (EVE)
// providing hierarchy management and selection and rendering control.
//==============================================================================

class REveElement
{
   friend class REveManager;
   friend class REveScene;

   REveElement& operator=(const REveElement&) = delete;

public:
   typedef std::list<REveElement*>              List_t;

   typedef std::set<REveElement*>               Set_t;

   typedef std::list<REveAunt*>                 AuntList_t;

private:
   ElementId_t      fElementId{0};        // Unique ID of an element.

protected:
   REveElement     *fMother {nullptr};
   REveScene       *fScene  {nullptr};
   REveElement     *fSelectionMaster {nullptr};

   ElementId_t get_mother_id() const;
   ElementId_t get_scene_id()  const;

   void assign_element_id_recurisvely();
   void assign_scene_recursively(REveScene* s);

protected:
   std::string      fName;                 //  Element name
   std::string      fTitle;                //  Element title / tooltip
   AuntList_t       fAunts;                //  List of aunts.
   List_t           fChildren;             //  List of children.
   TClass          *fChildClass {nullptr}; //  Class of acceptable children, others are rejected.
   REveCompound    *fCompound   {nullptr}; //  Compound this object belongs to.
   REveElement     *fVizModel   {nullptr}; //! Element used as model from VizDB.
   TString          fVizTag;               //  Tag used to query VizDB for model element.

   Int_t            fDenyDestroy{0};          //! Deny-destroy count.
   Bool_t           fDestroyOnZeroRefCnt{kTRUE};  //  Auto-destruct when ref-count reaches zero.

   Bool_t           fRnrSelf{kTRUE};                 //  Render this element.
   Bool_t           fRnrChildren{kTRUE};             //  Render children of this element.
   Bool_t           fCanEditMainColor{kFALSE};        //  Allow editing of main color.
   Bool_t           fCanEditMainTransparency{kFALSE}; //  Allow editing of main transparency.
   Bool_t           fCanEditMainTrans{kFALSE};        //  Allow editing of main transformation.

   Char_t           fMainTransparency{0};      //  Main-transparency variable.
   Color_t          fDefaultColor{kPink};  //  Default color for sub-classes that enable it.
   Color_t         *fMainColorPtr{nullptr};//  Pointer to main-color variable.
   std::unique_ptr<REveTrans> fMainTrans;   //  Pointer to main transformation matrix.

   void            *fUserData{nullptr};     //! Externally assigned and controlled user data.

   std::unique_ptr<REveRenderData> fRenderData;//! Vertex / normal / triangle index information for rendering.

   virtual void PreDeleteElement();
   virtual void RemoveElementsInternal();
   virtual void AnnihilateRecursively();

   static const std::string& ToString(Bool_t b);

public:
   REveElement(const std::string &name = "", const std::string &title = "");
   REveElement(const REveElement& e);
   virtual ~REveElement();

   ElementId_t GetElementId() const { return fElementId; }

   virtual REveElement* CloneElement() const;
   virtual REveElement* CloneElementRecurse(Int_t level = 0) const;
   virtual void         CloneChildrenRecurse(REveElement *dest, Int_t level = 0) const;

   const std::string &GetName()   const { return fName;  }
   const char* GetCName()  const { return fName.c_str();  }
   const std::string &GetTitle()  const { return fTitle; }
   const char* GetCTitle() const { return fTitle.c_str();  }

   virtual std::string GetHighlightTooltip(const std::set<int>&) const;

   void SetName (const std::string &name);
   void SetTitle(const std::string &title);
   void SetNameTitle(const std::string &name, const std::string &title);
   virtual void NameTitleChanged();

   const TString& GetVizTag() const               { return fVizTag; }
   void           SetVizTag(const TString& tag)   { fVizTag = tag;  }

   REveElement   *GetVizModel() const             { return fVizModel; }
   void           SetVizModel(REveElement* model);
   Bool_t         SetVizModelByTag();

   Bool_t         ApplyVizTag(const TString& tag, const TString& fallback_tag="");

   virtual void   PropagateVizParamsToProjecteds();
   virtual void   PropagateVizParamsToChildren(REveElement* el = nullptr);
   virtual void   CopyVizParams(const REveElement* el);
   virtual void   CopyVizParamsFromDB();
   void           SaveVizParams (std::ostream &out, const TString &tag, const TString &var);
   virtual void   WriteVizParams(std::ostream &out, const TString &var);

   REveCompound*  GetCompound()                { return fCompound; }
   void           SetCompound(REveCompound* c) { fCompound = c;    }

   bool         HasScene()  { return fScene  != nullptr; }
   bool         HasMother() { return fMother != nullptr; }

   REveScene*   GetScene()  { return fScene;  }
   REveElement* GetMother() { return fMother; }

   virtual void AddAunt(REveAunt *au);
   virtual void RemoveAunt(REveAunt *au);
   virtual void CheckReferenceCount(const std::string &from = "<unknown>");

   AuntList_t       &RefAunts()       { return fAunts; }
   const AuntList_t &RefAunts() const { return fAunts; }
   Int_t             NumAunts() const { return fAunts.size(); }
   Bool_t            HasAunts() const { return !fAunts.empty(); }

   TClass* GetChildClass() const { return fChildClass; }
   void    SetChildClass(TClass* c) { fChildClass = c; }

   List_t       &RefChildren()       { return fChildren; }
   const List_t &RefChildren() const { return fChildren; }
   Int_t         NumChildren() const { return fChildren.size(); }
   Bool_t        HasChildren() const { return !fChildren.empty(); }

   Bool_t       HasChild(REveElement *el);
   REveElement *FindChild(const TString &name, const TClass *cls = nullptr);
   REveElement *FindChild(TPRegexp &regexp, const TClass *cls = nullptr);
   Int_t        FindChildren(List_t &matches, const TString&  name, const TClass *cls = nullptr);
   Int_t        FindChildren(List_t &matches, TPRegexp& regexp, const TClass* cls = nullptr);
   REveElement *FirstChild() const;
   REveElement *LastChild () const;

   void EnableListElements(Bool_t rnr_self = kTRUE, Bool_t rnr_children = kTRUE);    // *MENU*
   void DisableListElements(Bool_t rnr_self = kFALSE, Bool_t rnr_children = kFALSE); // *MENU*

   Bool_t GetDestroyOnZeroRefCnt() const;
   void   SetDestroyOnZeroRefCnt(Bool_t d);

   Int_t  GetDenyDestroy() const;
   void   IncDenyDestroy();
   void   DecDenyDestroy();

   // --------------------------------

   TClass *IsA() const;

   virtual void ExportToCINT(const char *var_name); // *MENU*

   virtual Bool_t AcceptElement(REveElement *el);

   virtual void AddElement(REveElement *el);
   virtual void RemoveElement(REveElement *el);
   virtual void RemoveElementLocal(REveElement *el);
   virtual void RemoveElements();
   virtual void RemoveElementsLocal();

   virtual void AnnihilateElements();
   virtual void Annihilate();

   virtual void ProjectChild(REveElement *el, Bool_t same_depth = kTRUE);
   virtual void ProjectAllChildren(Bool_t same_depth = kTRUE);

   virtual void Destroy();                      // *MENU*
   virtual void DestroyOrWarn();
   virtual void DestroyElements();              // *MENU*

   virtual Bool_t CanEditElement() const { return kTRUE;    }
   virtual Bool_t SingleRnrState() const { return kFALSE;   }
   virtual Bool_t GetRnrSelf()     const { return fRnrSelf; }
   virtual Bool_t GetRnrChildren() const { return fRnrChildren; }
   virtual Bool_t GetRnrState()    const { return fRnrSelf && fRnrChildren; }
   virtual Bool_t GetRnrAnything() const { return fRnrSelf || (fRnrChildren && HasChildren()); }
   virtual Bool_t SetRnrSelf(Bool_t rnr);
   virtual Bool_t SetRnrChildren(Bool_t rnr);
   virtual Bool_t SetRnrSelfChildren(Bool_t rnr_self, Bool_t rnr_children);
   virtual Bool_t SetRnrState(Bool_t rnr);
   virtual void   PropagateRnrStateToProjecteds();

   void           SetupDefaultColorAndTransparency(Color_t col, Bool_t can_edit_color, Bool_t can_edit_transparency);

   virtual Bool_t CanEditMainColor() const   { return fCanEditMainColor; }
   void           SetEditMainColor(Bool_t x) { fCanEditMainColor = x;    }
   Color_t       *GetMainColorPtr()        const   { return fMainColorPtr;     }
   void           SetMainColorPtr(Color_t *colptr) { fMainColorPtr = colptr;   }

   virtual Bool_t  HasMainColor() const { return fMainColorPtr != nullptr; }
   virtual Color_t GetMainColor() const { return fMainColorPtr ? *fMainColorPtr : 0; }
   virtual void    SetMainColor(Color_t color);
   void            SetMainColorPixel(Pixel_t pixel);
   void            SetMainColorRGB(UChar_t r, UChar_t g, UChar_t b);
   void            SetMainColorRGB(Float_t r, Float_t g, Float_t b);
   virtual void    PropagateMainColorToProjecteds(Color_t color, Color_t old_color);

   virtual Bool_t  CanEditMainTransparency() const   { return fCanEditMainTransparency; }
   void            SetEditMainTransparency(Bool_t x) { fCanEditMainTransparency = x; }
   virtual Char_t  GetMainTransparency()     const { return fMainTransparency; }
   virtual void    SetMainTransparency(Char_t t);
   void            SetMainAlpha(Float_t alpha);
   virtual void    PropagateMainTransparencyToProjecteds(Char_t t, Char_t old_t);

   virtual Bool_t     CanEditMainTrans() const { return fCanEditMainTrans; }
   virtual Bool_t     HasMainTrans()     const { return fMainTrans.get() != nullptr; }
   virtual REveTrans* PtrMainTrans(Bool_t create=kTRUE);
   virtual REveTrans& RefMainTrans();
   virtual void       InitMainTrans(Bool_t can_edit=kTRUE);
   virtual void       DestroyMainTrans();

   virtual void SetTransMatrix(Double_t *carr);
   virtual void SetTransMatrix(const TGeoMatrix &mat);

   virtual Int_t WriteCoreJson(nlohmann::json &cj, Int_t rnr_offset);
   virtual void  BuildRenderData();

   void* GetUserData() const   { return fUserData; }
   void  SetUserData(void* ud) { fUserData = ud;   }

   REveRenderData *GetRenderData() const { return fRenderData.get(); }


   // Selection state and management
   //--------------------------------

protected:

   enum ECompoundSelectionColorBits
   {
      kCSCBImplySelectAllChildren           = BIT(0), // compound will select all children
      kCSCBTakeMotherAsMaster               = BIT(1), // element will take its mother as master
      kCSCBApplyMainColorToAllChildren      = BIT(2), // compound will apply color change to all children
      kCSCBApplyMainColorToMatchingChildren = BIT(3), // compound will apply color change to all children with matching color
      kCSCBApplyMainTransparencyToAllChildren      = BIT(4), // compound will apply transparency change to all children
      kCSCBApplyMainTransparencyToMatchingChildren = BIT(5)  // compound will apply transparency change to all children with matching color
   };

   enum EDestruct
   {
      kNone,
      kStandard,
      kAnnihilate
   };

   Short_t fImpliedSelected{0};   // How many times the element is implied selected -- needed during destruction.
   Bool_t  fPickable{0};          // Can element be selected.
   UChar_t fCSCBits{0};           // Compound Selection Color flags.

public:
   Bool_t IsPickable()    const { return fPickable; }
   void   SetPickable(Bool_t p) { fPickable = p; }
   void   SetPickableRecursively(Bool_t p);

   virtual REveElement* GetSelectionMaster();
   void         SetSelectionMaster(REveElement *el) { fSelectionMaster = el; }

   virtual void FillImpliedSelectedSet(Set_t& impSelSet);

   void   IncImpliedSelected() { ++fImpliedSelected; }
   void   DecImpliedSelected() { --fImpliedSelected; }
   int    GetImpliedSelected() { return fImpliedSelected; }

   void   RecheckImpliedSelections();

   void   SetCSCBits(UChar_t f)   { fCSCBits |=  f; }
   void   ResetCSCBits(UChar_t f) { fCSCBits &= ~f; }
   Bool_t TestCSCBits(UChar_t f) const { return (fCSCBits & f) != 0; }

   void   ResetAllCSCBits()                     { fCSCBits  =  0; }
   void   CSCImplySelectAllChildren()           { fCSCBits |= kCSCBImplySelectAllChildren; }
   void   CSCTakeMotherAsMaster()               { fCSCBits |= kCSCBTakeMotherAsMaster;  }
   void   CSCApplyMainColorToAllChildren()      { fCSCBits |= kCSCBApplyMainColorToAllChildren; }
   void   CSCApplyMainColorToMatchingChildren() { fCSCBits |= kCSCBApplyMainColorToMatchingChildren; }
   void   CSCApplyMainTransparencyToAllChildren()      { fCSCBits |= kCSCBApplyMainTransparencyToAllChildren; }
   void   CSCApplyMainTransparencyToMatchingChildren() { fCSCBits |= kCSCBApplyMainTransparencyToMatchingChildren; }

   virtual bool RequiresExtraSelectionData() const { return false; }
   virtual void FillExtraSelectionData(nlohmann::json&, const std::set<int>&) const {}

   // Change-stamping and change bits
   //---------------------------------

   enum EChangeBits
   {
      kCBColorSelection =  BIT(0), // Main color or select/hilite state changed.
      kCBTransBBox      =  BIT(1), // Transformation matrix or bounding-box changed.
      kCBObjProps       =  BIT(2), // Object changed, requires dropping its display-lists.
      kCBVisibility     =  BIT(3),  // Rendering of self/children changed.
      kCBElementAdded   =  BIT(4) // Element was added to a new parent.
      // kCBElementRemoved = BIT()  // Element was removed from a parent.

      // Deletions are handled in a special way in REveManager::PreDeleteElement().
   };

protected:
   UChar_t      fChangeBits{0};  //!
   Char_t       fDestructing{kNone}; //!

public:
   void StampColorSelection() { AddStamp(kCBColorSelection); }
   void StampTransBBox()      { AddStamp(kCBTransBBox); }
   void StampObjProps()       { AddStamp(kCBObjProps); }
   void StampObjPropsPreChk() { if ( ! (fChangeBits & kCBObjProps)) AddStamp(kCBObjProps); }
   void StampVisibility()     { AddStamp(kCBVisibility); }
   void StampElementAdded()   { AddStamp(kCBElementAdded); }
   // void StampElementRemoved() { AddStamp(kCBElementRemoved); }
   virtual void AddStamp(UChar_t bits);
   virtual void ClearStamps() { fChangeBits = 0; }

   UChar_t GetChangeBits() const { return fChangeBits; }

   // Menu entries for VizDB communication (here so they are last in the menu).

   void VizDB_Apply(const std::string& tag);    // *MENU*
   void VizDB_Reapply();                        // *MENU*
   void VizDB_UpdateModel(Bool_t update=kTRUE); // *MENU*
   void VizDB_Insert(const std::string& tag, Bool_t replace=kTRUE, Bool_t update=kTRUE); // *MENU*
};


//==============================================================================
// REveAunt
//==============================================================================

class REveAunt
{
public:
   virtual ~REveAunt() {}

   virtual bool HasNiece(REveElement *el) const = 0;
   virtual bool HasNieces() const = 0;

   virtual bool AcceptNiece(REveElement *) { return true; }

   virtual void AddNiece(REveElement *el)
   {
      // XXXX Check AcceptNiece() -- throw if not !!!!
      el->AddAunt(this);
      AddNieceInternal(el);
   }
   virtual void AddNieceInternal(REveElement *el) = 0;

   virtual void RemoveNiece(REveElement *el)
   {
      RemoveNieceInternal(el);
      el->RemoveAunt(this);
   }
   virtual void RemoveNieceInternal(REveElement *el) = 0;

   virtual void RemoveNieces() = 0;
};


//==============================================================================
// REveAuntAsList
//==============================================================================

class REveAuntAsList : public REveAunt
{
protected:
   REveElement::List_t fNieces;

public:
   virtual ~REveAuntAsList()
   {
      for (auto &n : fNieces) n->RemoveAunt(this);
   }

   bool HasNiece(REveElement *el) const override
   {
      return std::find(fNieces.begin(), fNieces.end(), el) != fNieces.end();
   }

   bool HasNieces() const override
   {
      return ! fNieces.empty();
   }

   void AddNieceInternal(REveElement *el) override
   {
      fNieces.push_back(el);
   }

   void RemoveNieceInternal(REveElement *el) override
   {
      fNieces.remove(el);
   }

   void RemoveNieces() override
   {
      for (auto &n : fNieces) n->RemoveAunt(this);
      fNieces.clear();
   }
};

} // namespace Experimental
} // namespace ROOT

#endif
