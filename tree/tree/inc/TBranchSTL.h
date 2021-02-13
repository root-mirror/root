// @(#)root/tree:$Id$
// author: Lukasz Janyst <ljanyst@cern.ch>

//------------------------------------------------------------------------------
// file:   TBranchSTL.h
//------------------------------------------------------------------------------

#ifndef ROOT_TBranchSTL
#define ROOT_TBranchSTL

#include "TBranch.h"
#include "TIndArray.h"

#include <map>
#include <vector>

class TTree;
class TVirtualCollectionProxy;
class TStreamerInfo;
class TBranchElement;

class TBranchSTL: public TBranch {
   public:
      TBranchSTL();
      TBranchSTL( TTree* tree, const char* name,
                  TVirtualCollectionProxy* collProxy,
                  Int_t buffsize, Int_t splitlevel );
      TBranchSTL( TBranch* parent, const char* name,
                  TVirtualCollectionProxy* collProxy,
                  Int_t buffsize, Int_t splitlevel,
                  TStreamerInfo* info, Int_t id );
      virtual ~TBranchSTL();
      virtual void           Browse( TBrowser *b );
      virtual Bool_t         IsFolder() const;
      virtual const char    *GetClassName() const { return fClassName.Data(); }
      virtual Int_t          GetExpectedType(TClass *&clptr,EDataType &type);
      virtual Int_t          GetEntry( Long64_t entry = 0, Int_t getall = 0 );
      virtual TStreamerInfo *GetInfo() const;
      virtual void           Print(Option_t*) const;
      virtual void           SetAddress( void* addr );

      ClassDef( TBranchSTL, 1 ) //Branch handling STL collection of pointers

   private:

      void ReadLeavesImpl( TBuffer& b );
      void FillLeavesImpl( TBuffer& b );
      virtual Int_t          FillImpl(ROOT::Internal::TBranchIMTHelper *);

      struct ElementBranchHelper_t
      {
         ElementBranchHelper_t():
            fBranch( 0 ), fPointers( 0 ), fId( 0 ),
            fBaseOffset( 0 ), fPosition( 0 ) {}

         TBranchElement*     fBranch;
         std::vector<void*>* fPointers;
         UChar_t             fId;
         UInt_t              fBaseOffset;
         Int_t               fPosition;
      };

      typedef std::map<TClass*, ElementBranchHelper_t> BranchMap_t;
      BranchMap_t fBranchMap;                           ///<! Branch map
      std::vector<ElementBranchHelper_t> fBranchVector; ///<! Branch vector

      TVirtualCollectionProxy* fCollProxy;    ///<! Collection proxy
      TBranch*                 fParent;       ///<! Parent of this branch
      TClass*                  fIndArrayCl;   ///<! Class of the ind array
      TIndArray                fInd;          ///<! Indices
      TString                  fContName;     ///<  Class name of referenced object
      TString                  fClassName;    ///<  Name of the parent class, if we're the data member
      mutable Int_t            fClassVersion; ///<  Version number of the class
      UInt_t                   fClCheckSum;   ///<  Class checksum
      mutable TStreamerInfo   *fInfo;         ///<! The streamer info
      char*                    fObject;       ///<! Pointer to object at address or the
      Int_t                    fID;           ///<  Element serial number in the streamer info
};

#endif // ROOT_TBranchSTL
