// @(#)root/mathcore:$Id$
// Authors: C. Gumpert    09/2011
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2011 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
//
// Header file for KDTree class
//


#ifndef ROOT_Math_KDTree
#define ROOT_Math_KDTree

//STL header
#include <cassert>
#include <vector>
#include <cmath>
#include <utility>

// ROOT include(s)
#include "RtypesCore.h"

namespace ROOT
{
  namespace Math
  {

     //______________________________________________________________________________
     //Begin_Html
     //End_Html
     template<class _DataPoint>
     class KDTree
     {
     public:

        typedef _DataPoint                        point_type;
        typedef typename _DataPoint::value_type   value_type;
        static UInt_t Dimension() {return _DataPoint::Dimension();}
        enum eSplitOption {
           kEffective = 0,                         //split according to effective entries
           kBinContent                             //split according to bin content
        };

     private:

        class ComparePoints
        {
        public:
           Bool_t   operator()(const point_type* pFirst,const point_type* pSecond) const;

           UInt_t   GetAxis() const       {return fAxis;}
           void     SetAxis(UInt_t iAxis) {fAxis = iAxis;}

        private:
           UInt_t   fAxis; //axis at which the points are compared
        };

        class Cut
        {
        public:
           Cut():fAxis(0),fCutValue(0) {}
           Cut(UInt_t iAxis,Double_t fNewCutValue):fAxis(iAxis),fCutValue(fNewCutValue) {}
           ~Cut() {}

           UInt_t       GetAxis() const                   {return fAxis;}
           value_type   GetCutValue() const               {return fCutValue;}
           void         SetAxis(UInt_t iAxis)             {fAxis = iAxis;}
           void         SetCutValue(Double_t fNewCutValue) {fCutValue = fNewCutValue;}

           Bool_t       operator<(const point_type& rPoint) const;
           Bool_t       operator>(const point_type& rPoint) const;

        private:
           UInt_t    fAxis;       //axis at which the splitting is done
           Double_t  fCutValue;   //split value
        };

        //forward declarations
        class BaseNode;
        class HeadNode;
        class SplitNode;
        class BinNode;
        class TerminalNode;

        class BaseNode
        {
        public:
           //constructor and destructor
           BaseNode(BaseNode* pParent = 0);
           virtual ~BaseNode();

           //providing usual functionality of a tree
           virtual BaseNode*        Clone() = 0;
           virtual const BinNode*   FindNode(const point_type& rPoint) const = 0;
           virtual void             GetClosestPoints(const point_type& rRef,UInt_t nPoints,std::vector<std::pair<const _DataPoint*,Double_t> >& vFoundPoints) const = 0;
           virtual void             GetPointsWithinDist(const point_type& rRef,value_type fDist,std::vector<const point_type*>& vFoundPoints) const = 0;
           virtual Bool_t           Insert(const point_type& rPoint) = 0;
           virtual void             Print(int iRow = 0) const = 0;

           //navigating the tree
           BaseNode*&                LeftChild()        {return fLeftChild;}
           const BaseNode*           LeftChild() const  {return fLeftChild;}
           BaseNode*&                Parent()           {return fParent;}
           const BaseNode*           Parent() const     {return fParent;}
           BaseNode*&                RightChild()       {return fRightChild;}
           const BaseNode*           RightChild() const {return fRightChild;}

           //information about relative position of current node
           BaseNode*&                GetParentPointer();
           virtual Bool_t            IsHeadNode() const {return false;}
           Bool_t                    IsLeftChild() const;

        private:
           // node should never be copied or assigned
           BaseNode(const BaseNode& ) {}
           BaseNode& operator=(const BaseNode& ) {return *this;}

           //links to adjacent nodes
           BaseNode*                 fParent;     //!pointer to parent node
           BaseNode*                 fLeftChild;  //!pointer to left child
           BaseNode*                 fRightChild; //!pointer to right child
        };

        class HeadNode : public BaseNode
        {
        public:
           //constructor and destructor
           HeadNode(BaseNode& rNode):BaseNode(&rNode) {}
           virtual ~HeadNode() {delete Parent();}

           //delegate everything to the actual root node of the tree
           virtual const BinNode*   FindNode(const point_type& rPoint) const {return Parent()->FindNode(rPoint);}
           virtual void             GetClosestPoints(const point_type& rRef,UInt_t nPoints,std::vector<std::pair<const _DataPoint*,Double_t> >& vFoundPoints) const;
           virtual void             GetPointsWithinDist(const point_type& rRef,value_type fDist,std::vector<const _DataPoint*>& vFoundPoints) const;
           virtual Bool_t           Insert(const point_type& rPoint) {return Parent()->Insert(rPoint);}
           virtual void             Print(Int_t) const {Parent()->Print();}

        private:
           // node should never be copied
           HeadNode(const HeadNode& ) {}
           HeadNode& operator=(const HeadNode& ) {return *this;}

           virtual HeadNode*        Clone();
           virtual bool             IsHeadNode() const {return true;}

           // only delegate everything else is private and should not be used
           using BaseNode::Parent;
           using BaseNode::LeftChild;
           using BaseNode::RightChild;

           using BaseNode::GetParentPointer;
           using BaseNode::IsLeftChild;
        };

        class SplitNode : public BaseNode
        {
        public:
           // constructors and destructors
           SplitNode(UInt_t iAxis,Double_t fCutValue,BaseNode& rLeft,BaseNode& rRight,BaseNode* pParent = 0);
           virtual ~SplitNode();

           //accessing information about this split node
           const Cut*               GetCut() const {return fCut;}
           virtual void             Print(Int_t iRow = 0) const;

        private:
           // node should never be copied
           SplitNode(const SplitNode& ) {}
           SplitNode& operator=(const SplitNode& ) {return *this;}

           virtual SplitNode*       Clone();
           virtual const BinNode*   FindNode(const point_type& rPoint) const;
           virtual void             GetClosestPoints(const point_type& rRef,UInt_t nPoints,std::vector<std::pair<const _DataPoint*,Double_t> >& vFoundPoints) const;
           virtual void             GetPointsWithinDist(const point_type& rRef,value_type fDist,std::vector<const _DataPoint*>& vFoundPoints) const;
           virtual Bool_t           Insert(const point_type& rPoint);

           const Cut*               fCut;     //pointer to cut object owned by this node
        };

        class BinNode : public BaseNode
        {
        protected:
           //save some typing
           typedef std::pair<value_type,value_type> tBoundary;
        public:
           // constructors and destructors
           BinNode(BaseNode* pParent = 0);
           BinNode(const BinNode& copy);
           virtual ~BinNode() {}

           // usual bin operations
           virtual void                            EmptyBin();
           virtual const BinNode*                  FindNode(const point_type& rPoint) const;
           point_type                              GetBinCenter() const;
           Double_t                                GetBinContent() const {return GetSumw();}
#ifndef _AIX
           virtual const std::vector<tBoundary>&   GetBoundaries() const {return fBoundaries;}
#else
           virtual void GetBoundaries() const { }
#endif
           Double_t                                GetDensity() const {return GetBinContent()/GetVolume();}
           Double_t                                GetEffectiveEntries() const {return (GetSumw2()) ? std::pow(GetSumw(),2)/GetSumw2() : 0;}
           UInt_t                                  GetEntries() const {return fEntries;}
           Double_t                                GetVolume() const;
           Double_t                                GetSumw() const {return fSumw;}
           Double_t                                GetSumw2() const {return fSumw2;}
           virtual Bool_t                          Insert(const point_type& rPoint);
           Bool_t                                  IsInBin(const point_type& rPoint) const;
           virtual void                            Print(int iRow = 0) const;

        protected:
           virtual BinNode*                        Clone();

           // intrinsic bin properties
           std::vector<tBoundary>                  fBoundaries;    //bin boundaries
           Double_t                                fSumw;          //sum of weights
           Double_t                                fSumw2;         //sum of weights^2
           UInt_t                                  fEntries;       //number of entries

        private:
           BinNode& operator=(const BinNode& rhs);

           // bin does not contain any point like information
           virtual void                    GetClosestPoints(const point_type&,UInt_t,std::vector<std::pair<const _DataPoint*,Double_t> >&) const {}
           virtual void                    GetPointsWithinDist(const point_type&,value_type,std::vector<const point_type*>&) const {}

           // a bin does not have children
           using BaseNode::LeftChild;
           using BaseNode::RightChild;
        };

        class TerminalNode : public BinNode
        {
           friend class KDTree<_DataPoint>;
           //save some typing
           typedef std::pair<value_type,value_type> tBoundary;

        public:
           //constructor and desctructor
           TerminalNode(Double_t iBucketSize,BaseNode* pParent = 0);
           virtual ~TerminalNode();

           virtual void                            EmptyBin();
#ifndef _AIX
           virtual const std::vector<tBoundary>&   GetBoundaries() const;
#else
           virtual void GetBoundaries() const;
#endif
           virtual void                            GetClosestPoints(const point_type& rRef,UInt_t nPoints,std::vector<std::pair<const _DataPoint*,Double_t> >& vFoundPoints) const;
           const std::vector<const point_type*>&   GetPoints() const {return fDataPoints;}
           virtual void                            GetPointsWithinDist(const point_type& rRef,value_type fDist,std::vector<const _DataPoint*>& vFoundPoints) const;
           virtual void                            Print(int iRow = 0) const;

        private:
           // node should never be copied
           TerminalNode(const TerminalNode& ) {}
           TerminalNode& operator=(const TerminalNode& ) {return *this;}

           // save some typing
           typedef typename std::vector<const point_type* >::iterator         data_it;
           typedef typename std::vector<const point_type* >::const_iterator   const_data_it;

           // creating new Terminal Node when splitting, copying elements in the given range
           TerminalNode(Double_t iBucketSize,UInt_t iSplitAxis,data_it first,data_it end);

           //tree operations
           virtual BinNode*                        Clone() {return ConvertToBinNode();}
           BinNode*                                ConvertToBinNode();
           virtual const BinNode*                  FindNode(const point_type&) const {return this;}
           virtual Bool_t                          Insert(const point_type& rPoint);
           void                                    Split();
           void                                    SetOwner(Bool_t bIsOwner = true) {fOwnData = bIsOwner;}
           void                                    SetSplitOption(eSplitOption opt) {fSplitOption = opt;}
           data_it                                 SplitEffectiveEntries();
           data_it                                 SplitBinContent();
           void                                    UpdateBoundaries();

           Bool_t                                  fOwnData;       // terminal node owns the data objects (default = false)
           eSplitOption                            fSplitOption;   // according to which figure of merit the node is splitted
           Double_t                                fBucketSize;    // target number of entries per bucket
           UInt_t                                  fSplitAxis;     // axis at which the next split will occur
           std::vector<const _DataPoint*>          fDataPoints;    // data points in this bucket
        };

     public:
        //////////////////////////////////////////////////////////////////////
        //
        // template<class _DataPoint> class KDTree<_DataPoint>::iterator
        //
        //////////////////////////////////////////////////////////////////////
        typedef BinNode Bin;
        class iterator
        {
           friend class KDTree<_DataPoint>;
        public:
           iterator(): fBin(0) {}
           iterator(const iterator& copy): fBin(copy.fBin) {}
           ~iterator() {}

           iterator&         operator++();
           const iterator&   operator++() const;
           iterator          operator++(int);
           const iterator    operator++(int) const;
           iterator&         operator--();
           const iterator&   operator--() const;
           iterator          operator--(int);
           const iterator    operator--(int) const;
           bool              operator==(const iterator& rIterator) const {return (fBin == rIterator.fBin);}
           bool              operator!=(const iterator& rIterator) const {return !(*this == rIterator);}
           iterator&         operator=(const iterator& rhs);
           Bin&              operator*() {return *fBin;}
           const Bin&        operator*() const {return *fBin;}
           Bin*              operator->() {return fBin;}
           const Bin*        operator->() const {return fBin;}

           TerminalNode*     TN() {assert(dynamic_cast<TerminalNode*>(fBin)); return (TerminalNode*)fBin;}

        private:
           iterator(BinNode* pNode): fBin(pNode) {}

           Bin*     Next() const;
           Bin*     Previous() const;

           mutable Bin* fBin;
        };

        //constructor and destructor
        KDTree(UInt_t iBucketSize);
        ~KDTree();

        //public member functions
        void            EmptyBins();
        iterator        End();
        const iterator  End() const;
        const Bin*      FindBin(const point_type& rPoint) const {return fHead->FindNode(rPoint);}
        iterator        First();
        const iterator  First() const;
        void            Freeze();
        Double_t        GetBucketSize() const {return fBucketSize;}
        void            GetClosestPoints(const point_type& rRef,UInt_t nPoints,std::vector<std::pair<const _DataPoint*,Double_t> >& vFoundPoints) const;
        Double_t         GetEffectiveEntries() const;
        KDTree<_DataPoint>* GetFrozenCopy();
        UInt_t          GetNBins() const;
        UInt_t          GetEntries() const;
        void            GetPointsWithinDist(const point_type& rRef,value_type fDist,std::vector<const point_type*>& vFoundPoints) const;
        Double_t        GetTotalSumw() const;
        Double_t        GetTotalSumw2() const;
        Bool_t          Insert(const point_type& rData) {return fHead->Parent()->Insert(rData);}
        Bool_t          IsFrozen() const {return fIsFrozen;}
        iterator        Last();
        const iterator  Last() const;
        void            Print() {fHead->Parent()->Print();}
        void            Reset();
        void            SetOwner(Bool_t bIsOwner = true);
        void            SetSplitOption(eSplitOption opt);

     private:
        KDTree();
        KDTree(const KDTree<point_type>& ) {}
        KDTree<point_type>& operator=(const KDTree<point_type>& ) {return *this;}

        BaseNode*  fHead;
        Double_t   fBucketSize;
        Bool_t     fIsFrozen;
     };


  }//namespace Math
}//namespace ROOT

#include "Math/KDTree.icc"

#endif // ROOT_Math_KDTree
