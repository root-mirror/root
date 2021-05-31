#ifndef ROOT_TKDTree
#define ROOT_TKDTree

#include "TObject.h"

#include "TMath.h"
#include <vector>
// Class to pass information for multi-threading
template <typename Index, typename Value>
class TKDTreeInfo{
public:
    Int_t crow;
    Int_t cnode;
    Int_t npoints;
    Int_t cpos;
    Int_t nthread;
};

template <typename Index, typename Value> class TKDTree : public TObject
{
public:

   TKDTree();
   TKDTree(Index npoints, Index ndim, UInt_t bsize);
   TKDTree(Index npoints, Index ndim, UInt_t bsize, Value **data);
   ~TKDTree();

   void            Build();  // Old version
   void            Build(Int_t nthread);  // build the tree with multi-thread
   void            Build_relay(TKDTreeInfo<Index, Value> info, std::vector<TKDTreeInfo<Index, Value> > *threads);

   Double_t        Distance(const Value *point, Index ind, Int_t type=2) const;
   void            DistanceToNode(const Value *point, Index inode, Value &min, Value &max, Int_t type=2);

   // Get indexes of left and right daughter nodes
   Int_t   GetLeft(Int_t inode)  const    {return inode*2+1;}
   Int_t   GetRight(Int_t inode) const    {return (inode+1)*2;}
   Int_t   GetParent(Int_t inode) const  {return (inode-1)/2;}
   //
   // Other getters
   Index*  GetPointsIndexes(Int_t node) const;
   void    GetNodePointsIndexes(Int_t node, Int_t &first1, Int_t &last1, Int_t &first2, Int_t &last2) const;
   UChar_t GetNodeAxis(Int_t id) const {return (id < 0 || id >= fNNodes) ? 0 : fAxis[id];}
   Value   GetNodeValue(Int_t id) const {return (id < 0 || id >= fNNodes) ? 0 : fValue[id];}
   Int_t   GetNNodes() const {return fNNodes;}
   Int_t   GetTotalNodes() const {return fTotalNodes;}
   Value*  GetBoundaries();
   Value*  GetBoundariesExact();
   Value*  GetBoundary(const Int_t node);
   Value*  GetBoundaryExact(const Int_t node);
   Index   GetNPoints() { return fNPoints; }
   Index   GetNDim()    { return fNDim; }
   Index   GetNPointsNode(Int_t node) const;

   //Getters for internal variables.
   Int_t   GetRowT0() {return fRowT0;}      //! smallest terminal row
   Int_t   GetCrossNode() {return fCrossNode;}  //! cross node
   Int_t   GetOffset() {return fOffset;}     //! offset in fIndPoints
   Index*  GetIndPoints() {return fIndPoints;}
   Index   GetBucketSize() {return fBucketSize;}

   void    FindNearestNeighbors(const Value *point, Int_t k, Index *ind, Value *dist);
   Index   FindNode(const Value * point) const;
   void    FindPoint(Value * point, Index &index, Int_t &iter);
   void    FindInRange(Value *point, Value range, std::vector<Index> &res);
   void    FindBNodeA(Value * point, Value * delta, Int_t &inode);

   Bool_t  IsTerminal(Index inode) const {return (inode>=fNNodes);}
   Int_t   IsOwner() { return fDataOwner; }
   Value   KOrdStat(Index ntotal, Value *a, Index k, Index *index) const;


   void    MakeBoundaries(Value *range = 0x0);
   void    MakeBoundariesExact();
   void    SetData(Index npoints, Index ndim, UInt_t bsize, Value **data);
   Int_t   SetData(Index idim, Value *data);
   void    SetOwner(Int_t owner) { fDataOwner = owner; }
   void    Spread(Index ntotal, Value *a, Index *index, Value &min, Value &max) const;

 private:
   TKDTree(const TKDTree &); // not implemented
   TKDTree<Index, Value>& operator=(const TKDTree<Index, Value>&); // not implemented
   void CookBoundaries(const Int_t node, Bool_t left);

   void UpdateNearestNeighbors(Index inode, const Value *point, Int_t kNN, Index *ind, Value *dist);
   void UpdateRange(Index inode, Value *point, Value range, std::vector<Index> &res);

 protected:
   Int_t   fDataOwner;  //! 0 - not owner, 2 - owner of the pointer array, 1 - owner of the whole 2-d array
   Int_t   fNNodes;     // size of node array
   Int_t   fTotalNodes; // total number of nodes (fNNodes + terminal nodes)
   Index   fNDim;       // number of dimensions
   Index   fNDimm;      // dummy 2*fNDim
   Index   fNPoints;    // number of multidimensional points
   Index   fBucketSize; // size of the terminal nodes
   UChar_t *fAxis;      //[fNNodes] nodes cutting axis
   Value   *fValue;     //[fNNodes] nodes cutting value
   //
   Value   *fRange;     //[fNDimm] range of data for each dimension
   Value   **fData;     //! data points
   Value   *fBoundaries;//! nodes boundaries


   Index   *fIndPoints; //! array of points indexes
   Int_t   fRowT0;      //! smallest terminal row - first row that contains terminal nodes
   Int_t   fCrossNode;  //! cross node - node that begins the last row (with terminal nodes only)
   Int_t   fOffset;     //! offset in fIndPoints - if there are 2 rows, that contain terminal nodes
                        //  fOffset returns the index in the fIndPoints array of the first point
                        //  that belongs to the first node on the second row.


   ClassDef(TKDTree, 1)  // KD tree
};

typedef TKDTree<Int_t, Double_t> TKDTreeID;
typedef TKDTree<Int_t, Float_t> TKDTreeIF;

// Test functions:
TKDTreeIF *  TKDTreeTestBuild();

#endif

