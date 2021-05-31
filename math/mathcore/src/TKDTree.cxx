// @(#)root/mathcore:$Id$
// Authors: Marian Ivanov and Alexandru Bercuci 04/03/2005

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TKDTree.h"
#include "TRandom.h"

#include "TString.h"
#include <string.h>
#include <limits>

#include <ROOT/TThreadExecutor.hxx>

templateClassImp(TKDTree);


/**
\class TKDTree
\brief Class implementing a kd-tree

Contents:
1. What is kd-tree
2. How to cosntruct kdtree - Pseudo code
3. Using TKDTree
   a. Creating the kd-tree and setting the data
   b. Navigating the kd-tree
4. TKDTree implementation - technical details
   a. The order of nodes in internal arrays
   b. Division algorithm
   c. The order of nodes in boundary related arrays



### 1. What is kdtree ? ( [http://en.wikipedia.org/wiki/Kd-tree] )

In computer science, a kd-tree (short for k-dimensional tree) is a space-partitioning data structure
for organizing points in a k-dimensional space. kd-trees are a useful data structure for several
applications, such as searches involving a multidimensional search key (e.g. range searches and
nearest neighbour searches). kd-trees are a special case of BSP trees.

A kd-tree uses only splitting planes that are perpendicular to one of the coordinate system axes.
This differs from BSP trees, in which arbitrary splitting planes can be used.
In addition, in the typical definition every node of a kd-tree, from the root to the leaves, stores a point.
This differs from BSP trees, in which leaves are typically the only nodes that contain points
(or other geometric primitives). As a consequence, each splitting plane must go through one of
the points in the kd-tree. kd-trees are a variant that store data only in leaf nodes.

### 2. Constructing a classical kd-tree ( Pseudo code)

Since there are many possible ways to choose axis-aligned splitting planes, there are many different ways
to construct kd-trees. The canonical method of kd-tree construction has the following constraints:

* As one moves down the tree, one cycles through the axes used to select the splitting planes.
  (For example, the root would have an x-aligned plane, the root's children would both have y-aligned
  planes, the root's grandchildren would all have z-aligned planes, and so on.)
* At each step, the point selected to create the splitting plane is the median of the points being
  put into the kd-tree, with respect to their coordinates in the axis being used. (Note the assumption
  that we feed the entire set of points into the algorithm up-front.)

This method leads to a balanced kd-tree, in which each leaf node is about the same distance from the root.
However, balanced trees are not necessarily optimal for all applications.
The following pseudo-code illustrates this canonical construction procedure (NOTE, that the procedure used
by the TKDTree class is a bit different, the following pseudo-code is given as a simple illustration of the
concept):

~~~~
function kdtree (list of points pointList, int depth)
{
    if pointList is empty
        return nil;
    else
    {
        // Select axis based on depth so that axis cycles through all valid values
        var int axis := depth mod k;

        // Sort point list and choose median as pivot element
        select median from pointList;

        // Create node and construct subtrees
        var tree_node node;
        node.location := median;
        node.leftChild := kdtree(points in pointList before median, depth+1);
        node.rightChild := kdtree(points in pointList after median, depth+1);
        return node;
    }
}
~~~~

Our construction method is optimized to save memory, and differs a bit from the constraints above.
In particular, the division axis is chosen as the one with the biggest spread, and the point to create the
splitting plane is chosen so, that one of the two subtrees contains exactly 2^k terminal nodes and is a
perfectly balanced binary tree, and, while at the same time, trying to keep the number of terminal nodes
in the 2 subtrees as close as possible. The following section gives more details about our implementation.

### 3. Using TKDTree

#### 3a. Creating the tree and setting the data
    The interface of the TKDTree, that allows to set input data, has been developped to simplify using it
    together with TTree::Draw() functions. That's why the data has to be provided column-wise. For example:

\code{.cpp}
    {
    TTree *datatree = ...
    ...
    datatree->Draw("x:y:z", "selection", "goff");
    //now make a kd-tree on the drawn variables
    TKDTreeID *kdtree = new TKDTreeID(npoints, 3, 1);
    kdtree->SetData(0, datatree->GetV1());
    kdtree->SetData(1, datatree->GetV2());
    kdtree->SetData(2, datatree->GetV3());
    kdtree->Build();
    }
    NOTE, that this implementation of kd-tree doesn't support adding new points after the tree has been built
    Of course, it's not necessary to use TTree::Draw(). What is important, is to have data columnwise.
    An example with regular arrays:
    {
    Int_t npoints = 100000;
    Int_t ndim = 3;
    Int_t bsize = 1;
    Double_t xmin = -0.5;
    Double_t xmax = 0.5;
    Double_t *data0 = new Double_t[npoints];
    Double_t *data1 = new Double_t[npoints];
    Double_t *data2 = new Double_t[npoints];
    Double_t *y     = new Double_t[npoints];
    for (Int_t i=0; i<npoints; i++){
       data0[i]=gRandom->Uniform(xmin, xmax);
       data1[i]=gRandom->Uniform(xmin, xmax);
       data2[i]=gRandom->Uniform(xmin, xmax);
    }
    TKDTreeID *kdtree = new TKDTreeID(npoints, ndim, bsize);
    kdtree->SetData(0, data0);
    kdtree->SetData(1, data1);
    kdtree->SetData(2, data2);
    kdtree->Build();
    }
\endcode

    By default, the kd-tree doesn't own the data and doesn't delete it with itself. If you want the
    data to be deleted together with the kd-tree, call TKDTree::SetOwner(kTRUE).

    Most functions of the kd-tree don't require the original data to be present after the tree
    has been built. Check the functions documentation for more details.

#### 3b. Navigating the kd-tree

    Nodes of the tree are indexed top to bottom, left to right. The root node has index 0. Functions
    TKDTree::GetLeft(Index inode), TKDTree::GetRight(Index inode) and TKDTree::GetParent(Index inode)
    allow to find the children and the parent of a given node.

    For a given node, one can find the indexes of the original points, contained in this node,
    by calling the GetNodePointsIndexes(Index inode) function. Additionally, for terminal nodes,
    there is a function GetPointsIndexes(Index inode) that returns a pointer to the relevant
    part of the index array. To find the number of point in the node
    (not only terminal), call TKDTree::GetNpointsNode(Index inode).

### 4.  TKDtree implementation details - internal information, not needed to use the kd-tree.

####  4a. Order of nodes in the node information arrays:

TKDtree is optimized to minimize memory consumption.
Nodes of the TKDTree do not store pointers to the left and right children or to the parent node,
but instead there are several 1-d arrays of size fNNodes with information about the nodes.
The order of the nodes information in the arrays is described below. It's important to understand
it, if one's class needs to store some kind of additional information on the per node basis, for
example, the fit function parameters.

- Drawback:   Insertion to the TKDtree is not supported.
- Advantage:  Random access is supported

As noted above, the construction of the kd-tree involves choosing the axis and the point on
that axis to divide the remaining points approximately in half. The exact algorithm for choosing
the division point is described in the next section. The sequence of divisions is
recorded in the following arrays:
~~~~
fAxix[fNNodes]  - Division axis (0,1,2,3 ...)
fValue[fNNodes] - Division value
~~~~

Given the index of a node in those arrays, it's easy to find the indices, corresponding to
children nodes or the parent node:
Suppose, the parent node is stored under the index inode. Then:
- Left child `index  = inode*2+1`
- Right child `index =  (inode+1)*2`

Suppose, that the child node is stored under the index inode. Then:
- Parent `index = inode/2`

Number of division nodes and number of terminals :
`fNNodes = (fNPoints/fBucketSize)`

The nodes are filled always from left side to the right side:
Let inode be the index of a node, and irow - the index of a row
The TKDTree looks the following way:
Ideal case:
~~~~
Number of _terminal_ nodes = 2^N,  N=3

           INode
irow 0     0                                                                   -  1 inode
irow 1     1                              2                                    -  2 inodes
irow 2     3              4               5               6                    -  4 inodes
irow 3     7       8      9      10       11     12       13      14           -  8 inodes
~~~~

Non ideal case:
~~~~
Number of _terminal_ nodes = 2^N+k,  N=3  k=1

          INode
irow 0     0                                                                   - 1 inode
irow 1     1                              2                                    - 2 inodes
irow 2     3              4               5               6                    - 3 inodes
irow 3     7       8      9      10       11     12       13      14           - 8 inodes
irow 4     15  16                                                              - 2 inodes
~~~~

#### 4b. The division algorithm:

As described above, the kd-tree is built by repeatingly dividing the given set of points into
2 smaller sets. The cut is made on the axis with the biggest spread, and the value on the axis,
on which the cut is performed, is chosen based on the following formula:
Suppose, we want to divide n nodes into 2 groups, left and right. Then the left and right
will have the following number of nodes:

~~~~
n=2^k+rest

Left  = 2^k-1 +  ((rest>2^k-2) ?  2^k-2      : rest)
Right = 2^k-1 +  ((rest>2^k-2) ?  rest-2^k-2 : 0)
~~~~

For example, let `n_nodes=67`. Then, the closest `2^k=64, 2^k-1=32, 2^k-2=16`.
Left node gets `32+3=35` sub-nodes, and the right node gets 32 sub-nodes

The division process continues until all the nodes contain not more than a predefined number
of points.

#### 4c. The order of nodes in boundary-related arrays

Some kd-tree based algorithms need to know the boundaries of each node. This information can
be computed by calling the TKDTree::MakeBoundaries() function. It fills the following arrays:

- `fRange` : array containing the boundaries of the domain:
      `| 1st dimension (min + max) | 2nd dimension (min + max) | ...`
`fBoundaries` : nodes boundaries
      `| 1st node {1st dim * 2 elements | 2nd dim * 2 elements | ...} | 2nd node {...} | ...`


The nodes are arranged in the order described in section 3a.


- **Note**: the storage of the TKDTree in a file which include also the contained data is not
      supported. One must store the data separatly in a file (e.g. using a TTree) and then
      re-creating the TKDTree from the data, after having read them from the file

@ingroup Random


*/
////////////////////////////////////////////////////////////////////////////////
/// Default constructor. Nothing is built

template <typename  Index, typename Value>
TKDTree<Index, Value>::TKDTree() :
   TObject()
   ,fDataOwner(kFALSE)
   ,fNNodes(0)
   ,fTotalNodes(0)
   ,fNDim(0)
   ,fNDimm(0)
   ,fNPoints(0)
   ,fBucketSize(0)
   ,fAxis(0x0)
   ,fValue(0x0)
   ,fRange(0x0)
   ,fData(0x0)
   ,fBoundaries(0x0)
   ,fIndPoints(0x0)
   ,fRowT0(0)
   ,fCrossNode(0)
   ,fOffset(0)
{
}

template <typename  Index, typename Value>
TKDTree<Index, Value>::TKDTree(Index npoints, Index ndim, UInt_t bsize) :
   TObject()
   ,fDataOwner(0)
   ,fNNodes(0)
   ,fTotalNodes(0)
   ,fNDim(ndim)
   ,fNDimm(2*ndim)
   ,fNPoints(npoints)
   ,fBucketSize(bsize)
   ,fAxis(0x0)
   ,fValue(0x0)
   ,fRange(0x0)
   ,fData(0x0)
   ,fBoundaries(0x0)
   ,fIndPoints(0x0)
   ,fRowT0(0)
   ,fCrossNode(0)
   ,fOffset(0)
{
// Create the kd-tree of npoints from ndim-dimensional space. Parameter bsize stands for the
// maximal number of points in the terminal nodes (buckets).
// Proceed by calling one of the SetData() functions and then the Build() function
// Note, that updating the tree with new data after the Build() function has been called is
// not possible!

}

////////////////////////////////////////////////////////////////////////////////
/// Create a kd-tree from the provided data array. This function only sets the data,
/// call Build() to build the tree!!!
/// Parameteres:
/// - npoints - total number of points. Adding points after the tree is built is not supported
/// - ndim    - number of dimensions
/// - bsize   - maximal number of points in the terminal nodes
/// - data    - the data array
///
/// The data should be placed columnwise (like in a TTree).
/// The columnwise orientation is chosen to simplify the usage together with TTree::GetV1() like functions.
/// An example of filling such an array for a 2d case:
/// Double_t **data = new Double_t*[2];
/// data[0] = new Double_t[npoints];
/// data[1] = new Double_t[npoints];
/// for (Int_t i=0; i<npoints; i++){
///    data[0][i]=gRandom->Uniform(-1, 1); //fill the x-coordinate
///    data[1][i]=gRandom->Uniform(-1, 1); //fill the y-coordinate
/// }
///
/// By default, the kd-tree doesn't own the data. If you want the kd-tree to delete the data array, call
/// kdtree->SetOwner(kTRUE).
///

template <typename  Index, typename Value>
TKDTree<Index, Value>::TKDTree(Index npoints, Index ndim, UInt_t bsize, Value **data) :
   TObject()
   ,fDataOwner(0)
   ,fNNodes(0)
   ,fTotalNodes(0)
   ,fNDim(ndim)
   ,fNDimm(2*ndim)
   ,fNPoints(npoints)
   ,fBucketSize(bsize)
   ,fAxis(0x0)
   ,fValue(0x0)
   ,fRange(0x0)
   ,fData(data) //Columnwise!!!!!
   ,fBoundaries(0x0)
   ,fIndPoints(0x0)
   ,fRowT0(0)
   ,fCrossNode(0)
   ,fOffset(0)
{

   //Build();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor
/// By default, the original data is not owned by kd-tree and is not deleted with it.
/// If you want to delete the data along with the kd-tree, call SetOwner(kTRUE).

template <typename  Index, typename Value>
TKDTree<Index, Value>::~TKDTree()
{
   if (fAxis) delete [] fAxis;
   if (fValue) delete [] fValue;
   if (fIndPoints) delete [] fIndPoints;
   if (fRange) delete [] fRange;
   if (fBoundaries) delete [] fBoundaries;
   if (fData) {
      if (fDataOwner==1){
         //the tree owns all the data
         for(int idim=0; idim<fNDim; idim++) delete [] fData[idim];
      }
      if (fDataOwner>0) {
         //the tree owns the array of pointers
         delete [] fData;
      }
   }
//   if (fDataOwner && fData){
//      for(int idim=0; idim<fNDim; idim++) delete [] fData[idim];
//      delete [] fData;
//   }
}


////////////////////////////////////////////////////////////////////////////////
///
/// Build the kd-tree
///
/// 1. calculate number of nodes
/// 2. calculate first terminal row
/// 3. initialize index array
/// 4. non recursive building of the binary tree
///
///
/// The tree is divided recursively. See class description, section 4b for the details
/// of the division alogrithm

template <typename Index, typename Value>
void TKDTree<Index, Value>::Build_relay(TKDTreeInfo<Index, Value> info, std::vector<TKDTreeInfo<Index, Value> > *threads) {
   Int_t crow = info.crow;
   Int_t cnode = info.cnode;
   Int_t npoints = info.npoints;
   Int_t cpos = info.cpos;
   Int_t nthread = info.nthread;

   //    Distribute thread before reach the setting
   if( (1<<(crow+1)) <= nthread )
   {
      // Same as normal, just different at the end.
      // Can be replaced by build_atom if choose a decent way to return nleft and nright together.(not done yet)
      if (npoints <= fBucketSize) {
         return;
      }

      Int_t nbuckets0 = npoints/fBucketSize;           //current number of  buckets
      if (npoints%fBucketSize) nbuckets0++;            //
      Int_t restRows = fRowT0-crow;  // rest of fully occupied node row
      if (restRows<0) restRows =0;
      for (;nbuckets0>(2<<restRows); restRows++) {}
      Int_t nfull = 1<<restRows;
      Int_t nrest = nbuckets0-nfull;
      Int_t nleft =0, nright =0;
      //
      if (nrest>(nfull/2)){
         nleft  = nfull*fBucketSize;
         nright = npoints-nleft;
      }else{
         nright = nfull*fBucketSize/2;
         nleft  = npoints-nright;
      }

      //
      //find the axis with biggest spread
      Value maxspread=0;
      Value tempspread, min, max;
      Index axspread=0;
      Value *array;
      for (Int_t idim=0; idim<fNDim; idim++){
         array = fData[idim];
         Spread(npoints, array, fIndPoints+cpos, min, max);
         tempspread = max - min;
         if (maxspread < tempspread) {
            maxspread=tempspread;
            axspread = idim;
         }
         if(cnode) continue;
         //printf("set %d %6.3f %6.3f\n", idim, min, max);
         fRange[2*idim] = min; fRange[2*idim+1] = max;
      }
      array = fData[axspread];
      KOrdStat(npoints, array, nleft, fIndPoints+cpos);
      fAxis[cnode]  = axspread;
      fValue[cnode] = array[fIndPoints[cpos+nleft]];

      //cout << "threading: cnode " << cnode << " npoints " << npoints <<" cpos "<< cpos << endl;
      TKDTreeInfo<Index, Value> left={crow+1, cnode*2+1, nleft, cpos, nthread};
      TKDTreeInfo<Index, Value> right={crow+1, cnode*2+2, nright, cpos+nleft, nthread};

      Build_relay(left, threads);
      Build_relay(right, threads);
   }else
   {
      //Record information of nodes for workers to do.
      // printf("%d: %d end.\n",nthread,crow);
      (*threads).push_back(info);
   }
   return;
}

template <typename  Index, typename Value>
void TKDTree<Index, Value>::Build(Int_t nthread)//nthread = 2^n n=0,1,2,3...
{
   //1.
   fNNodes = fNPoints/fBucketSize-1;
   if (fNPoints%fBucketSize) fNNodes++;
   fTotalNodes = fNNodes + fNPoints/fBucketSize + ((fNPoints%fBucketSize)?1:0);
   //2.
   fRowT0=0;
   for ( ;(fNNodes+1)>(1<<fRowT0);fRowT0++) {}
   fRowT0-=1;
   //         2 = 2**0 + 1
   //         3 = 2**1 + 1
   //         4 = 2**1 + 2
   //         5 = 2**2 + 1
   //         6 = 2**2 + 2
   //         7 = 2**2 + 3
   //         8 = 2**2 + 4

   //3.
   // allocate space for boundaries
   fRange = new Value[2*fNDim];
   fIndPoints= new Index[fNPoints];
   for (Index i=0; i<fNPoints; i++) fIndPoints[i] = i;
   fAxis  = new UChar_t[fNNodes];
   fValue = new Value[fNNodes];
   //
   fCrossNode = (1<<(fRowT0+1))-1;
   if (fCrossNode<fNNodes) fCrossNode = 2*fCrossNode+1;
   //
   //  fOffset = (((fNNodes+1)-(1<<fRowT0)))*2;
   Int_t   over   = (fNNodes+1)-(1<<fRowT0);
   Int_t   filled = ((1<<fRowT0)-over)*fBucketSize;
   fOffset = fNPoints-filled;
   //Pre-setting finished

   //Distribute
   std::vector<TKDTreeInfo<Index, Value> > threads;

   TKDTreeInfo<Index, Value> root = {0, 0, fNPoints, 0, nthread};
   Build_relay(root, &threads);

   //This lambda is worker doing Build_atom down-ward to the leaves
   auto Build_sub = [&](TKDTreeInfo<Index, Value> info)
   {
      Int_t crow = info.crow;
      Int_t cnode = info.cnode;
      Int_t npoints = info.npoints;
      Int_t cpos = info.cpos;
      //    Function describe:
      //        Start expanse the KDTree node with given status(node).
      //
      //    stack for non recursive build - size 128 bytes enough
      Int_t rowStack[128];
      Int_t nodeStack[128];
      Int_t npointStack[128];
      Int_t posStack[128];
      Int_t currentIndex = 0;

      //    Insert current status node as the first one
      rowStack[0]    = crow;
      nodeStack[0]   = cnode;
      npointStack[0] = npoints;
      posStack[0]   = cpos;

      while (currentIndex>=0){

         npoints  = npointStack[currentIndex];
         if (npoints<=fBucketSize) {
            currentIndex--;
            continue; // terminal node
         }
         crow     = rowStack[currentIndex];
         cpos     = posStack[currentIndex];
         cnode    = nodeStack[currentIndex];
         //printf("currentIndex %d npoints %d node %d\n", currentIndex, npoints, cnode);
         //
         // divide points
         Int_t nbuckets0 = npoints/fBucketSize;           //current number of  buckets
         if (npoints%fBucketSize) nbuckets0++;            //
         Int_t restRows = fRowT0-rowStack[currentIndex];  // rest of fully occupied node row
         if (restRows<0) restRows =0;
         for (;nbuckets0>(2<<restRows); restRows++) {}
         Int_t nfull = 1<<restRows;
         Int_t nrest = nbuckets0-nfull;
         Int_t nleft =0, nright =0;
         //
         if (nrest>(nfull/2)){
            nleft  = nfull*fBucketSize;
            nright = npoints-nleft;
         }else{
            nright = nfull*fBucketSize/2;
            nleft  = npoints-nright;
         }

         //
         //find the axis with biggest spread
         Value maxspread=0;
         Value tempspread, min, max;
         Index axspread=0;
         Value *array;
         for (Int_t idim=0; idim<fNDim; idim++){
            array = fData[idim];
            Spread(npoints, array, fIndPoints+cpos, min, max);
            tempspread = max - min;
            if (maxspread < tempspread) {
               maxspread=tempspread;
               axspread = idim;
            }
            if(cnode) continue;
            //printf("set %d %6.3f %6.3f\n", idim, min, max);
            fRange[2*idim] = min; fRange[2*idim+1] = max;
         }
         array = fData[axspread];
         KOrdStat(npoints, array, nleft, fIndPoints+cpos);
         fAxis[cnode]  = axspread;
         fValue[cnode] = array[fIndPoints[cpos+nleft]];
         //printf("Set node %d : ax %d val %f\n", cnode, node->fAxis, node->fValue);
         //
         //
         npointStack[currentIndex] = nleft;
         rowStack[currentIndex]    = crow+1;
         posStack[currentIndex]    = cpos;
         nodeStack[currentIndex]   = cnode*2+1;
         currentIndex++;
         npointStack[currentIndex] = nright;
         rowStack[currentIndex]    = crow+1;
         posStack[currentIndex]    = cpos+nleft;
         nodeStack[currentIndex]   = (cnode*2)+2;
      }
      return 0;
   };

   ROOT::TThreadExecutor pool;
   //printf("%d distrubting %d.\n",nthread,threads.size());
   auto results = pool.Map(Build_sub, threads);
}

template <typename  Index, typename Value>
void TKDTree<Index, Value>::Build()
{
   //1.
   fNNodes = fNPoints/fBucketSize-1;
   if (fNPoints%fBucketSize) fNNodes++;
   fTotalNodes = fNNodes + fNPoints/fBucketSize + ((fNPoints%fBucketSize)?1:0);
   //2.
   fRowT0=0;
   for ( ;(fNNodes+1)>(1<<fRowT0);fRowT0++) {}
   fRowT0-=1;
   //         2 = 2**0 + 1
   //         3 = 2**1 + 1
   //         4 = 2**1 + 2
   //         5 = 2**2 + 1
   //         6 = 2**2 + 2
   //         7 = 2**2 + 3
   //         8 = 2**2 + 4

   //3.
   // allocate space for boundaries
   fRange = new Value[2*fNDim];
   fIndPoints= new Index[fNPoints];
   for (Index i=0; i<fNPoints; i++) fIndPoints[i] = i;
   fAxis  = new UChar_t[fNNodes];
   fValue = new Value[fNNodes];
   //
   fCrossNode = (1<<(fRowT0+1))-1;
   if (fCrossNode<fNNodes) fCrossNode = 2*fCrossNode+1;
   //
   //  fOffset = (((fNNodes+1)-(1<<fRowT0)))*2;
   Int_t   over   = (fNNodes+1)-(1<<fRowT0);
   Int_t   filled = ((1<<fRowT0)-over)*fBucketSize;
   fOffset = fNPoints-filled;

   //
   //    printf("Row0      %d\n", fRowT0);
   //    printf("CrossNode %d\n", fCrossNode);
   //    printf("Offset    %d\n", fOffset);
   //
   //
   //4.
   //    stack for non recursive build - size 128 bytes enough
   Int_t rowStack[128];
   Int_t nodeStack[128];
   Int_t npointStack[128];
   Int_t posStack[128];
   Int_t currentIndex = 0;
   Int_t iter =0;
   rowStack[0]    = 0;
   nodeStack[0]   = 0;
   npointStack[0] = fNPoints;
   posStack[0]   = 0;
   //
   Int_t nbucketsall =0;
   while (currentIndex>=0){
      iter++;
      //
      Int_t npoints  = npointStack[currentIndex];
      if (npoints<=fBucketSize) {
         //printf("terminal node : index %d iter %d\n", currentIndex, iter);
         currentIndex--;
         nbucketsall++;
         continue; // terminal node
      }
      Int_t crow     = rowStack[currentIndex];
      Int_t cpos     = posStack[currentIndex];
      Int_t cnode    = nodeStack[currentIndex];
      //printf("currentIndex %d npoints %d node %d\n", currentIndex, npoints, cnode);
      //
      // divide points
      Int_t nbuckets0 = npoints/fBucketSize;           //current number of  buckets
      if (npoints%fBucketSize) nbuckets0++;            //
      Int_t restRows = fRowT0-rowStack[currentIndex];  // rest of fully occupied node row
      if (restRows<0) restRows =0;
      for (;nbuckets0>(2<<restRows); restRows++) {}
      Int_t nfull = 1<<restRows;
      Int_t nrest = nbuckets0-nfull;
      Int_t nleft =0, nright =0;
      //
      if (nrest>(nfull/2)){
         nleft  = nfull*fBucketSize;
         nright = npoints-nleft;
      }else{
         nright = nfull*fBucketSize/2;
         nleft  = npoints-nright;
      }

      //
      //find the axis with biggest spread
      Value maxspread=0;
      Value tempspread, min, max;
      Index axspread=0;
      Value *array;
      for (Int_t idim=0; idim<fNDim; idim++){
         array = fData[idim];
         Spread(npoints, array, fIndPoints+cpos, min, max);
         tempspread = max - min;
         if (maxspread < tempspread) {
            maxspread=tempspread;
            axspread = idim;
         }
         if(cnode) continue;
         //printf("set %d %6.3f %6.3f\n", idim, min, max);
         fRange[2*idim] = min; fRange[2*idim+1] = max;
      }
      array = fData[axspread];
      KOrdStat(npoints, array, nleft, fIndPoints+cpos);
      fAxis[cnode]  = axspread;
      fValue[cnode] = array[fIndPoints[cpos+nleft]];
      //printf("Set node %d : ax %d val %f\n", cnode, node->fAxis, node->fValue);
      //
      //
      npointStack[currentIndex] = nleft;
      rowStack[currentIndex]    = crow+1;
      posStack[currentIndex]    = cpos;
      nodeStack[currentIndex]   = cnode*2+1;
      currentIndex++;
      npointStack[currentIndex] = nright;
      rowStack[currentIndex]    = crow+1;
      posStack[currentIndex]    = cpos+nleft;
      nodeStack[currentIndex]   = (cnode*2)+2;
      //
      if (0){
         // consistency check
         Info("Build()", "%s", Form("points %d left %d right %d", npoints, nleft, nright));
         if (nleft<nright) Warning("Build", "Problem Left-Right");
         if (nleft<0 || nright<0) Warning("Build()", "Problem Negative number");
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
///Find kNN nearest neighbors to the point in the first argument
///Returns 1 on success, 0 on failure
///Arrays ind and dist are provided by the user and are assumed to be at least kNN elements long

template <typename  Index, typename Value>
void TKDTree<Index, Value>::FindNearestNeighbors(const Value *point, const Int_t kNN, Index *ind, Value *dist)
{

   if (!ind || !dist) {
      Error("FindNearestNeighbors", "Working arrays must be allocated by the user!");
      return;
   }
   for (Int_t i=0; i<kNN; i++){
      dist[i]=std::numeric_limits<Value>::max();
      ind[i]=-1;
   }
   MakeBoundariesExact();
   UpdateNearestNeighbors(0, point, kNN, ind, dist);

}

////////////////////////////////////////////////////////////////////////////////
///Update the nearest neighbors values by examining the node inode

template <typename Index, typename Value>
void TKDTree<Index, Value>::UpdateNearestNeighbors(Index inode, const Value *point, Int_t kNN, Index *ind, Value *dist)
{
   Value min=0;
   Value max=0;
   DistanceToNode(point, inode, min, max);
   if (min > dist[kNN-1]){
      //there are no closer points in this node
      return;
   }
   if (IsTerminal(inode)) {
      //examine points one by one
      Index f1, l1, f2, l2;
      GetNodePointsIndexes(inode, f1, l1, f2, l2);
      for (Int_t ipoint=f1; ipoint<=l1; ipoint++){
         Double_t d = Distance(point, fIndPoints[ipoint]);
         if (d<dist[kNN-1]){
            //found a closer point
            Int_t ishift=0;
            while(ishift<kNN && d>dist[ishift])
               ishift++;
            //replace the neighbor #ishift with the found point
            //and shift the rest 1 index value to the right
            for (Int_t i=kNN-1; i>ishift; i--){
               dist[i]=dist[i-1];
               ind[i]=ind[i-1];
            }
            dist[ishift]=d;
            ind[ishift]=fIndPoints[ipoint];
         }
      }
      return;
   }
   if (point[fAxis[inode]]<fValue[inode]){
      //first examine the node that contains the point
      UpdateNearestNeighbors(GetLeft(inode), point, kNN, ind, dist);
      UpdateNearestNeighbors(GetRight(inode), point, kNN, ind, dist);
   } else {
      UpdateNearestNeighbors(GetRight(inode), point, kNN, ind, dist);
      UpdateNearestNeighbors(GetLeft(inode), point, kNN, ind, dist);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Find the distance between point of the first argument and the point at index value ind
///Type argument specifies the metric: type=2 - L2 metric, type=1 - L1 metric

template <typename Index, typename Value>
Double_t TKDTree<Index, Value>::Distance(const Value *point, Index ind, Int_t type) const
{
   Double_t dist = 0;
   if (type==2){
      for (Int_t idim=0; idim<fNDim; idim++){
         dist+=(point[idim]-fData[idim][ind])*(point[idim]-fData[idim][ind]);
      }
      return TMath::Sqrt(dist);
   } else {
      for (Int_t idim=0; idim<fNDim; idim++){
         dist+=TMath::Abs(point[idim]-fData[idim][ind]);
      }

      return dist;
   }
   return -1;

}

////////////////////////////////////////////////////////////////////////////////
///Find the minimal and maximal distance from a given point to a given node.
///Type argument specifies the metric: type=2 - L2 metric, type=1 - L1 metric
///If the point is inside the node, both min and max are set to 0.

template <typename Index, typename Value>
void TKDTree<Index, Value>::DistanceToNode(const Value *point, Index inode, Value &min, Value &max, Int_t type)
{
   Value *bound = GetBoundaryExact(inode);
   min = 0;
   max = 0;
   Double_t dist1, dist2;

   if (type==2){
      for (Int_t idim=0; idim<fNDimm; idim+=2){
         dist1 = (point[idim/2]-bound[idim])*(point[idim/2]-bound[idim]);
         dist2 = (point[idim/2]-bound[idim+1])*(point[idim/2]-bound[idim+1]);
         //min+=TMath::Min(dist1, dist2);
         if (point[idim/2]<bound[idim] || point[idim/2]>bound[idim+1])
            min+= (dist1>dist2)? dist2 : dist1;
         // max+=TMath::Max(dist1, dist2);
         max+= (dist1>dist2)? dist1 : dist2;
      }
      min = TMath::Sqrt(min);
      max = TMath::Sqrt(max);
   } else {
      for (Int_t idim=0; idim<fNDimm; idim+=2){
         dist1 = TMath::Abs(point[idim/2]-bound[idim]);
         dist2 = TMath::Abs(point[idim/2]-bound[idim+1]);
         //min+=TMath::Min(dist1, dist2);
         min+= (dist1>dist2)? dist2 : dist1;
         // max+=TMath::Max(dist1, dist2);
         max+= (dist1>dist2)? dist1 : dist2;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// returns the index of the terminal node to which point belongs
/// (index in the fAxis, fValue, etc arrays)
/// returns -1 in case of failure

template <typename  Index, typename Value>
Index TKDTree<Index, Value>::FindNode(const Value * point) const
{
   Index stackNode[128], inode;
   Int_t currentIndex =0;
   stackNode[0] = 0;
   while (currentIndex>=0){
      inode    = stackNode[currentIndex];
      if (IsTerminal(inode)) return inode;

      currentIndex--;
      if (point[fAxis[inode]]<=fValue[inode]){
         currentIndex++;
         stackNode[currentIndex]=(inode<<1)+1; //GetLeft()
      }
      if (point[fAxis[inode]]>=fValue[inode]){
         currentIndex++;
         stackNode[currentIndex]=(inode+1)<<1; //GetRight()
      }
   }

   return -1;
}



////////////////////////////////////////////////////////////////////////////////
///
/// find the index of point
/// works only if we keep fData pointers

template <typename  Index, typename Value>
void TKDTree<Index, Value>::FindPoint(Value * point, Index &index, Int_t &iter){
   Int_t stackNode[128];
   Int_t currentIndex =0;
   stackNode[0] = 0;
   iter =0;
   //
   while (currentIndex>=0){
      iter++;
      Int_t inode    = stackNode[currentIndex];
      currentIndex--;
      if (IsTerminal(inode)){
         // investigate terminal node
         Int_t indexIP  = (inode >= fCrossNode) ? (inode-fCrossNode)*fBucketSize : (inode-fNNodes)*fBucketSize+fOffset;
         printf("terminal %d indexP %d\n", inode, indexIP);
         for (Int_t ibucket=0;ibucket<fBucketSize;ibucket++){
            Bool_t isOK    = kTRUE;
            indexIP+=ibucket;
            printf("ibucket %d index %d\n", ibucket, indexIP);
            if (indexIP>=fNPoints) continue;
            Int_t index0   = fIndPoints[indexIP];
            for (Int_t idim=0;idim<fNDim;idim++) if (fData[idim][index0]!=point[idim]) isOK = kFALSE;
            if (isOK) index = index0;
         }
         continue;
      }

      if (point[fAxis[inode]]<=fValue[inode]){
         currentIndex++;
         stackNode[currentIndex]=(inode*2)+1;
      }
      if (point[fAxis[inode]]>=fValue[inode]){
         currentIndex++;
         stackNode[currentIndex]=(inode*2)+2;
    }
  }
  //
  //  printf("Iter\t%d\n",iter);
}

////////////////////////////////////////////////////////////////////////////////
///Find all points in the sphere of a given radius "range" around the given point
///1st argument - the point
///2nd argument - radius of the shere
///3rd argument - a vector, in which the results will be returned

template <typename  Index, typename Value>
void TKDTree<Index, Value>::FindInRange(Value * point, Value range, std::vector<Index> &res)
{
   MakeBoundariesExact();
   UpdateRange(0, point, range, res);
}

////////////////////////////////////////////////////////////////////////////////
///Internal recursive function with the implementation of range searches

template <typename  Index, typename Value>
void TKDTree<Index, Value>::UpdateRange(Index inode, Value* point, Value range, std::vector<Index> &res)
{
   Value min, max;
   DistanceToNode(point, inode, min, max);
   if (min>range) {
      //all points of this node are outside the range
      return;
   }
   if (max<range && max>0) {
      //all points of this node are inside the range

      Index f1, l1, f2, l2;
      GetNodePointsIndexes(inode, f1, l1, f2, l2);

      for (Int_t ipoint=f1; ipoint<=l1; ipoint++){
         res.push_back(fIndPoints[ipoint]);
      }
      for (Int_t ipoint=f2; ipoint<=l2; ipoint++){
         res.push_back(fIndPoints[ipoint]);
      }
      return;
   }

   //this node intersects with the range
   if (IsTerminal(inode)){
      //examine the points one by one
      Index f1, l1, f2, l2;
      Double_t d;
      GetNodePointsIndexes(inode, f1, l1, f2, l2);
      for (Int_t ipoint=f1; ipoint<=l1; ipoint++){
         d = Distance(point, fIndPoints[ipoint]);
         if (d <= range){
            res.push_back(fIndPoints[ipoint]);
         }
      }
      return;
   }
   if (point[fAxis[inode]]<=fValue[inode]){
      //first examine the node that contains the point
      UpdateRange(GetLeft(inode),point, range, res);
      UpdateRange(GetRight(inode),point, range, res);
   } else {
      UpdateRange(GetLeft(inode),point, range, res);
      UpdateRange(GetRight(inode),point, range, res);
   }
}

////////////////////////////////////////////////////////////////////////////////
///return the indices of the points in that terminal node
///for all the nodes except last, the size is fBucketSize
///for the last node it's fOffset%fBucketSize

template <typename Index, typename Value>
Index*  TKDTree<Index, Value>::GetPointsIndexes(Int_t node) const
{
   if (!IsTerminal(node)){
      printf("GetPointsIndexes() only for terminal nodes, use GetNodePointsIndexes() instead\n");
      return 0;
   }
   Int_t offset = (node >= fCrossNode) ? (node-fCrossNode)*fBucketSize : fOffset+(node-fNNodes)*fBucketSize;
   return &fIndPoints[offset];
}

////////////////////////////////////////////////////////////////////////////////
///Return the indices of points in that node
///Indices are returned as the first and last value of the part of indices array, that belong to this node
///Sometimes points are in 2 intervals, then the first and last value for the second one are returned in
///third and fourth parameter, otherwise first2 is set to 0 and last2 is set to -1
///To iterate over all the points of the node #inode, one can do, for example:
///Index *indices = kdtree->GetPointsIndexes();
///Int_t first1, last1, first2, last2;
///kdtree->GetPointsIndexes(inode, first1, last1, first2, last2);
///for (Int_t ipoint=first1; ipoint<=last1; ipoint++){
///   point = indices[ipoint];
///   //do something with point;
///}
///for (Int_t ipoint=first2; ipoint<=last2; ipoint++){
///   point = indices[ipoint];
///   //do something with point;
///}

template <typename Index, typename Value>
void  TKDTree<Index, Value>::GetNodePointsIndexes(Int_t node, Int_t &first1, Int_t &last1, Int_t &first2, Int_t &last2) const
{

   if (IsTerminal(node)){
      //the first point in the node is computed by the following formula:
      Index offset = (node >= fCrossNode) ? (node-fCrossNode)*fBucketSize : fOffset+(node-fNNodes)*fBucketSize;
      first1 = offset;
      last1 = offset + GetNPointsNode(node)-1;
      first2 = 0;
      last2 = -1;
      return;
   }

   Index firsttermnode = fNNodes;
   Index ileft = node;
   Index iright = node;
   Index f1, l1, f2, l2;
//this is the left-most node
   while (ileft<firsttermnode)
      ileft = GetLeft(ileft);
//this is the right-most node
   while (iright<firsttermnode)
      iright = GetRight(iright);

   if (ileft>iright){
//      first1 = firsttermnode;
//      last1 = iright;
//      first2 = ileft;
//      last2 = fTotalNodes-1;
      GetNodePointsIndexes(firsttermnode, f1, l1, f2, l2);
      first1 = f1;
      GetNodePointsIndexes(iright, f1, l1, f2, l2);
      last1 = l1;
      GetNodePointsIndexes(ileft, f1, l1, f2, l2);
      first2 = f1;
      GetNodePointsIndexes(fTotalNodes-1, f1, l1, f2, l2);
      last2 = l1;

   }  else {
      GetNodePointsIndexes(ileft, f1, l1, f2, l2);
      first1 = f1;
      GetNodePointsIndexes(iright, f1, l1, f2, l2);
      last1 = l1;
      first2 = 0;
      last2 = -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Get number of points in this node
///for all the terminal nodes except last, the size is fBucketSize
///for the last node it's fOffset%fBucketSize, or if fOffset%fBucketSize==0, it's also fBucketSize

template <typename Index, typename Value>
Index TKDTree<Index, Value>::GetNPointsNode(Int_t inode) const
{
   if (IsTerminal(inode)){

      if (inode!=fTotalNodes-1) return fBucketSize;
      else {
         if (fOffset%fBucketSize==0) return fBucketSize;
         else return fOffset%fBucketSize;
      }
   }

   Int_t f1, l1, f2, l2;
   GetNodePointsIndexes(inode, f1, l1, f2, l2);
   Int_t sum = l1-f1+1;
   sum += l2-f2+1;
   return sum;
}


////////////////////////////////////////////////////////////////////////////////
/// Set the data array. See the constructor function comments for details

template <typename  Index, typename Value>
void TKDTree<Index, Value>::SetData(Index npoints, Index ndim, UInt_t bsize, Value **data)
{
// TO DO
//
// Check reconstruction/reallocation of memory of data. Maybe it is not
// necessary to delete and realocate space but only to use the same space

   Clear();

   //Columnwise!!!!
   fData = data;
   fNPoints = npoints;
   fNDim = ndim;
   fBucketSize = bsize;

   Build();
}

////////////////////////////////////////////////////////////////////////////////
///Set the coordinate #ndim of all points (the column #ndim of the data matrix)
///After setting all the data columns, proceed by calling Build() function
///Note, that calling this function after Build() is not possible
///Note also, that no checks on the array sizes is performed anywhere

template <typename  Index, typename Value>
Int_t TKDTree<Index, Value>::SetData(Index idim, Value *data)
{
   if (fAxis || fValue) {
      Error("SetData", "The tree has already been built, no updates possible");
      return 0;
   }

   if (!fData) {
      fData = new Value*[fNDim];
   }
   fData[idim]=data;
   fDataOwner = 2;
   return 1;
}


////////////////////////////////////////////////////////////////////////////////
///Calculate spread of the array a

template <typename  Index, typename Value>
void TKDTree<Index, Value>::Spread(Index ntotal, Value *a, Index *index, Value &min, Value &max) const
{
   Index i;
   min = a[index[0]];
   max = a[index[0]];
   for (i=0; i<ntotal; i++){
      if (a[index[i]]<min) min = a[index[i]];
      if (a[index[i]]>max) max = a[index[i]];
   }
}


////////////////////////////////////////////////////////////////////////////////
///
///copy of the TMath::KOrdStat because I need an Index work array

template <typename  Index, typename Value>
Value TKDTree<Index, Value>::KOrdStat(Index ntotal, Value *a, Index k, Index *index) const
{
   Index i, ir, j, l, mid;
   Index arr;
   Index temp;

   Index rk = k;
   l=0;
   ir = ntotal-1;
   for(;;) {
      if (ir<=l+1) { //active partition contains 1 or 2 elements
         if (ir == l+1 && a[index[ir]]<a[index[l]])
         {temp = index[l]; index[l]=index[ir]; index[ir]=temp;}
         Value tmp = a[index[rk]];
         return tmp;
      } else {
         mid = (l+ir) >> 1; //choose median of left, center and right
         {temp = index[mid]; index[mid]=index[l+1]; index[l+1]=temp;}//elements as partitioning element arr.
          if (a[index[l]]>a[index[ir]])  //also rearrange so that a[l]<=a[l+1]
          {temp = index[l]; index[l]=index[ir]; index[ir]=temp;}

          if (a[index[l+1]]>a[index[ir]])
          {temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;}

          if (a[index[l]]>a[index[l+1]])
          {temp = index[l]; index[l]=index[l+1]; index[l+1]=temp;}

          i=l+1;        //initialize pointers for partitioning
          j=ir;
          arr = index[l+1];
          for (;;) {
             do i++; while (a[index[i]]<a[arr]);
             do j--; while (a[index[j]]>a[arr]);
             if (j<i) break;  //pointers crossed, partitioning complete
             {temp=index[i]; index[i]=index[j]; index[j]=temp;}
          }
          index[l+1]=index[j];
          index[j]=arr;
          if (j>=rk) ir = j-1; //keep active the partition that
          if (j<=rk) l=i;      //contains the k_th element
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Build boundaries for each node. Note, that the boundaries here are built
/// based on the splitting planes of the kd-tree, and don't necessarily pass
/// through the points of the original dataset. For the latter functionality
/// see function MakeBoundariesExact()
/// Boundaries can be retrieved by calling GetBoundary(inode) function that would
/// return an array of boundaries for the specified node, or GetBoundaries() function
/// that would return the complete array.

template <typename Index, typename Value>
void TKDTree<Index, Value>::MakeBoundaries(Value *range)
{

   if(range) memcpy(fRange, range, fNDimm*sizeof(Value));
   // total number of nodes including terminal nodes
   Int_t totNodes = fNNodes + fNPoints/fBucketSize + ((fNPoints%fBucketSize)?1:0);
   fBoundaries = new Value[totNodes*fNDimm];
   //Info("MakeBoundaries(Value*)", Form("Allocate boundaries for %d nodes", totNodes));


   // loop
   Value *tbounds = 0x0, *cbounds = 0x0;
   Int_t cn;
   for(int inode=fNNodes-1; inode>=0; inode--){
      tbounds = &fBoundaries[inode*fNDimm];
      memcpy(tbounds, fRange, fNDimm*sizeof(Value));

      // check left child node
      cn = (inode<<1)+1;
      if(IsTerminal(cn)) CookBoundaries(inode, kTRUE);
      cbounds = &fBoundaries[fNDimm*cn];
      for(int idim=0; idim<fNDim; idim++) tbounds[idim<<1] = cbounds[idim<<1];

      // check right child node
      cn = (inode+1)<<1;
      if(IsTerminal(cn)) CookBoundaries(inode, kFALSE);
      cbounds = &fBoundaries[fNDimm*cn];
      for(int idim=0; idim<fNDim; idim++) tbounds[(idim<<1)+1] = cbounds[(idim<<1)+1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// define index of this terminal node

template <typename Index, typename Value>
void TKDTree<Index, Value>::CookBoundaries(const Int_t node, Bool_t LEFT)
{
   Int_t index = (node<<1) + (LEFT ? 1 : 2);
   //Info("CookBoundaries()", Form("Node %d", index));

   // build and initialize boundaries for this node
   Value *tbounds = &fBoundaries[fNDimm*index];
   memcpy(tbounds, fRange, fNDimm*sizeof(Value));
   Bool_t flag[256];  // cope with up to 128 dimensions
   memset(flag, kFALSE, fNDimm);
   Int_t nvals = 0;

   // recurse parent nodes
   Int_t pn = node;
   while(pn >= 0 && nvals < fNDimm){
      if(LEFT){
         index = (fAxis[pn]<<1)+1;
         if(!flag[index]) {
            tbounds[index] = fValue[pn];
            flag[index] = kTRUE;
            nvals++;
         }
      } else {
         index = fAxis[pn]<<1;
         if(!flag[index]) {
            tbounds[index] = fValue[pn];
            flag[index] = kTRUE;
            nvals++;
         }
      }
      LEFT = pn&1;
      pn =  (pn - 1)>>1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Build boundaries for each node. Unlike MakeBoundaries() function
/// the boundaries built here always pass through a point of the original dataset
/// So, for example, for a terminal node with just one point minimum and maximum for each
/// dimension are the same.
/// Boundaries can be retrieved by calling GetBoundaryExact(inode) function that would
/// return an array of boundaries for the specified node, or GetBoundaries() function
/// that would return the complete array.

template <typename Index, typename Value>
void TKDTree<Index, Value>::MakeBoundariesExact()
{

   // total number of nodes including terminal nodes
   //Int_t totNodes = fNNodes + fNPoints/fBucketSize + ((fNPoints%fBucketSize)?1:0);
   if (fBoundaries){
      //boundaries were already computed for this tree
      return;
   }
   fBoundaries = new Value[fTotalNodes*fNDimm];
   Value *min = new Value[fNDim];
   Value *max = new Value[fNDim];
   for (Index inode=fNNodes; inode<fTotalNodes; inode++){
      //go through the terminal nodes
      for (Index idim=0; idim<fNDim; idim++){
         min[idim]= std::numeric_limits<Value>::max();
         max[idim]=-std::numeric_limits<Value>::max();
      }
      Index *points = GetPointsIndexes(inode);
      Index npoints = GetNPointsNode(inode);
      //find max and min in each dimension
      for (Index ipoint=0; ipoint<npoints; ipoint++){
         for (Index idim=0; idim<fNDim; idim++){
            if (fData[idim][points[ipoint]]<min[idim])
               min[idim]=fData[idim][points[ipoint]];
            if (fData[idim][points[ipoint]]>max[idim])
               max[idim]=fData[idim][points[ipoint]];
         }
      }
      for (Index idim=0; idim<fNDimm; idim+=2){
         fBoundaries[inode*fNDimm + idim]=min[idim/2];
         fBoundaries[inode*fNDimm + idim+1]=max[idim/2];
      }
   }

   delete [] min;
   delete [] max;

   Index left, right;
   for (Index inode=fNNodes-1; inode>=0; inode--){
      //take the min and max of left and right
      left = GetLeft(inode)*fNDimm;
      right = GetRight(inode)*fNDimm;
      for (Index idim=0; idim<fNDimm; idim+=2){
         //take the minimum on each dimension
         fBoundaries[inode*fNDimm+idim]=TMath::Min(fBoundaries[left+idim], fBoundaries[right+idim]);
         //take the maximum on each dimension
         fBoundaries[inode*fNDimm+idim+1]=TMath::Max(fBoundaries[left+idim+1], fBoundaries[right+idim+1]);

      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///
/// find the smallest node covering the full range - start
///

template <typename  Index, typename Value>
   void TKDTree<Index, Value>::FindBNodeA(Value *point, Value *delta, Int_t &inode){
   inode =0;
   for (;inode<fNNodes;){
      if (TMath::Abs(point[fAxis[inode]] - fValue[inode])<delta[fAxis[inode]]) break;
      inode = (point[fAxis[inode]] < fValue[inode]) ? (inode*2)+1: (inode*2)+2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the boundaries

template <typename  Index, typename Value>
   Value* TKDTree<Index, Value>::GetBoundaries()
{
   if(!fBoundaries) MakeBoundaries();
   return fBoundaries;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the boundaries

template <typename  Index, typename Value>
   Value* TKDTree<Index, Value>::GetBoundariesExact()
{
   if(!fBoundaries) MakeBoundariesExact();
   return fBoundaries;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a boundary

template <typename  Index, typename Value>
   Value* TKDTree<Index, Value>::GetBoundary(const Int_t node)
{
   if(!fBoundaries) MakeBoundaries();
   return &fBoundaries[node*2*fNDim];
}

////////////////////////////////////////////////////////////////////////////////
/// Get a boundary

template <typename  Index, typename Value>
Value* TKDTree<Index, Value>::GetBoundaryExact(const Int_t node)
{
   if(!fBoundaries) MakeBoundariesExact();
   return &fBoundaries[node*2*fNDim];
}


////////////////////////////////////////////////////////////////////////////////
///
/// Example function to
///

TKDTreeIF *TKDTreeTestBuild(const Int_t npoints, const Int_t bsize){
   Float_t *data0 =  new Float_t[npoints*2];
   Float_t *data[2];
   data[0] = &data0[0];
   data[1] = &data0[npoints];
   for (Int_t i=0;i<npoints;i++) {
      data[1][i]= gRandom->Rndm();
      data[0][i]= gRandom->Rndm();
   }
   TKDTree<Int_t, Float_t> *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
   return kdtree;
}



template class TKDTree<Int_t, Float_t>;
template class TKDTree<Int_t, Double_t>;
