// @(#)root/eve7:$Id$
// Author: Sergey Linev, 14.12.2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveGeomData
#define ROOT7_REveGeomData

#include <ROOT/REveRenderData.hxx>

#include <vector>
#include <string>
#include <functional>
#include <memory>

class TGeoNode;
class TGeoManager;
class TGeoShape;
class TGeoMatrix;
class TGeoVolume;

// do not use namespace to avoid too long JSON

namespace ROOT {
namespace Experimental {

class REveRenderData;
class RGeomBrowserIter;

/** Base description of geometry node, required only to build hierarchy */

class REveGeomNodeBase {
public:
   enum EVis { vis_off = 0, vis_this = 1, vis_chlds = 2, vis_lvl1 = 4 };

   int id{0};               ///< node id, index in array
   std::string name;        ///< node name
   std::vector<int> chlds;  ///< list of childs id
   int vis{vis_off};        ///< visibility flag, combination of EVis flags
   std::string color;       ///< rgb code without rgb() prefix
   int sortid{0};           ///<! place in sorted array, to check cuts, or id of original node when used search structures

   REveGeomNodeBase(int _id = 0) : id(_id) {}

   bool IsVisible() const { return vis & vis_this; }

   int GetVisDepth() const { return (vis & vis_chlds) ? 999999 : ((vis & vis_lvl1) ? 1 : 0); }

   /** Set indication if node really rendered - depends from selection */
   // void SetDisplayed(bool on) { vis = on ? (vis | vis_displayed) : (vis & ~vis_displayed); }
   // bool IsDisplayed() const { return vis & vis_displayed; }
};

/** Full node description including matrices and other attributes */

class REveGeomNode : public REveGeomNodeBase  {
public:
   std::vector<float> matr; ///< matrix for the node, can have reduced number of elements
   double vol{0};           ///<! volume estimation
   int nfaces{0};           ///<! number of shape faces
   int numvischld{0};       ///<! number of visible childs, if all can be jump over
   int idshift{0};          ///<! used to jump over then scan all geom hierarchy
   bool useflag{false};     ///<! extra flag, used for selection
   float opacity{1.};       ///<! opacity of the color

   REveGeomNode(int _id = 0) : REveGeomNodeBase(_id) {}

   /** True when there is shape and it can be displayed */
   bool CanDisplay() const { return (vol > 0.) && (nfaces > 0); }
};


/** Information block for render data, stored in binary buffer */

class REveShapeRenderInfo {
public:
   // render data, equivalent of REveElement::WriteCoreJson
   int rnr_offset{-1};        ///< rnr_offset;
   std::string rnr_func;      ///< fRenderData->GetRnrFunc();
   int vert_size{0};          ///< fRenderData->SizeV();
   int norm_size{0};          ///< fRenderData->SizeN();
   int index_size{0};         ///< fRenderData->SizeI();
   TGeoShape *shape{nullptr}; ///< original shape - can be much less than binary data
   std::vector<unsigned char> raw;  ///< raw shape data with render information, JSON_base64
   // int trans_size{0};      ///< fRenderData->SizeT(); not used in GeomViewer
};

/** REveGeomVisible contains description of visible node
 * It is path to the node plus reference to shape rendering data */

class REveGeomVisible {
public:
   int nodeid{0};                    ///< selected node id,
   std::vector<int> stack;           ///< path to the node, index in list of childs
   std::string color;                ///< color in rgb format
   double opacity{1};                ///< opacity
   REveShapeRenderInfo *ri{nullptr}; ///< render information for the shape, can be same for different nodes

   REveGeomVisible() = default;
   REveGeomVisible(int id, const std::vector<int> &_stack) : nodeid(id), stack(_stack) {}
};

/** Object with full description for drawing geometry
 * It includes list of visible items and list of nodes required to build them */

class REveGeomDrawing {
public:
   int numnodes{0};                         ///< total number of nodes in description
   std::string drawopt;                     ///< draw options for TGeoPainter
   int nsegm{0};                            ///< number of segments for cylindrical shapes
   int binlen{0};                           ///< extra binary data for that drawing
   std::vector<REveGeomNode*> nodes;        ///< all used nodes to display visible items and not known for client
   std::vector<REveGeomVisible> visibles;   ///< all visible items
};


/** Request object send from client for different operations */
class REveGeomRequest {
public:
   std::string oper;  ///< operation like HIGHL or HOVER
   std::string path;  ///< path parameter, used with HOVER
   std::vector<int> stack; ///< stack parameter, used with HIGHL
};

class REveGeomNodeInfo {
public:
   std::string fullpath;  ///< full path to node
   std::string node_type;  ///< node class name
   std::string node_name;  ///< node name
   std::string shape_type; ///< shape type (if any)
   std::string shape_name; ///< shape class name (if any)

   REveShapeRenderInfo *ri{nullptr}; ///< rendering information (if applicable)

   std::vector<unsigned char> rndr_binary; ///<  binary render data (if available)
};

using REveGeomScanFunc_t = std::function<bool(REveGeomNode &, std::vector<int> &, bool)>;


class REveGeomDescription {

   friend class RGeomBrowserIter;

   class ShapeDescr {
   public:
      int id{0};                                   ///<! sequential id
      TGeoShape *fShape{nullptr};                  ///<! original shape
      int nfaces{0};                               ///<! number of faces in render data
      std::unique_ptr<REveRenderData> fRenderData; ///<! binary render data
      REveShapeRenderInfo fRenderInfo;             ///<! render information for client
      ShapeDescr(TGeoShape *s) : fShape(s) {}

      /// Provide render info for visible item
      REveShapeRenderInfo *rndr_info() { return (nfaces>0) && (fRenderInfo.rnr_offset>=0) ? &fRenderInfo : nullptr; }
   };

   std::vector<TGeoNode *> fNodes;  ///<! flat list of all nodes
   std::string fDrawOptions;        ///< default draw options for client
   std::vector<REveGeomNode> fDesc; ///< converted description, send to client

   int fTopDrawNode{0};             ///<! selected top node
   std::vector<int> fSortMap;       ///<! nodes in order large -> smaller volume
   int fNSegments{0};               ///<! number of segments for cylindrical shapes
   std::vector<ShapeDescr> fShapes; ///<! shapes with created descriptions
   std::vector<REveRenderData*> fRndrShapes; ///<! list of shapes which should be packet into binary
   int fRndrOffest{0};              ///<! current render offset

   std::string fDrawJson;           ///<! JSON with main nodes drawn by client
   std::vector<unsigned char> fDrawBinary;   ///<! binary data for main draw nodes
   int fDrawIdCut{0};               ///<! sortid used for selection of most-significant nodes
   int fFacesLimit{0};              ///<! maximal number of faces to be selected for drawing
   int fNodesLimit{0};              ///<! maximal number of nodes to be selected for drawing
   bool fPreferredOffline{false};   ///<! indicates that full description should be provided to client
   bool fBuildShapes{true};         ///<! if TGeoShape build already on server (default) or send as is to client

   int fJsonComp{0};                ///<! default JSON compression

   void PackMatrix(std::vector<float> &arr, TGeoMatrix *matr);

   void ScanNode(TGeoNode *node, std::vector<int> &numbers, int offset);

   int MarkVisible(bool on_screen = false);

   void ScanNodes(bool only_visible, REveGeomScanFunc_t func);

   void ResetRndrInfos();

   ShapeDescr &FindShapeDescr(TGeoShape *shape);

   ShapeDescr &MakeShapeDescr(TGeoShape *shape, bool acc_rndr = false);

   void BuildRndrBinary(std::vector<unsigned char> &buf);

   void CopyMaterialProperties(TGeoVolume *vol, REveGeomNode &node);

   void CollectNodes(REveGeomDrawing &drawing);

public:
   REveGeomDescription() = default;

   void Build(TGeoManager *mgr);

   /** Number of unique nodes in the geometry */
   int GetNumNodes() const { return fDesc.size(); }

   bool IsBuild() const { return GetNumNodes() > 0; }

   /** Set maximal number of nodes which should be selected for drawing */
   void SetMaxVisNodes(int cnt) { fNodesLimit = cnt; }

   /** Returns maximal visible number of nodes, ignored when non-positive */
   int GetMaxVisNodes() const { return fNodesLimit; }

   /** Set maximal number of faces which should be selected for drawing */
   void SetMaxVisFaces(int cnt) { fFacesLimit = cnt; }

   /** Returns maximal visible number of faces, ignored when non-positive */
   int GetMaxVisFaces() const { return fFacesLimit; }

   /** Set preference of offline operations.
    * Server provides more info to client from the begin on to avoid communication */
   void SetPreferredOffline(bool on) { fPreferredOffline = on; }

   /** Is offline operations preferred.
    * After get full description, client can do most operations without extra requests */
   bool IsPreferredOffline() const { return fPreferredOffline; }

   bool CollectVisibles();

   bool IsPrincipalEndNode(int nodeid);

   std::string ProcessBrowserRequest(const std::string &req = "");

   bool HasDrawData() const { return (fDrawJson.length() > 0) && (fDrawBinary.size() > 0) && (fDrawIdCut > 0); }
   const std::string &GetDrawJson() const { return fDrawJson; }
   const std::vector<unsigned char> &GetDrawBinary() const { return fDrawBinary; }
   void ClearRawData();

   int SearchVisibles(const std::string &find, std::string &hjson, std::string &json, std::vector<unsigned char> &binary);

   int FindNodeId(const std::vector<int> &stack);

   std::string ProduceModifyReply(int nodeid);

   std::vector<int> MakeStackByIds(const std::vector<int> &ids);

   std::vector<int> MakeIdsByStack(const std::vector<int> &stack);

   std::vector<int> MakeStackByPath(const std::string &path);

   std::string MakePathByStack(const std::vector<int> &stack);

   bool ProduceDrawingFor(int nodeid, std::string &json, std::vector<unsigned char> &binary, bool check_volume = false);

   bool ChangeNodeVisibility(int nodeid, bool selected);

   void SelectVolume(TGeoVolume *);

   void SelectNode(TGeoNode *);

   /** Set number of segments for cylindrical shapes, if 0 - default value will be used */
   void SetNSegments(int n = 0) { fNSegments = n; }
   /** Return of segments for cylindrical shapes, if 0 - default value will be used */
   int GetNSegments() const { return fNSegments; }

   /** Set JSON compression level for data transfer */
   void SetJsonComp(int comp = 0) { fJsonComp = comp; }
   /** Returns JSON compression level for data transfer */
   int GetJsonComp() const  { return fJsonComp; }

   /** Set draw options as string for JSROOT TGeoPainter */
   void SetDrawOptions(const std::string &opt = "") { fDrawOptions = opt; }
   /** Returns draw options, used for JSROOT TGeoPainter */
   std::string GetDrawOptions() const { return fDrawOptions; }

   /** Instruct to build binary 3D model already on the server (true) or send TGeoShape as is to client, which can build model itself */
   void SetBuildShapes(bool on = true) { fBuildShapes = on; }
   /** Retuns true if binary 3D model build already by C++ server (default) */
   bool IsBuildShapes() const { return fBuildShapes; }

   std::unique_ptr<REveGeomNodeInfo> MakeNodeInfo(const std::string &path);
};


} // namespace Experimental
} // namespace ROOT

#endif
