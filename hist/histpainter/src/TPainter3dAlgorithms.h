// @(#)root/histpainter:$Id$
// Author: Rene Brun, Evgueni Tcherniaev, Olivier Couet   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPainter3dAlgorithms
#define ROOT_TPainter3dAlgorithms

#include "TAttLine.h"

#include "TAttFill.h"

const Int_t kCARTESIAN   = 1;
const Int_t kPOLAR       = 2;
const Int_t kCYLINDRICAL = 3;
const Int_t kSPHERICAL   = 4;
const Int_t kRAPIDITY    = 5;

class TF3;
class TView;

class TPainter3dAlgorithms : public TAttLine, public TAttFill {

private:
   Double_t     fRmin[3];          /// Lower limits of lego
   Double_t     fRmax[3];          /// Upper limits of lego
   Double_t     *fAphi;            ///
   Int_t        fNaphi;            /// Size of fAphi
   Int_t        fSystem;           /// Coordinate system
   Int_t       *fColorMain;        ///
   Int_t       *fColorDark;        ///
   Int_t        fColorTop;         ///
   Int_t        fColorBottom;      ///
   Int_t       *fEdgeColor;        ///
   Int_t       *fEdgeStyle;        ///
   Int_t       *fEdgeWidth;        ///
   Int_t        fEdgeIdx;          ///
   Int_t        fMesh;             /// (=1 if mesh to draw, o otherwise)
   Int_t        fNStack;           /// Number of histograms in the stack to be painted
   Double_t     fFmin;             /// IsoSurface minimum function value
   Double_t     fFmax;             /// IsoSurface maximum function value
   Int_t        fNcolor;           /// Number of colours per Iso surface
   Int_t        fIc1;              /// Base colour for the 1st Iso Surface
   Int_t        fIc2;              /// Base colour for the 2nd Iso Surface
   Int_t        fIc3;              /// Base colour for the 3rd Iso Surface

public:
   typedef void (TPainter3dAlgorithms::*DrawFaceFunc_t)(Int_t *, Double_t *, Int_t, Int_t *, Double_t *);
   typedef void (TPainter3dAlgorithms::*LegoFunc_t)(Int_t,Int_t,Int_t&,Double_t*,Double_t*,Double_t*);
   typedef void (TPainter3dAlgorithms::*SurfaceFunc_t)(Int_t,Int_t,Double_t*,Double_t*);

private:
   DrawFaceFunc_t  fDrawFace;        /// pointer to face drawing function
   LegoFunc_t      fLegoFunction;    /// pointer to lego function
   SurfaceFunc_t   fSurfaceFunction; /// pointer to surface function

public:
   TPainter3dAlgorithms();
   TPainter3dAlgorithms(Double_t *rmin, Double_t *rmax, Int_t system=1);
   virtual ~TPainter3dAlgorithms();
   void    BackBox(Double_t ang);
   void    FrontBox(Double_t ang);
   void    DrawFaceGouraudShaded(Int_t *icodes, Double_t xyz[][3], Int_t np, Int_t *iface, Double_t *t);
   void    DrawFaceMode1(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *t);
   void    DrawFaceMode2(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *t);
   void    DrawFaceMode3(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *t);
   void    DrawFaceMove1(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *tt);
   void    DrawFaceMove2(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *tt);
   void    DrawFaceMove3(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *tt);
   void    DrawLevelLines(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *tt);
   void    DrawFaceRaster1(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *tt);
   void    DrawFaceRaster2(Int_t *icodes, Double_t *xyz, Int_t np, Int_t *iface, Double_t *tt);
   void    GouraudFunction(Int_t ia, Int_t ib, Double_t *f, Double_t *t);
   void    ImplicitFunction(TF3 *f3, Double_t *rmin, Double_t *rmax, Int_t nx, Int_t ny, Int_t nz, const char *chopt);
   void    IsoSurface (Int_t ns, Double_t *s, Int_t nx, Int_t ny, Int_t nz, Double_t *x, Double_t *y, Double_t *z, const char *chopt);
   void    LegoCartesian(Double_t ang, Int_t nx, Int_t ny, const char *chopt);
   void    LegoFunction(Int_t ia, Int_t ib, Int_t &nv, Double_t *ab, Double_t *vv, Double_t *t);
   void    LegoPolar(Int_t iordr, Int_t na, Int_t nb, const char *chopt);
   void    LegoCylindrical(Int_t iordr, Int_t na, Int_t nb, const char *chopt);
   void    LegoSpherical(Int_t ipsdr, Int_t iordr, Int_t na, Int_t nb, const char *chopt);
   void    SetDrawFace(DrawFaceFunc_t pointer);
   void    SetIsoSurfaceParameters(Double_t fmin, Double_t fmax, Int_t ncolor, Int_t ic1, Int_t ic2, Int_t ic3){fFmin=fmin; fFmax=fmax; fNcolor=ncolor; fIc1=ic1; fIc2=ic2; fIc3=ic3;}
   void    SetLegoFunction(LegoFunc_t pointer);
   void    SetMesh(Int_t mesh=1) {fMesh=mesh;}
   void    SetSurfaceFunction(SurfaceFunc_t pointer);
   void    SetColorDark(Color_t color, Int_t n=0);
   void    SetColorMain(Color_t color, Int_t n=0);
   void    SetEdgeAtt(Color_t color=1, Style_t style=1, Width_t width=1, Int_t n=0);
   void    SideVisibilityDecode(Double_t val, Int_t &iv1, Int_t &iv2, Int_t &iv3, Int_t &iv4, Int_t &iv5, Int_t &iv6, Int_t &ir);
   void    SideVisibilityEncode(Int_t iopt, Double_t phi1, Double_t phi2, Double_t &val);
   void    SurfaceCartesian(Double_t ang, Int_t nx, Int_t ny, const char *chopt);
   void    SurfacePolar(Int_t iordr, Int_t na, Int_t nb, const char *chopt);
   void    SurfaceCylindrical(Int_t iordr, Int_t na, Int_t nb, const char *chopt);
   void    SurfaceFunction(Int_t ia, Int_t ib, Double_t *f, Double_t *t);
   void    SurfaceSpherical(Int_t ipsdr, Int_t iordr, Int_t na, Int_t nb, const char *chopt);

//       Color and function levels
//
public:
   void    DefineGridLevels(Int_t ndivz);
   void    ColorFunction(Int_t nl, Double_t *fl, Int_t *icl, Int_t &irep);
   void    Spectrum(Int_t nl, Double_t fmin, Double_t fmax, Int_t ic, Int_t idc, Int_t &irep);
   void    FindLevelLines(Int_t np, Double_t *f, Double_t *t);
   void    FindPartEdge(Double_t *p1, Double_t *p2, Double_t f1, Double_t f2, Double_t fmin, Double_t fmax, Int_t &kpp, Double_t *pp);
   void    FillPolygon(Int_t n, Double_t *p, Double_t *f);

private:
   static const Int_t    NumOfColorLevels = 256;
   Int_t        fNlevel;                         // Number of color levels
   Double_t     fFunLevel[NumOfColorLevels+1];   // Function levels corresponding to colors
   Int_t        fColorLevel[NumOfColorLevels+2]; // Color levels corresponding to functions

   static const Int_t    NumOfLevelLines = 200;
   Int_t        fNlines;                         // Number of lines
   Int_t        fLevelLine[NumOfLevelLines];     // Corresponding levels
   Double_t     fPlines[NumOfLevelLines*6];      // End points of lines

   void    Luminosity(TView *view, Double_t *anorm, Double_t &flum);

//       Light and surface properties
//
public:
   void    LightSource(Int_t nl, Double_t yl, Double_t xscr, Double_t yscr, Double_t zscr, Int_t &irep);
   void    SurfaceProperty(Double_t qqa, Double_t qqd, Double_t qqs, Int_t nnqs, Int_t &irep);

private:
   static const Int_t    NumOfLights = 4;
   Int_t        fLoff;
   Double_t     fYdl;
   Double_t     fYls[NumOfLights];
   Double_t     fVls[NumOfLights*3];
   Double_t     fQA;
   Double_t     fQD;
   Double_t     fQS;
   Int_t        fNqs;

//        Moving screen - specialized hidden line removal algorithm
//        for rendering 2D histograms
//
public:
   void    InitMoveScreen(Double_t xmin, Double_t xmax);
   void    FindVisibleDraw(Double_t *r1, Double_t *r2);
   void    ModifyScreen(Double_t *r1, Double_t *r2);

private:
   static const Int_t    MaxNT = 100;
   Double_t     fX0;         // minimal x
   Double_t     fDX;         // x size
   Int_t        fNT;         // number of t-segments
   Double_t     fT[MaxNT*2]; // t-segments

   static const Int_t    NumOfSlices = 2000;
   Double_t    fU[NumOfSlices*2];
   Double_t    fD[NumOfSlices*2];

//        Raster screen - specialized hidden line removal algorithm
//        for rendering 3D polygons ordered by depth (front-to-back)
//
public:
   void    InitRaster(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Int_t nx, Int_t ny);
   void    ClearRaster();
   void    FindVisibleLine(Double_t *p1, Double_t *p2, Int_t ntmax, Int_t &nt, Double_t *t);
   void    FillPolygonBorder(Int_t nn, Double_t *xy);

private:
   Double_t    fXrast;     // minimal x
   Double_t    fYrast;     // minimal y
   Double_t    fDXrast;    // x size
   Double_t    fDYrast;    // y size
   Int_t       fNxrast;    // number of pixels in x
   Int_t       fNyrast;    // number of pixels in y
   Int_t       fIfrast;    // flag, if it not zero them the algorithm is off
   Int_t       *fRaster;   // pointer to raster buffer
   Int_t       fJmask[30]; // indices of subsets of n-bit masks (n is from 1 to 30)
   Int_t       fMask[465]; // set of masks (30+29+28+...+1)=465

//        Marching Cubes 33 - construction of iso-surfaces, see publication CERN-CN-95-17
//
public:
   void    MarchingCube(Double_t fiso, Double_t p[8][3], Double_t f[8], Double_t g[8][3], Int_t &nnod, Int_t &ntria, Double_t xyz[][3], Double_t grad[][3], Int_t itria[][3]);

protected:
   void    MarchingCubeCase00(Int_t k1, Int_t k2, Int_t k3, Int_t k4, Int_t k5, Int_t k6, Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase03(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase04(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase06(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase07(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase10(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase12(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeCase13(Int_t &nnod, Int_t &ntria, Double_t xyz[52][3], Double_t grad[52][3], Int_t itria[48][3]);
   void    MarchingCubeSetTriangles(Int_t ntria, Int_t it[][3], Int_t itria[48][3]);
   void    MarchingCubeMiddlePoint(Int_t nnod, Double_t xyz[52][3], Double_t grad[52][3], Int_t it[][3], Double_t *pxyz, Double_t *pgrad);
   void    MarchingCubeSurfacePenetration(Double_t a00, Double_t a10, Double_t a11, Double_t a01, Double_t b00, Double_t b10, Double_t b11, Double_t b01, Int_t &irep);
   void    MarchingCubeFindNodes(Int_t nnod, Int_t *ie, Double_t xyz[52][3], Double_t grad[52][3]);

private:
   Double_t    fP8[8][3]; // vertices
   Double_t    fF8[8];    // function values
   Double_t    fG8[8][3]; // function gradients

//        Z-depth sorting algorithm for set of triangles
//
public:
   void    ZDepth(Double_t xyz[52][3], Int_t &nface, Int_t iface[48][3], Double_t dface[48][6], Double_t abcd[48][4], Int_t *iorder);

protected:
   void    TestEdge(Double_t del, Double_t xyz[52][3], Int_t i1, Int_t i2, Int_t iface[3], Double_t abcd[4], Int_t &irep);

};

#endif
