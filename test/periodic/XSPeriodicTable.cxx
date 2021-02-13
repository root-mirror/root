/*
 * $Header$
 * $Log$
 *
 * Implements the periodic table of elements
 */

#include <TGFrame.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <snprintf.h>

#include "XSVarious.h"
#include "XSElementList.h"
#include "XSPeriodicTable.h"

//ClassImp(XSTblElement);

/* =========== XSTblElement ============== */
XSTblElement::XSTblElement(const TGWindow *p, Int_t z, UInt_t color)
   : TGButton(p,z)
{
   Z = z;

   char str[5];
   snprintf(str,5,"%d",Z);
   lZ    = new TGLabel(this,str);
   lName = new TGLabel(this,XSelements->Mnemonic(Z), blueBoldGC);

   ChangeBackground(color);

   tpZ = NULL;
   tpName = NULL;
} // XSTblElement

/* ----- ~XSTblElement ----- */
XSTblElement::~XSTblElement()
{
   delete lZ;
   delete lName;
   delete tpZ;
   delete tpName;
} // ~XSTblElement

/* ----- Layout ----- */
void
XSTblElement::Layout()
{
   int   w = GetWidth() - 2*fBorderWidth;
   int   h = (GetHeight() - 2*fBorderWidth) / 2;

   lZ->MoveResize(fBorderWidth, fBorderWidth, w, h);
   lZ->MapWindow();

   lName->MoveResize(fBorderWidth, fBorderWidth+h, w, h);
   lName->MapWindow();

   if (tpZ == NULL) {
      tpZ = new TGToolTip(fClient->GetRoot(),
                          lZ, XSelements->Name(Z), 1000);
      tpName = new TGToolTip(fClient->GetRoot(),
                             lName, XSelements->Name(Z), 1000);
   }

   TGButton::Layout();
} // Layout

/* ------ SetState ------ */
void
XSTblElement::SetState(EButtonState state, Bool_t emit)
{
   if (state != fState) {
      if (state==kButtonDown /*|| state==kButtonEngaged*/) {
         lZ->Move(lZ->GetX()+1, lZ->GetY()+1);
         lName->Move(lName->GetX()+1, lName->GetY()+1);
      } else
      if (state == kButtonUp) {
         lZ->Move(lZ->GetX()-1, lZ->GetY()-1);
         lName->Move(lName->GetX()-1, lName->GetY()-1);
      }
   }

   TGButton::SetState(state, emit);
} // SetState

/* ------ SetState ------ */
void
XSTblElement::ChangeBackground( ULong_t color )
{
   lZ->ChangeBackground(color);
   lName->ChangeBackground(color);
} // ChangeBackground

/////////////////////////////////////////////////////////////////////

//ClassImp(XSPeriodicTable);

static   Int_t   colwidth[XSPTBL_COLS] =
      { 8, 5, 5, 2, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5 };
static   Int_t   colsum;
static   Int_t   rowheight[XSPTBL_ROWS] =
      { 5, 5, 10, 10, 10, 10, 10, 10, 10, 3, 10, 10 };
static   Int_t   rowsum;

#define   A(x)   (Long_t)(char *)(x)
static   Long_t   ptable[XSPTBL_ROWS][XSPTBL_COLS] = {
   { A("Group"), A("1"), A("2"), 0, A("3"), A("4"), A("5"), A("6"), A("7"),
   A("8"), A("9"), A("10"), A("11"), A("12"), A("13"), A("14"), A("15"),
   A("16"), A("17"), A("18") },
   { A("Period"), 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0},
   { A("1"), 1, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 2},
   { A("2"), 3, 4, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 5,  6, 7, 8, 9,10},
   { A("3"),11,12, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0,13, 14,15,16,17,18},
   { A("4"),19,20, 0,21, 22,23,24,25,26, 27,28,29,30,31, 32,33,34,35,36},
   { A("5"),37,38, 0,39, 40,41,42,43,44, 45,46,47,48,49, 50,51,52,53,54},
   { A("6"),55,56,A("*"),71, 72,73,74,75,76, 77,78,79,80,81, 82,83,84,85,86},
   { A("7"),87,88,A("**"),103, 104,105,106,107,108, 109,110,111,112,113, 114,115,116,117,118},
   { 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0 },
   { A("* Lanthanoids"),0,0,A("*"),57, 58,59,60,61,62, 63,64,65,66,67, 68,69,70,0,0},
   { A("** Actinoids"),0,0,A("**"),89, 90,91,92,93,94, 95,96,97,98,99, 100,101,102,0,0}
};
#undef A


#define COLORS   7

//#define w   0xFFFFFFL   // White
//#define A   0xCCCCCCL   // Gray
//#define B   0xCCCCFFL   // Blue
//#define Y   0xFFFFCCL   // Yellow
//#define O   0xFFCCCCL   // Orange
//#define G   0xCCFFCCL   // Green

#define w   0
#define   A   1
#define B   2
#define Y   3
#define   O   4
#define G   5
static const char *colorName[COLORS] = {
         "White",
         "#CCCCCC",
         "#CCCCFF",
         "#FFFFCC",
         "#FFCCCC",
         "#CCFFCC",
         "yellow",
/******
         "White",
         "LightGrey",
         "LightSteelBlue1",
         "lemon chiffon",
         "pink",
         "DarkSeaGreen1",
         "yellow",
*******/
      };

static   ULong_t   colorPixels[COLORS];

static Byte_t   ecolor[XSPTBL_ROWS][XSPTBL_COLS] = {
  { A,A,A,w,A, A,A,A,A,A, A,A,A,A,A, A,A,A,A,A },
  { A,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w },
  { A,B,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,Y },
  { A,B,B,w,w, w,w,w,w,w, w,w,w,w,Y, Y,Y,Y,Y,Y },
  { A,B,B,w,w, w,w,w,w,w, w,w,w,w,Y, Y,Y,Y,Y,Y },
  { A,B,B,w,O, O,O,O,O,O, O,O,O,O,Y, Y,Y,Y,Y,Y },
  { A,B,B,w,O, O,O,O,O,O, O,O,O,O,Y, Y,Y,Y,Y,Y },
  { A,B,B,w,O, O,O,O,O,O, O,O,O,O,Y, Y,Y,Y,Y,Y },
  { A,B,B,w,O, O,O,O,O,O, O,O,O,O,A, Y,A,Y,A,Y },
  { w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w },
  { A,A,A,w,G, G,G,G,G,G, G,G,G,G,G, G,G,G,w,w },
  { A,A,A,w,G, G,G,G,G,G, G,G,G,G,G, G,G,G,w,w }
};
#undef w
#undef A
#undef B
#undef Y
#undef O
#undef G

/* ----- XSPeriodicTable ----- */
XSPeriodicTable::XSPeriodicTable(const TGWindow *msgWnd, const TGWindow* p,
         UInt_t w, UInt_t h)
   : TGCompositeFrame(p,w,h, kFitWidth | kFitHeight | kSunkenFrame)
{
   int      i, j;

   rowsum = 0;
   for (j=0; j<XSPTBL_ROWS; j++)
      rowsum += rowheight[j];

   colsum = 0;
   for (j=0; j<XSPTBL_COLS; j++)
      colsum += colwidth[j];

   /** Initialise colors **/
   for (i=0; i<COLORS; i++)
      gClient->GetColorByName(colorName[i], colorPixels[i]);

   for (j=0; j<XSPTBL_ROWS; j++)
      for (i=0; i<XSPTBL_COLS; i++) {
         long   val = ptable[j][i];
         int   col = colorPixels[ecolor[j][i]];
         if (val) {
            if (val < 120) {
               elem[j][i] = new XSTblElement(this,
                     val,col);
               ((XSTblElement*)elem[j][i])->Associate(msgWnd);
            } else {
               elem[j][i] = new TGLabel(this,
                     (char*)val);
               elem[j][i]->ChangeBackground(col);
            }
         } else
            elem[j][i] = NULL;
      }
   ChangeBackground(colorPixels[0]);

   width = w;
   height = h;
} // XSPeriodicTable

/* ----- ~XSPeriodicTable ----- */
XSPeriodicTable::~XSPeriodicTable()
{
   for (int j=0; j<XSPTBL_ROWS; j++)
      for (int i=0; i<XSPTBL_COLS; i++)
         if (elem[j][i])
            delete elem[j][i];
} // ~XSPeriodicTable

/* ----- Layout ----- */
void
XSPeriodicTable::Layout()
{
   int   w, h;
   Float_t   cwidth[XSPTBL_COLS];
   Float_t   rheight[XSPTBL_ROWS];

   TGCompositeFrame::Layout();

   w = GetWidth() - 2*fBorderWidth - 4;
   h = GetHeight() - 2*fBorderWidth - 4;

   // --- Normalise height and widths ---
   for (int i=0; i<XSPTBL_COLS; i++) {
      cwidth[i] = (Float_t)(colwidth[i]*w)/(Float_t)colsum;
   }

   for (int i=0; i<XSPTBL_ROWS; i++) {
      rheight[i] = (Float_t)(rowheight[i]*h)/(Float_t)rowsum;
   }

   Float_t   yt = fBorderWidth+2;
   for (int j=0; j<XSPTBL_ROWS; j++) {
      Float_t xl = fBorderWidth+2;
      for (int i=0; i<XSPTBL_COLS; i++) {
         if (j<10) {
            if (elem[j][i]) {
               elem[j][i]->MoveResize(
                  (int)xl,
                  (int)yt,
                  (int)cwidth[i]-1,
                  (int)rheight[j]-1);
               elem[j][i]->MapWindow();
            }
            xl += cwidth[i];
         } else {
            if (i<3) {
               Float_t w = cwidth[0]+cwidth[1]+
                     cwidth[2];
               if (elem[j][0]) {
                  elem[j][0]->MoveResize(
                     (int)xl,
                     (int)yt,
                     (int)w-1,
                     (int)rheight[j]-1);
                  elem[j][0]->MapWindow();
               }
               xl += w;
               i = 2;
            } else {
               if (elem[j][i]) {
                  elem[j][i]->MoveResize(
                     (int)xl,
                     (int)yt,
                     (int)cwidth[i]-1,
                     (int)rheight[j]-1);
                  elem[j][i]->MapWindow();
               }
               xl += cwidth[i];
            }
         }
      }
      yt += rheight[j];
   }
} // Layout

/* ----- SelectZ ----- */
void
XSPeriodicTable::SelectZ( ULong_t Z )
{
   // Find the selection and change the color
   for (int j=0; j<XSPTBL_ROWS; j++)
      for (int i=0; i<XSPTBL_COLS; i++)
         if ((ULong_t)ptable[j][i] == Z) {
            // Setting the background
            elem[j][i]->ChangeBackground(colorPixels[6]);
         }
} // SelectZ
