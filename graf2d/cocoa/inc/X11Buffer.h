// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov   29/02/2012

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_X11Buffer
#define ROOT_X11Buffer

#include <vector>
#include <string>

#include <Cocoa/Cocoa.h>

#include "CocoaGuiTypes.h"
#include "GuiTypes.h"

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// Unfortunately, TGCocoa's drawing methods can be called in a                  //
// "wrong" time and place: not from QuartzView -drawRect.                       //
// For example, on mouse move. This is bad and unnatural for Cocoa application, //
// since I expect GUI to draw only when I'm ready == ... called from drawRect.  //
// In X11 commands are buffered and this buffer is flushed at some points.      //
// I'm trying to emulate this, just to make GUI happy.                          //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

@class QuartzView;

namespace ROOT {
namespace MacOSX {

namespace Details {
class CocoaPrivate;
}

namespace X11 {

class Command {
   friend class CommandBuffer;

protected:
   const Drawable_t fID;
   const GCValues_t fGC;

public:
   Command(Drawable_t wid);
   Command(Drawable_t wid, const GCValues_t &gc);
   virtual ~Command();

   virtual bool HasOperand(Drawable_t drawable)const;
   virtual bool IsGraphicsCommand()const;//By-default - false.

   virtual void Execute()const = 0;
   virtual void Execute(CGContextRef /*ctx*/)const;

private:
   Command(const Command &rhs);
   Command &operator = (const Command &rhs);
};

class DrawLine : public Command {
private:
   const Point fP1;
   const Point fP2;

public:
   DrawLine(Drawable_t wid, const GCValues_t &gc, const Point &p1, const Point &p2);
   void Execute()const;
   bool IsGraphicsCommand()const
   {
      return true;
   }
};

class DrawSegments : public Command {
private:
   std::vector<Segment_t> fSegments;

public:
   DrawSegments(Drawable_t wid, const GCValues_t &gc, const Segment_t *segments, Int_t nSegments);
   void Execute()const;
   bool IsGraphicsCommand()const
   {
      return true;
   }
};

class ClearArea : public Command {
private:
   const Rectangle_t fArea;//to be replaced with X11::Rectangle

public:
   ClearArea(Window_t wid, const Rectangle_t &area);
   void Execute()const;
   bool IsGraphicsCommand()const
   {
      return true;
   }
};

class CopyArea : public Command {
private:
   const Drawable_t  fSrc;
   const Rectangle_t fArea;//to be replaced with X11::Rectangle
   const Point     fDstPoint;

public:
   CopyArea(Drawable_t src, Drawable_t dst, const GCValues_t &gc, const Rectangle_t &area, const Point &dstPoint);

   bool HasOperand(Drawable_t drawable)const;
   bool IsGraphicsCommand()const
   {
      return true;
   }

   void Execute()const;

};

class DrawString : public Command {
private:
   const Point       fPoint;
   const std::string fText;

public:
   DrawString(Drawable_t wid, const GCValues_t &gc, const Point &point, const std::string &text);

   bool IsGraphicsCommand()const
   {
      return true;
   }

   void Execute()const;
};

class FillRectangle : public Command {
private:
   const Rectangle_t fRectangle;//to be replaced with X11::Rectangle

public:
   FillRectangle(Drawable_t wid, const GCValues_t &gc, const Rectangle_t &rectangle);

   bool IsGraphicsCommand()const
   {
      return true;
   }

   void Execute()const;
};

class FillPolygon : public Command {
private:
   std::vector<Point_t> fPolygon;

public:
   FillPolygon(Drawable_t wid, const GCValues_t &gc, const Point_t *points, Int_t nPoints);

   bool IsGraphicsCommand()const
   {
      return true;
   }

   void Execute()const;
};

class DrawRectangle : public Command {
private:
   Rectangle_t fRectangle;//to be replaced with X11::Rectangle

public:
   DrawRectangle(Drawable_t wid, const GCValues_t &gc, const Rectangle_t &rectangle);

   bool IsGraphicsCommand()const
   {
      return true;
   }

   void Execute()const;
};

class UpdateWindow : public Command {
private:
   QuartzView *fView;

public:
   UpdateWindow(QuartzView *view);

   bool IsGraphicsCommand()const
   {
      return true;
   }

   void Execute()const;
};

class DeletePixmap : public Command {
public:
   DeletePixmap(Pixmap_t pixmap);
   void Execute()const;
};

//Set of 'xor' operations, required by TCanvas and ExecuteEvent's machinery.
class DrawBoxXor : public Command {
private:
   Point fP1;
   Point fP2;

public:
   DrawBoxXor(Window_t windowID, const Point &p1, const Point &p2);

   void Execute()const;
   void Execute(CGContextRef ctx)const;
};

class DrawLineXor : public Command {
private:
   Point fP1;
   Point fP2;

public:
   DrawLineXor(Window_t windowID, const Point &p1, const Point &p2);

   void Execute()const;
   void Execute(CGContextRef ctx)const;

   Point start() const {return fP1;}
   Point end() const {return fP2;}
};

class CommandBuffer {
private:
   CommandBuffer(const CommandBuffer &rhs);
   CommandBuffer &operator = (const CommandBuffer &rhs);

   std::vector<Command *> fCommands;
   std::vector<QuartzView *> fViewBranch;

   std::vector<Command *> fXorOps;
public:
   typedef std::vector<Command *>::size_type size_type;

   CommandBuffer();
   ~CommandBuffer();

   void AddDrawLine(Drawable_t wid, const GCValues_t &gc, Int_t x1, Int_t y1, Int_t x2, Int_t y2);
   void AddDrawSegments(Drawable_t wid, const GCValues_t &gc, const Segment_t *segments, Int_t nSegments);
   void AddClearArea(Window_t wid, Int_t x, Int_t y, UInt_t w, UInt_t h);
   void AddCopyArea(Drawable_t src, Drawable_t dst, const GCValues_t &gc,  Int_t srcX, Int_t srcY, UInt_t width, UInt_t height, Int_t dstX, Int_t dstY);
   void AddDrawString(Drawable_t wid, const GCValues_t &gc, Int_t x, Int_t y, const char *text, Int_t len);
   void AddFillRectangle(Drawable_t wid, const GCValues_t &gc, Int_t x, Int_t y, UInt_t w, UInt_t h);
   void AddFillPolygon(Drawable_t wid, const GCValues_t &gc, const Point_t *polygon, Int_t nPoints);
   void AddDrawRectangle(Drawable_t wid, const GCValues_t &gc, Int_t x, Int_t y, UInt_t w, UInt_t h);
   void AddUpdateWindow(QuartzView *view);
   void AddDeletePixmap(Pixmap_t pixmap);

   //'XOR' graphics for canvas.
   void AddDrawBoxXor(Window_t windowID, Int_t x1, Int_t y1, Int_t x2, Int_t y2);
   void AddDrawLineXor(Window_t windowID, Int_t x1, Int_t y1, Int_t x2, Int_t y2);

   void Flush(Details::CocoaPrivate *impl);
   void FlushXOROps(Details::CocoaPrivate *impl);
   void RemoveOperationsForDrawable(Drawable_t wid);
   void RemoveGraphicsOperationsForWindow(Window_t wid);
   void RemoveXORGraphicsOperationsForWindow(Window_t wid);

   size_type BufferSize()const
   {
      return fCommands.size();
   }

   void ClearXOROperations();
private:
   void ClearCommands();
   //Clip related stuff.
   struct WidgetRect {
      int fX1;
      int fY1;
      int fX2;
      int fY2;

      WidgetRect()
         : fX1(0), fY1(0), fX2(0), fY2(0)
      {
      }

      WidgetRect(int leftX, int bottomY, int rightX, int topY)
         : fX1(leftX), fY1(bottomY), fX2(rightX), fY2(topY)
      {
      }
   };

   void ClipOverlaps(QuartzView *view);
   void BuildClipRegion(const WidgetRect &rect);

   std::vector<WidgetRect> fRectsToClip;
   std::vector<CGRect> fClippedRegion;
   std::vector<int> fXBounds;
   std::vector<int> fYBounds;
   std::vector<bool> fGrid;
};

}//X11
}//MacOSX
}//ROOT

#endif
