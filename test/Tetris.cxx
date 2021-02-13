// @(#)root/test:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   04/10/98

///////////////////////////////////////////////////////////////////
//  ROOT implementation of the simple Tetris game
//  Layout and some hints were taken from Qt /examples/tetris
//
//  To run this game do the following:
//  $ root
//  root [0] gSystem.Load("libGpad")
//  root [1] gSystem.Load("Tetris")
//  root [2] Tetris t
//  <play game>
//  root [2] .q
//
///////////////////////////////////////////////////////////////////

#include <TVirtualX.h>
#include <TGClient.h>
#include <KeySymbols.h>
#include <TRootCanvas.h>
#include <TApplication.h>
#include <TList.h>
#include <snprintf.h>
#include "Tetris.h"

static Tetris *gTetris;                    // game manager

static const UInt_t gBoxPixelSize = 20;   // size of TetrisBox in pixels


ClassImp(Tetris);

///////////////////////////////////////////////////////////////////
//  TetrisBox - the main brick of the game
///////////////////////////////////////////////////////////////////
TetrisBox::TetrisBox(Int_t x, Int_t y, UInt_t type, TPad* pad) :
   TWbox(0,0,1,1,33,2,1)
{
   // Create brick

   fType = type;
   fPad  = pad;

   //-------  append box to pad
   SetBit(kMustCleanup);
   SetBit(kCanDelete);
   pad->GetListOfPrimitives()->Add(this);
   SetXY(x,y);
}

void TetrisBox::SetX(Int_t x)
{
   // Set X measured in boxes units

   // width in pixels of pad
   Float_t width  = (Float_t)fPad->XtoPixel(fPad->GetX2());

   Coord_t x1 = ((Float_t)x)*gBoxPixelSize/width;
   Coord_t x2 = ((Float_t)x+1)*gBoxPixelSize/width;

   SetX1(x1);
   SetX2(x2);
   fX = x;
}

void TetrisBox::SetY(Int_t y)
{
   // Set Y measured in boxes units

   // height in pixels of pad
   Float_t height = (Float_t)fPad->YtoPixel(fPad->GetY1());

   Coord_t y1 = ((Float_t)y)*gBoxPixelSize/height;
   Coord_t y2 = ((Float_t)y+1)*gBoxPixelSize/height;

   SetY1(y1);
   SetY2(y2);
   fY = y;
}

void TetrisBox::Paint(Option_t *option)
{
   // Paint box if it's not hidden.

   if (!IsHidden() && fPad) TWbox::Paint(option);
}

void TetrisBox::Erase()
{
   // erase box

   Double_t fX1sav = fX1;
   Double_t fY1sav = fY1;
   Double_t fX2sav = fX2;
   Double_t fY2sav = fY2;

   // erase 2 pix extra
   fX1 = fX1-fPad->PixeltoX(2);
   fY1 = fY1+fPad->PixeltoY(2);
   fX2 = fX2+fPad->PixeltoX(2);
   fY2 = fY2-fPad->PixeltoY(2);

   SetFillColor(fPad->GetFillColor());
   SetBorderMode(0);
   Paint();

   fX1 = fX1sav;
   fY1 = fY1sav;
   fX2 = fX2sav;
   fY2 = fY2sav;
}


///////////////////////////////////////////////////////////////////
//  TetrisPiece - tetris piece is set of up to 4 TetrisBoxes
///////////////////////////////////////////////////////////////////
static Int_t gPieceTypes[10][4][2] = {
                                      {{-1,1},   //   *
                                       {-1,0},   //   * *
                                       { 0,0},   //     *
                                       { 0,-1}},

                                      {{ 1, 1},  //    *
                                       { 1, 0},  //  * *
                                       { 0, 0},  //  *
                                       { 0,-1}},

                                      {{ 0, 1},  //   *
                                       { 0, 0},  //   *
                                       { 0,-1},  //   *
                                       { 0,-2}}, //   *

                                      {{ 0, 1},  //    *
                                       { 1, 0},  //  * * *
                                       { 0, 0},
                                       { -1,0}},

                                      {{ 1, 1},  //  * *
                                       { 0, 1},  //  * *
                                       { 1, 0},  //
                                       { 0, 0}}, //

                                      {{ 0, 1},   //   *
                                       { 0, 0},   //   *
                                       { 0,-1},   // * *
                                       { -1,-1}},

                                      {{ 0, 1},   //  *
                                       { 0, 0},   //  *
                                       { 0,-1},   //  * *
                                       { 1,-1}},

                                      {{ 0, 1},    // *
                                       { 0, 1},    //      hidden
                                       { 0, 1},    //      hidden
                                       { 0, 1}},   //      hidden

                                      {{ 0, 1},   //  *
                                       { 0, 0},   //  *
                                       { 0, 1},   //       hidden
                                       { 0, 0}},  //       hidden

                                      {{ 0, 1},   //   *
                                       { 0, 0},   //   *
                                       { 0,-1},   //   *
                                       { 0, 0}}}; //       hidden

static Color_t gPieceColors[10] = { 2,3,4,5,6,7,13,9,28,41 };


TetrisPiece::~TetrisPiece()
{
   // Clean up piece.

   for (int i = 0; i < 4; i++) delete fBoxes[i];
}

void TetrisPiece::Initialize(UInt_t type, TPad* pad)
{
   for (int i = 0; i < 4; i++) fBoxes[i] = new TetrisBox(0,0,0,pad);
   fX = 0;
   fY = 0;
   SetType(type);
}

void TetrisPiece::SetType(UInt_t type)
{
   // re-initialization

   if (type < 1 || type > 10)
      type = 9;

   for (int i = 0 ; i < 4 ; i++) {
      fBoxes[i]->SetType(type);
      fBoxes[i]->SetFillColor(gPieceColors[type-1]);
      fBoxes[i]->SetX(gPieceTypes[type-1][i][0]+fX);
      fBoxes[i]->SetY(gPieceTypes[type-1][i][1]+fY);
   }
   HideSomeBoxes(type);
   fType = type;
}

void TetrisPiece::HideSomeBoxes(UInt_t type)
{
   // Make invisible some boxes for 1,2,3 length pieces

   switch(type) {
      case 8:  fBoxes[1]->Hide();
      case 9:  fBoxes[2]->Hide();
      case 10: fBoxes[3]->Hide();
      default: return;
   }
}

Bool_t TetrisPiece::RotateRight()
{
   // Rotate anticlockwise  around (fX,fY) point

   // don't rotate square pieces
   if (GetType() == 5 || GetType() == 8) return kFALSE;

   Int_t tmp;

   for (int i = 0 ; i < 4 ; i++) {
      tmp = GetXx(i);
      SetXx(i,-GetYy(i));
      SetYy(i,tmp);
   }
   return kTRUE;
}

Bool_t TetrisPiece::RotateLeft()
{
   // Rotate clockwise  around (fX,fY) point

   // don't rotate square pieces
   if (GetType() == 5 || GetType() == 8) return kFALSE;

   Int_t tmp;

   for (int i = 0 ; i < 4 ; i++) {
      tmp = GetXx(i);
      SetXx(i,GetYy(i));
      SetYy(i,-tmp);
   }
   return kTRUE;
}

void  TetrisPiece::Hide()
{
   // Hide this piece

   for (int i = 0 ; i < 4 ; i++) {
      fBoxes[i]->Hide();
   }
}

void  TetrisPiece::Show()
{
   // Show this piece

   for (int i = 0 ; i < 4 ; i++) {
      fBoxes[i]->SetType(fType);
   }
   HideSomeBoxes(fType);
}

void  TetrisPiece::SetX(Int_t x)
{
   // Change of X position of the whole piece

   for (int i = 0 ; i < 4 ; i++) SetX(i,x+GetXx(i));
   fX = x;
}

void TetrisPiece::SetY(Int_t y)
{
   // Change of Y position of the whole piece

   for (int i = 0 ; i < 4 ; i++) SetY(i,y+GetYy(i));
   fY = y;
}

void TetrisPiece::SetXY(Int_t x,Int_t y)
{
   // Change of X,Y position of the whole piece

   for (int i = 0 ; i < 4 ; i++) {
      SetX(i,x+GetXx(i));
      SetY(i,y+GetYy(i));
   }
   fX = x;
   fY = y;
}


///////////////////////////////////////////////////////////////////
//  CurrentPiece =  TetrisPiece + TTimer = live TetrisPiece
///////////////////////////////////////////////////////////////////
CurrentPiece::CurrentPiece(UInt_t type,TetrisBoard* board) :
   TetrisPiece(type,board), TTimer(1000,kTRUE)
{
   // Initialize new piece

   fBoard = board;

   Int_t line      = fBoard->GetHeight()-2;
   Int_t xPosition = fBoard->GetWidth()/2;

   SetXY(xPosition,line);

   if (!CanMoveTo(xPosition,line))    { gTetris->StopGame(); return; }

   fBoard->Modified();
   fBoard->Update();
   fBoard->SetDropped(kFALSE);
   SetSpeed();       // set speed of moving according to game level
   Start();          // add this timer to sytem timers list = start moving
}

void CurrentPiece::MoveTo(int x, int y)
{
   // Move this to (x,y) and  draw it there

   Erase();
   SetXY(x,y);            // set new coordinates
   fBoard->Modified();    // drawing
   fBoard->Update();
}

Bool_t CurrentPiece::CanMoveTo(int x, int y)
{
   // Can move this piece to (x,y)?

   Bool_t return_value;

   int  savX = fX;
   int  savY = fY;

   SetXY(x,y);                    // set new coordinates
   return_value = CanPosition();  // if inside board and no non-zero squares underneath
   SetXY(savX,savY);              // go back

   return return_value;
}

Bool_t CurrentPiece::CanPosition()
{
   // Check if piece position is allowed

   int x, y;

   for (int i = 0 ; i < 4 ; i++) {
      GetXY(i,x,y);      // coordinates of piece boxes at test position
      if (x < 0 || x >= fBoard->GetWidth()  ||
          y < 0 || y >= fBoard->GetHeight() ||
          !fBoard->IsEmpty(x,y))          return kFALSE;
   }
   return kTRUE;      // Inside board and no non-zero squares underneath.
}

Bool_t CurrentPiece::RotateRight()
{
   // Rotate clockwise. Returns kTRUE if succeeded.

   Bool_t return_value;

   Erase();
   TetrisPiece::RotateRight();
   return_value = CanPosition();
   if (!return_value) TetrisPiece::RotateLeft();    // rotate back

   fBoard->Modified();    // drawing
   fBoard->Update();
   return return_value;
}

Bool_t CurrentPiece::RotateLeft()
{
   // Rotate anti clockwise. Returns kTRUE if succeeded.

   Bool_t return_value;

   Erase();
   TetrisPiece::RotateLeft();
   return_value = CanPosition();
   if (!return_value) TetrisPiece::RotateRight();     // rotate back

   fBoard->Modified();    // drawing
   fBoard->Update();
   return return_value;
}

Bool_t CurrentPiece::OneLineDown()
{
   // Move one line down. Returns kTRUE if succeeded.

   int y = GetY();
   int x = GetX();

   y--;
   if (!CanMoveTo(x,y))  return kFALSE;

   MoveTo(x,y);
   return kTRUE;
}

Bool_t CurrentPiece::DropDown()
{
   // Move the piece to lowest allowed line. Returns kTRUE if succeeded.

   int y = GetY();
   int x = GetX();
   int dropHeight = 0;

   while (CanMoveTo(x,--y)) dropHeight++;  //  find lower allowed line

   y++;
   MoveTo(x,y);         //  .. and move to
   Stop();              //  stop moving
   fBoard->PieceDropped(this, dropHeight);
   return kTRUE;
}

Bool_t CurrentPiece::MoveLeft(int steps)
{
   // Move piece to the left. Return kTRUE if succeeded.

   int y = GetY();
   int x = GetX();

   while(steps) {
      if (!CanMoveTo(--x ,y)) return kFALSE;  // can't move
      MoveTo(x,y);
      steps--;
   }
   return kTRUE;
}

Bool_t CurrentPiece::MoveRight(int steps)
{
   // Move piece to the right. Return kTRUE if succeeded.

   int y = GetY();
   int x = GetX();

   while(steps) {
      if (!CanMoveTo(++x,y)) return kFALSE;  // can't move
      MoveTo(x,y);
      steps--;
   }
   return kTRUE;
}

Bool_t CurrentPiece::Notify()
{
   // Actions after time out.

   if (OneLineDown()) {
      TTimer::Reset();
      return kFALSE;
   } else {
      Stop();                        // stop moving
      fBoard->PieceDropped(this,0);  // piece can't move -> stop moving update state of TetrisBoard
      return kTRUE;
  }
}

void CurrentPiece::SetSpeed()
{
   // Set speed according to level.

   const Int_t factor = 2;
   SetTime(1000/(1 + gTetris->GetLevel()*factor));
}

void  CurrentPiece::Paint(Option_t*)
{
   // Paint it in fBoard.

   TPad* padsav = (TPad*)TVirtualPad::Pad();
   fBoard->cd();

   for (int i = 0 ; i < 4 ; i++) {
      fBoxes[i]->SetBorderMode(1);
      fBoxes[i]->SetFillColor(gPieceColors[fType-1]);
      fBoxes[i]->Paint();
   }
   padsav->cd();
}

void CurrentPiece::Erase()
{
   // Erase = paint with the same FillColor as fBoard has

   TPad* padsav = (TPad*)TVirtualPad::Pad();
   fBoard->cd();

   for (int i = 0 ; i < 4 ; i++) {
      fBoxes[i]->Erase();
   }
   padsav->cd();
}


///////////////////////////////////////////////////////////////////
//  Game  board
///////////////////////////////////////////////////////////////////
TetrisBoard::TetrisBoard(Float_t xlow, Float_t ylow,Float_t xup,Float_t yup) :
   TPad("tetris_board","Tetris Board",xlow,ylow,xup,yup,17,4,-1)
{
   // Game board constructor.

   fWidth  = (int)(fMother->XtoAbsPixel(GetX2())*(xup-xlow))/gBoxPixelSize;
   fHeight = (int)(fMother->YtoAbsPixel(GetY1())*(yup-ylow))/gBoxPixelSize;
   Double_t box = fMother->PixeltoX(gBoxPixelSize);
   Double_t xx = xlow + box*fWidth;

   if (xx<xup) {
      xx += fMother->PixeltoX(1);
      SetPad(xlow, ylow, xx, yup);
   }

   fBoard = new TetrisBoxPtr[fWidth*fHeight];
   fFilledLines = 0;
   Clear();
}

void TetrisBoard::Clear(Option_t *)
{
   // Delete/clear all objects.

   GetListOfPrimitives()->Delete();   // delete all object in this pad (including TetrisPiece)

   for (int i = 0; i < fWidth; i++)
      for (int j = 0; j < fHeight; j++)
         Board(i,j) = 0;              // clear board map

   fIsDropped = kTRUE;
}

void TetrisBoard::Hide()
{
   // Hide all objects.

   TetrisBox *box;
   TIter nextin(GetListOfPrimitives());

   while ((box = (TetrisBox*)nextin())) box->Hide();
   Modified();
   Update();
}

void TetrisBoard::Show()
{
   // Show all objects

   TetrisBox *box;
   TIter nextin(GetListOfPrimitives());

   while ((box = (TetrisBox*)nextin())) box->Show();
   Modified();
   Update();
}

Bool_t TetrisBoard::IsFullLine(Int_t line)
{
   // Check if line is full.

   Bool_t return_value = kTRUE;

   for (int i = 0; i < fWidth; i++)
      return_value = return_value && !IsEmpty(i,line);

   return return_value;
}

Bool_t TetrisBoard::IsEmptyLine(Int_t line)
{
   // Check if line is empty

   Bool_t return_value = kTRUE;

   for (int i = 0; i < fWidth; i++)
      return_value = return_value && IsEmpty(i,line);

   return return_value;
}

void TetrisBoard::RemoveLine(Int_t line)
{
   // Remove all TetrisBoxes of the line

   for (int i=0; i<fWidth; i++) {
      if (Board(i,line))  // when you delete TObject it's also removed from Pad
      delete Board(i,line);

      Board(i,line) = 0;
   }
}

void TetrisBoard::MoveOneLineDown(Int_t line)
{
   // All  boxes of this line move to (line-1)

   if (!line) return;   // don't move line==0

   for (int i = 0; i < fWidth; i++) {
      if (!IsEmpty(i,line)) {
         Board(i,line)->MoveOneLineDown();     // change  position of Boxes
         Board(i,line-1) = Board(i,line);      // remapping
      }
      Board(i,line)=0;   // this line become empty
   }
   Modified();
   Update();
}

void TetrisBoard::RemoveFullLines()
{
   // Remove full lines

   for (int i = 0; i < fFilledLines; i++) {
      while (IsFullLine(i)) {
         RemoveLine(i);
         gTetris->UpdateLinesRemoved();
         AllAboveLinesDown(i);
         fFilledLines--;
      }
   }
}

void  TetrisBoard::GluePiece(TetrisPiece* piece)
{
   // Add pointers to piece boxes to fBoard::fBoard

   int x,y;
   TetrisBox *box;

   for (int i = 0 ; i < 4 ; i++) {
     piece->GetXY(i,x,y);
     box = piece->GetTetrisBox(i);
     if (box->IsHidden()) { delete  box; continue;}     // delete hidden boxes
     Board(x,y) = piece->GetTetrisBox(i);     // add pointers to piece boxes to board map
     if (y>fFilledLines) fFilledLines = y+1;  // update number of non empty lines
   }
}

void TetrisBoard::PieceDropped(TetrisPiece* piece, int height)
{
   // Actions after piece was droped

   Int_t add2score = height*gTetris->GetLevel() + 10;    // update score policy (could be modified)

   fIsDropped = kTRUE;
   GluePiece(piece);
   RemoveFullLines();
   //Print();             // possible printig of board map on the terminal

   gTetris->UpdatePiecesDropped();
   gTetris->UpdateScore(add2score);
   gTetris->CreateNewPiece();        // create new CurrentPiece
   fIsDropped = kFALSE;
}

void TetrisBoard::Print(const char *) const
{
   // Used for testing

   printf("\n");

   for (int j = fHeight-1; j > -1; j--) {
      for (int i = 0; i < fWidth; i++)
         ((TetrisBoard*)this)->IsEmpty(i,j) ? printf("|   ") : printf("| * ") ;
      printf("|\n");
   }
}

void TetrisBoard::PaintModified()
{
   // Overload this method to acelerate graphics
   // (do not draw tens of heap boxes while current box is moving)

   if (!fIsDropped && gTetris->IsGameOn() && !gTetris->IsPaused())
      gTetris->fCurrentPiece->Paint();
   else
      TPad::PaintModified();
}

///////////////////////////////////////////////////////////////////
//  NextPiecePad
//  used to show next piece.
///////////////////////////////////////////////////////////////////
NextPiecePad::NextPiecePad(Float_t xlow, Float_t ylow, Float_t xup, Float_t yup)
   : TPad("next_piece","Next Piece Pad",xlow,ylow,xup,yup,17,4,-1)
{
   // Next piece pad ctor.

   fPiece = new TetrisPiece(this);
   fPiece->Hide();        // hide piece at start

   // (how to get pixel size?)
   Int_t x = (int)(fMother->XtoAbsPixel(GetX2())*(xup-xlow))/gBoxPixelSize/2;
   Int_t y = (int)(fMother->YtoAbsPixel(GetY1())*(yup-ylow))/gBoxPixelSize/2;

   fPiece->SetXY(x,y);    // move to the center of pad
   Modified(kTRUE);
   Update();
}


///////////////////////////////////////////////////////////////////
//  PauseButton - push button
//  ExecuteEvent mehtod used to pause the game
///////////////////////////////////////////////////////////////////
PauseButton::PauseButton(Float_t xlow, Float_t ylow, Float_t xup, Float_t yup) :
   TButton("Pause"," ",xlow,ylow,xup,yup)
{
   // Pause button constructor

   SetBorderSize(5);     //  decoration stuff....
   SetTextSize(0.45);
   SetFillColor(42);
}

void PauseButton::ExecuteEvent(Int_t event, Int_t, Int_t)
{
   // Action after mouse click

   if (event == kButton1Up) {
      IsPressed() ? gTetris->Continue() :  gTetris->Pause();
      Modified(kTRUE);
   }
}


///////////////////////////////////////////////////////////////////
//  QuitButton - push button
//  ExecuteEvent mehtod used to quit the game
///////////////////////////////////////////////////////////////////
QuitButton::QuitButton(Float_t xlow, Float_t ylow, Float_t xup, Float_t yup) :
   TButton("Quit"," ",xlow,ylow,xup,yup)
{
   // Quit button constructor

   SetBorderSize(5);   //  decoration stuff....
   SetTextSize(0.45);
   SetFillColor(42);
}

void QuitButton::ExecuteEvent(Int_t event, Int_t, Int_t)
{
   // Action after mouse click

   if (event == kButton1Up) gTetris->Quit();
}


///////////////////////////////////////////////////////////////////
//  NewGameButton - push button
//  ExecuteEvent mehtod used to start new game the game
///////////////////////////////////////////////////////////////////
NewGameButton::NewGameButton(Float_t xlow, Float_t ylow, Float_t xup, Float_t yup)
   : TButton("New Game"," ",xlow,ylow,xup,yup)
{
   // New game button constructor

   SetBorderSize(5);
   SetTextSize(0.45);
   SetFillColor(42);
}

void NewGameButton::ExecuteEvent(Int_t event, Int_t, Int_t)
{
   // Ation after mouse click

   if (event == kButton1Up) {
      gTetris->NewGame();      // always starts new game
      Modified(kTRUE);
   }
}


///////////////////////////////////////////////////////////////////
//  InfoPad -
///////////////////////////////////////////////////////////////////
InfoPad::InfoPad(const char* title, Float_t xlow, Float_t ylow, Float_t xup, Float_t yup)
   : TPad("info_pad",title,xlow,ylow,xup,yup,17,4,-1), TAttText(22,0,2,71,0.65)
{
   // InfoPad constructor

   SetBit(kCanDelete);

   TText *text = new TText(xlow,yup,title);   // draw title of the information pad
   text->SetTextSize(0.45*(yup-ylow));
   text->SetY(yup+0.2*text->GetTextSize());
   fMother->GetListOfPrimitives()->Add(text);

   text = new TText(0.5,0.5,"0");          // this text used to display fValue
   GetListOfPrimitives()->Add(text);

   fValue = 0;
   Modified(kTRUE);
   Update();
}

void InfoPad::PaintModified()
{
   // Actions after pad was modified (resize event, user's Modified(kTRUE) ...)

   char    str[40];

   snprintf(str,40,"%d",fValue);

   TObject *obj = GetListOfPrimitives()->First();

   if (obj && obj->InheritsFrom("TText")) {
      TText *text = (TText*)obj;
      text->SetTitle(str);                // set title according to fValue

      text->SetTextSize(GetTextSize());
      text->SetTextFont(GetTextFont());
      text->SetTextAlign(GetTextAlign());
      text->SetTextColor(GetTextColor());
      text->SetTextAngle(GetTextAngle());
      text->TAttText::Modify();

      text->SetX(0.5);        // draw centered
      text->SetY(0.5);
   }
   TPad::PaintModified();
}


///////////////////////////////////////////////////////////////////
//  KeyHandler - virtual frame used to catch and handle key events
///////////////////////////////////////////////////////////////////
KeyHandler::KeyHandler() : TGFrame(gClient->GetRoot(),0,0)
{
   // Key handler constructor.

   // get main frame of Tetris canvas
   TRootCanvas *main_frame = (TRootCanvas*)(gTetris->GetCanvasImp());

   // bind arrow keys and space-bar key
   main_frame->BindKey((const TGWindow*)this, gVirtualX->KeysymToKeycode(kKey_Up),    0);
   main_frame->BindKey((const TGWindow*)this, gVirtualX->KeysymToKeycode(kKey_Left),  0);
   main_frame->BindKey((const TGWindow*)this, gVirtualX->KeysymToKeycode(kKey_Right), 0);
   main_frame->BindKey((const TGWindow*)this, gVirtualX->KeysymToKeycode(kKey_Down),  0);
   main_frame->BindKey((const TGWindow*)this, gVirtualX->KeysymToKeycode(kKey_Space), 0);
}

KeyHandler::~KeyHandler()
{
   // Cleanup key handler.

   // get main frame of Tetris canvas
   TRootCanvas *main_frame = (TRootCanvas*)(gTetris->GetCanvasImp());

   // remove binding of arrow keys and space-bar key
   main_frame->RemoveBind(this, gVirtualX->KeysymToKeycode(kKey_Up),    0);
   main_frame->RemoveBind(this, gVirtualX->KeysymToKeycode(kKey_Left),  0);
   main_frame->RemoveBind(this, gVirtualX->KeysymToKeycode(kKey_Right), 0);
   main_frame->RemoveBind(this, gVirtualX->KeysymToKeycode(kKey_Down),  0);
   main_frame->RemoveBind(this, gVirtualX->KeysymToKeycode(kKey_Space), 0);
   // restore key auto repeat functionality, was turned off in TGMainFrame::HandleKey()
   gVirtualX->SetKeyAutoRepeat(kTRUE);
}


Bool_t KeyHandler::HandleKey(Event_t *event)
{
   // Handle arrow and spacebar keys

   char tmp[2];
   UInt_t keysym;

   gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);

   if (event->fType == kGKeyPress) {
      switch ((EKeySym)keysym) {
         case kKey_Left:
            gTetris->MoveLeft();
            break;
         case kKey_Right:
            gTetris->MoveRight();
            break;
         case kKey_Down:
            gTetris->RotateRight();
            break;
         case kKey_Up:
            gTetris->RotateLeft();
            break;
         case kKey_Space:
            gTetris->DropDown();
            break;
         default:
            return kTRUE;
      }
   }
   return kTRUE;
}

///////////////////////////////////////////////////////////////////
//  UpdateLevelTimer
///////////////////////////////////////////////////////////////////
UpdateLevelTimer::UpdateLevelTimer(ULong_t time) : TTimer(time,kTRUE)
{
   // Update level timer constructor

   SetBit(kCanDelete);   // delete this when gTetris is deleted
   gTetris->GetListOfPrimitives()->Add(this);
}

Bool_t UpdateLevelTimer::Notify()
{
   // Actions after time out

   if (!gTetris->IsGameOn()) {
      Remove();
      return kTRUE;
   }
   gTetris->UpdateLevel();
   TTimer::Reset();
   return kFALSE;
}


///////////////////////////////////////////////////////////////////
//   Tetris =  Game manager
///////////////////////////////////////////////////////////////////
Tetris::Tetris() :
   TCanvas("Tetris","Have a little fun with ROOT!",200,200,700,500)
{
   // Tetris constructor

   gTetris = this;

   fCurrentPiece = 0;

   //-----------  play board ------------
   fBoard            = new  TetrisBoard(0.35,0.05,0.7,0.95);
   fBoard->Draw();

   // ----------  info pads -------------
   fNextPiece        =  new  NextPiecePad(0.05,0.65,0.25,0.95);
   fNextPiece->Draw();

   fLinesRemoved     =  new  InfoPad("Lines Removed",0.75,0.8,0.95,0.9);
   fLinesRemoved->Draw();

   fLevel            =  new  InfoPad("Level",0.75,0.6,0.95,0.7);
   fLevel->Draw();

   fScore            =  new  InfoPad("Score",0.75,0.4,0.95,0.5);
   fScore->Draw();

   //------------ buttons ----------------
   fNewGame          =  new  NewGameButton(0.05,0.05,0.25,0.15);
   fNewGame->Draw();

   fQuit             =  new  QuitButton(0.05,0.2,0.25,0.3);
   fQuit->Draw();

   fPause            =  new  PauseButton(0.05,0.35,0.25,0.45);
   fPause->Draw();

   fPiecesDropped = 0;
   SetFillColor(21);

   fKeyHandler = new KeyHandler();
   fUpdateLevelTimer = new UpdateLevelTimer(60000);  // every  minute
   SetFixedSize();
   Update();
   PrintHelpInfo();
   fEditable = kFALSE;
}

void Tetris::PrintHelpInfo()
{
   // Prints help info

   printf("\n\n\n");
   printf("             Move   Piece Left     ---------     left-arrow\n");
   printf("             Move   Piece Right    ---------     right-arrow\n");
   printf("             Rotate Piece          ---------     up/down-arrow \n");
   printf("             Drop   Piece Down     ---------     space-bar\n");
   printf("\n\n\n");
}

void Tetris::CreateNewPiece()
{
   // Create  current and next pieces

   UInt_t type = fNextPiece->GetPiece()->GetType();
   fNextPiece->NewPiece();
   fCurrentPiece = new CurrentPiece(type,fBoard);
}

void Tetris::SetFixedSize()
{
   // Set size of canvas

   ((TRootCanvas*)fCanvasImp)->SetWMSizeHints(fCw,fCh+20,fCw,fCh,0,0);
}

void Tetris::Quit()
{
   // Stop game and delete canvas (i.e. tetris itself)

   delete fKeyHandler; fKeyHandler = 0;
   StopGame();
   ((TRootCanvas*)fCanvasImp)->CloseWindow();
}

void Tetris::NewGame()
{
   // Start new game

   gVirtualX->SetInputFocus(((TRootCanvas*)fCanvasImp)->GetId());

   if (IsGameOn()) StopGame();       // stop privious game
   fScore->Reset();
   fLinesRemoved->Reset();
   fPiecesDropped = 0;
   SetLevel(1);
   fUpdateLevelTimer->Start();
   fBoard->Clear();
   fNewGame->SetPressed(kTRUE);
   CreateNewPiece();                 // start game
}

void Tetris::StopGame()
{
   // Stop the game

   fUpdateLevelTimer->Stop();
   if (fCurrentPiece) fCurrentPiece->Stop();
   fNewGame->SetPressed(kFALSE);
   fPause->SetPressed(kFALSE);
}

void Tetris::Pause()
{
   // Pause the game

   if (!IsGameOn())  return;
   if (fCurrentPiece) fCurrentPiece->Stop();
   fPause->SetPressed(kTRUE);
   fBoard->Hide();
}

void  Tetris::Continue()
{
   // Continue the game

   if (!IsGameOn()) return;
   fBoard->Show();
   fPause->SetPressed(kFALSE);
   if (fCurrentPiece) fCurrentPiece->Start();
}

void Tetris::MoveLeft()
{
   // Move pice to the left

   if (!IsGameOn() || IsPaused() || IsWaiting())  return;
   fCurrentPiece->MoveLeft();
}

void Tetris::MoveRight()
{
   // Move piece to the right

   if (!IsGameOn() || IsPaused() || IsWaiting())  return;
   fCurrentPiece->MoveRight();
}

void Tetris::DropDown()
{
   // Drop piece down

   if (!IsGameOn() || IsPaused() || IsWaiting())  return;
   fCurrentPiece->DropDown();
}

void Tetris::RotateRight()
{
   // Rotate piece right

   if (!IsGameOn() || IsPaused() || IsWaiting())  return;
   fCurrentPiece->RotateRight();
 }

void Tetris::RotateLeft()
{
   // Rotate piece left

   if (!IsGameOn() || IsPaused() || IsWaiting())  return;
   fCurrentPiece->RotateLeft();
 }

void Tetris::SetLevel(int level)
{
   // Set difficulty level

   fLevel->SetValue(level);
   if (fCurrentPiece) fCurrentPiece->SetSpeed();
}
