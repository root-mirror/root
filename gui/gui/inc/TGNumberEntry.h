// @(#)root/gui:$Id$
// Author: Daniel Sigg   03/09/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TGNumberEntry
#define ROOT_TGNumberEntry

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGNumberEntry, TGNumberEntryField and TGNumberFormat                 //
//                                                                      //
// TGNumberEntry is a number entry input widget with up/down buttons.   //
// TGNumberEntryField is a number entry input widget.                   //
// TGNumberFormat contains enum types to specify the numeric format .   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"
#include "TGTextEntry.h"
#include "TGButton.h"


class TGNumberFormat {
public:
   enum EStyle {             // Style of number entry field
      kNESInteger = 0,       // Integer
      kNESRealOne = 1,       // Fixed fraction real, one digit
      kNESRealTwo = 2,       // Fixed fraction real, two digit
      kNESRealThree = 3,     // Fixed fraction real, three digit
      kNESRealFour = 4,      // Fixed fraction real, four digit
      kNESReal = 5,          // Real number
      kNESDegree = 6,        // Degree
      kNESMinSec = 7,        // Minute:seconds
      kNESHourMin = 8,       // Hour:minutes
      kNESHourMinSec = 9,    // Hour:minute:seconds
      kNESDayMYear = 10,     // Day/month/year
      kNESMDayYear = 11,     // Month/day/year
      kNESHex = 12           // Hex
   };

   enum EAttribute {         // Attributes of number entry field
      kNEAAnyNumber = 0,     // Any number
      kNEANonNegative = 1,   // Non-negative number
      kNEAPositive = 2       // Positive number
   };

   enum ELimit {             // Limit selection of number entry field
      kNELNoLimits = 0,      // No limits
      kNELLimitMin = 1,      // Lower limit only
      kNELLimitMax = 2,      // Upper limit only
      kNELLimitMinMax = 3    // Both lower and upper limits
   };

   enum EStepSize {          // Step for number entry field increase
      kNSSSmall = 0,         // Small step
      kNSSMedium = 1,        // Medium step
      kNSSLarge = 2,         // Large step
      kNSSHuge = 3           // Huge step
   };

   virtual ~TGNumberFormat() { }
   ClassDef(TGNumberFormat,0)  // Class defining namespace for several enums used by TGNumberEntry
};


class TGNumberEntryField : public TGTextEntry, public TGNumberFormat {

protected:
   Bool_t        fNeedsVerification; // Needs verification of input
   EStyle        fNumStyle;          // Number style
   EAttribute    fNumAttr;           // Number attribute
   ELimit        fNumLimits;         // Limit attributes
   Double_t      fNumMin;            // Lower limit
   Double_t      fNumMax;            // Upper limit
   Bool_t        fStepLog;           // Logarithmic steps for increase?

public:
   TGNumberEntryField(const TGWindow *p, Int_t id,
                      Double_t val, GContext_t norm,
                      FontStruct_t font = GetDefaultFontStruct(),
                      UInt_t option = kSunkenFrame | kDoubleBorder,
                      Pixel_t back = GetWhitePixel());
   TGNumberEntryField(const TGWindow *parent = 0,
                      Int_t id = -1, Double_t val = 0,
                      EStyle style = kNESReal,
                      EAttribute attr = kNEAAnyNumber,
                      ELimit limits = kNELNoLimits,
                      Double_t min = 0, Double_t max = 1);

   virtual void SetNumber(Double_t val);
   virtual void SetIntNumber(Long_t val);
   virtual void SetTime(Int_t hour, Int_t min, Int_t sec);
   virtual void SetDate(Int_t year, Int_t month, Int_t day);
   virtual void SetHexNumber(ULong_t val);
   void SetText(const char* text, Bool_t emit = kTRUE) override;

   virtual Double_t GetNumber() const;
   virtual Long_t   GetIntNumber() const;
   virtual void     GetTime(Int_t& hour, Int_t& min, Int_t& sec) const;
   virtual void     GetDate(Int_t& year, Int_t& month, Int_t& day) const;
   virtual ULong_t  GetHexNumber() const;

   virtual Int_t GetCharWidth(const char* text = "0") const;
   virtual void  IncreaseNumber(EStepSize step = kNSSSmall,
                                Int_t sign = 1, Bool_t logstep = kFALSE);
   virtual void  SetFormat(EStyle style,
                           EAttribute attr = kNEAAnyNumber);
   virtual void  SetLimits(ELimit limits = kNELNoLimits,
                           Double_t min = 0, Double_t max = 1);
   void  SetState(Bool_t state) override;
   virtual void  SetLogStep(Bool_t on = kTRUE) {
      // Set logarithmic steps
      fStepLog = on; }

   virtual EStyle GetNumStyle() const {
      // Get the numerical style
      return fNumStyle; }
   virtual EAttribute GetNumAttr() const {
      // Get the numerical attribute
      return fNumAttr; }
   virtual ELimit GetNumLimits() const {
      // Get the numerialc limit attribute
      return fNumLimits; }
   virtual Double_t GetNumMin() const {
      // Get the lower limit
      return fNumMin; }
   virtual Double_t GetNumMax() const {
      // Get the upper limit
      return fNumMax; }
   virtual Bool_t IsLogStep() const {
      // Is log step enabled?
      return fStepLog; }

   Bool_t HandleKey(Event_t* event) override;
   Bool_t HandleFocusChange (Event_t* event) override;
   void   TextChanged(const char *text = 0) override;
   void   ReturnPressed() override;
   void   Layout() override;
   Bool_t IsEditable() const override { return kFALSE; }
   virtual void   InvalidInput(const char *instr) { Emit("InvalidInput(char*)", instr); }   //*SIGNAL*
   void   SavePrimitive(std::ostream &out, Option_t * = "") override;

   ClassDef(TGNumberEntryField,0)  // A text entry field used by a TGNumberEntry
};



class TGNumberEntry : public TGCompositeFrame, public TGWidget,
   public TGNumberFormat {

   // dummy data members - just to say about options for context menu
   EStyle fNumStyle;//*OPTION={GetMethod="GetNumStyle";SetMethod="SetNumStyle";Items=(0="Int",5="Real",6="Degree",9="Hour:Min:Sec",10="Day/Month/Year",12="Hex")}*
   EAttribute fNumAttr; // *OPTION={GetMethod="GetNumAttr";SetMethod="SetNumAttr";Items=(0="&AnyNumber",1="&Non negative",2="&Positive")}*
   ELimit fNumLimits; // *OPTION={GetMethod="GetNumLimits";SetMethod="SetNumLimits";Items=(0="&No Limits",1="Limit M&in",2="Limit M&ax",2="Min &and Max")}*

private:
   const TGPicture  *fPicUp;      // Up arrow
   const TGPicture  *fPicDown;    // Down arrow

   TGNumberEntry(const TGNumberEntry&);             // not implemented
   TGNumberEntry& operator=(const TGNumberEntry&);  // not implemented

protected:
   TGNumberEntryField *fNumericEntry;  // Number text entry field
   TGButton           *fButtonUp;      // Button for increasing value
   TGButton           *fButtonDown;    // Button for decreasing value
   Bool_t              fButtonToNum;   // Send button messages to parent rather than number entry field

public:
   TGNumberEntry(const TGWindow *parent = 0, Double_t val = 0,
                 Int_t digitwidth = 5, Int_t id = -1,
                 EStyle style = kNESReal,
                 EAttribute attr = kNEAAnyNumber,
                 ELimit limits = kNELNoLimits,
                 Double_t min = 0, Double_t max = 1);
   ~TGNumberEntry() override;

   virtual void SetNumber(Double_t val) {
      // Set the numeric value (floating point representation)
      fNumericEntry->SetNumber(val); }
   virtual void SetIntNumber(Long_t val) {
      // Set the numeric value (integer representation)
      fNumericEntry->SetIntNumber(val); }
   virtual void SetTime(Int_t hour, Int_t min, Int_t sec) {
      // Set the numeric value (time format)
      fNumericEntry->SetTime(hour, min, sec); }
   virtual void SetDate(Int_t year, Int_t month, Int_t day) {
      // Set the numeric value (date format)
      fNumericEntry->SetDate(year, month, day); }
   virtual void SetHexNumber(ULong_t val) {
      // Set the numeric value (hex format)
      fNumericEntry->SetHexNumber(val); }
   virtual void SetText(const char* text) {
      // Set the value (text format)
      fNumericEntry->SetText(text); }
   virtual void SetState(Bool_t enable = kTRUE);

   virtual Double_t GetNumber() const {
      // Get the numeric value (floating point representation)
      return fNumericEntry->GetNumber(); }
   virtual Long_t GetIntNumber() const {
      // Get the numeric value (integer representation)
      return fNumericEntry->GetIntNumber (); }
   virtual void GetTime(Int_t& hour, Int_t& min, Int_t& sec) const {
      // Get the numeric value (time format)
      fNumericEntry->GetTime(hour, min, sec); }
   virtual void GetDate(Int_t& year, Int_t& month, Int_t& day) const {
      // Get the numeric value (date format)
      fNumericEntry->GetDate(year, month, day); }
   virtual ULong_t GetHexNumber() const {
      // Get the numeric value (hex format)
      return fNumericEntry->GetHexNumber(); }
   virtual void IncreaseNumber(EStepSize step = kNSSSmall,
                               Int_t sign = 1, Bool_t logstep = kFALSE) {
      // Increase the number value
      fNumericEntry->IncreaseNumber(step, sign, logstep); }
   virtual void SetFormat(EStyle style, EAttribute attr = TGNumberFormat::kNEAAnyNumber) {
      // Set the numerical format
      fNumericEntry->SetFormat(style, attr); }
   virtual void SetLimits(ELimit limits = TGNumberFormat::kNELNoLimits,
                          Double_t min = 0, Double_t max = 1) {
      // Set the numerical limits.
      fNumericEntry->SetLimits(limits, min, max); }

   virtual EStyle GetNumStyle() const {
      // Get the numerical style
      return fNumericEntry->GetNumStyle(); }
   virtual EAttribute GetNumAttr() const {
      // Get the numerical attribute
      return fNumericEntry->GetNumAttr(); }
   virtual ELimit GetNumLimits() const {
      // Get the numerical limit attribute
      return fNumericEntry->GetNumLimits(); }
   virtual Double_t GetNumMin() const {
      // Get the lower limit
      return fNumericEntry->GetNumMin(); }
   virtual Double_t GetNumMax() const {
      // Get the upper limit
      return fNumericEntry->GetNumMax(); }
   virtual Bool_t IsLogStep() const {
      // Is log step enabled?
      return fNumericEntry->IsLogStep(); }
   virtual void   SetButtonToNum(Bool_t state);

   void SetNumStyle(EStyle style) {
         SetFormat(style, GetNumAttr()); }                  //*SUBMENU*
   void SetNumAttr(EAttribute attr = kNEAAnyNumber) {
         SetFormat(GetNumStyle(), attr); }                  //*SUBMENU*
   void SetNumLimits(ELimit limits = kNELNoLimits) {
         SetLimits(limits, GetNumMin(), GetNumMax());  }    //*SUBMENU*
   void SetLimitValues(Double_t min = 0, Double_t max = 1) {
         SetLimits(GetNumLimits(), min, max);  }            //*MENU*
   virtual void SetLogStep(Bool_t on = kTRUE);              //*TOGGLE* *GETTER=IsLogStep

   void   Associate(const TGWindow *w) override;
   Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override;
   virtual void   ValueChanged(Long_t val);     //*SIGNAL*
   virtual void   ValueSet(Long_t val);         //*SIGNAL*

   TGNumberEntryField *GetNumberEntry() const {
      // Get the number entry field
      return fNumericEntry; }
   TGButton *GetButtonUp() const {
      // Get the up button
      return fButtonUp; }
   TGButton *GetButtonDown() const {
      // Get the down button
      return fButtonDown; }

   Bool_t IsEditable() const override { return kFALSE; }

   UInt_t GetDefaultHeight() const override { return fNumericEntry->GetDefaultHeight(); }
   void SavePrimitive(std::ostream &out, Option_t * = "") override;
   TGLayoutManager *GetLayoutManager() const override;

   ClassDef(TGNumberEntry,0)  // Entry field widget for several numeric formats
};


class TGNumberEntryLayout : public TGLayoutManager {
protected:
   TGNumberEntry *fBox;        // pointer to numeric control box

private:
   TGNumberEntryLayout(const TGNumberEntryLayout&);             // not implemented
   TGNumberEntryLayout& operator=(const TGNumberEntryLayout&);  // not implemented

public:
   TGNumberEntryLayout(TGNumberEntry *box): fBox(box) { }
   void Layout() override;
   TGDimension GetDefaultSize() const override;

   ClassDef(TGNumberEntryLayout,0)  // Layout manager for number entry widget
};


#endif
