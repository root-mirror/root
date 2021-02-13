/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RFrame
#define ROOT7_RFrame

#include "ROOT/RDrawable.hxx"

#include "ROOT/RDrawableRequest.hxx"
#include "ROOT/RAttrLine.hxx"
#include "ROOT/RAttrFill.hxx"
#include "ROOT/RAttrMargins.hxx"
#include "ROOT/RAttrAxis.hxx"
#include "ROOT/RAttrValue.hxx"
#include "ROOT/RPadUserAxis.hxx"

#include <memory>
#include <map>

class TRootIOCtor;

namespace ROOT {
namespace Experimental {


/** \class RFrame
\ingroup GpadROOT7
\brief Holds an area where drawing on user coordinate-system can be performed.
\authors Axel Naumann <axel@cern.ch> Sergey Linev <s.linev@gsi.de>
\date 2017-09-26
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RFrame : public RDrawable  {

   friend class RPadBase;

public:

   class RUserRanges {
      std::vector<double> values;  ///< min/max values for all dimensions
      std::vector<bool> flags;     ///< flag if values available

      void UpdateDim(unsigned ndim, const RUserRanges &src)
      {
         if (src.IsUnzoom(ndim)) {
            ClearMinMax(ndim);
         } else {
            if (src.HasMin(ndim))
               AssignMin(ndim, src.GetMin(ndim));
            if (src.HasMax(ndim))
               AssignMax(ndim, src.GetMax(ndim));
         }
      }

   public:
      // Default constructor - for I/O
      RUserRanges() = default;

      // Constructor for 1-d ranges
      RUserRanges(double xmin, double xmax)
      {
         AssignMin(0, xmin);
         AssignMax(0, xmax);
      }

      // Constructor for 2-d ranges
      RUserRanges(double xmin, double xmax, double ymin, double ymax)
      {
         AssignMin(0, xmin);
         AssignMax(0, xmax);
         AssignMin(1, ymin);
         AssignMax(1, ymax);
      }

      // Extend number of dimensions which can be stored in the object
      void Extend(unsigned ndim = 3)
      {
         if (ndim*2 > values.size()) {
            values.resize(ndim*2, 0.);
            flags.resize(ndim*2, false);
         }
      }

      bool HasMin(unsigned ndim) const { return (ndim*2 < flags.size()) && flags[ndim*2]; }
      double GetMin(unsigned ndim) const { return (ndim*2 < values.size()) ? values[ndim*2] : 0.; }

      // Assign minimum for specified dimension
      void AssignMin(unsigned ndim, double value)
      {
          Extend(ndim+1);
          values[ndim*2] = value;
          flags[ndim*2] = true;
      }

      bool HasMax(unsigned ndim) const { return (ndim*2+1 < flags.size()) && flags[ndim*2+1]; }
      double GetMax(unsigned ndim) const { return (ndim*2+1 < values.size()) ? values[ndim*2+1] : 0.; }

      // Assign maximum for specified dimension
      void AssignMax(unsigned ndim, double value)
      {
          Extend(ndim+1);
          values[ndim*2+1] = value;
          flags[ndim*2+1] = true;
      }

      void ClearMinMax(unsigned ndim)
      {
         if (ndim*2+1 < flags.size())
            flags[ndim*2] = flags[ndim*2+1] = false;

         if (ndim*2+1 < values.size())
            values[ndim*2] = values[ndim*2+1] = 0.;
      }

      /** Returns true if axis configured as unzoomed, can be specified from client */
      bool IsUnzoom(unsigned ndim) const
      {
         return (ndim*2+1 < flags.size()) && (ndim*2+1 < values.size()) &&
               !flags[ndim*2] && !flags[ndim*2+1] &&
               (values[ndim*2] < -0.5) && (values[ndim*2+1] < -0.5);
      }

      // Returns true if any value is specified
      bool IsAny() const
      {
         for (auto fl : flags)
            if (fl) return true;
         return false;
      }

      void Update(const RUserRanges &src)
      {
         UpdateDim(0, src);
         UpdateDim(1, src);
         UpdateDim(2, src);
      }

   };

private:
   RAttrMargins fMargins{this, "margin"};         ///<! frame margins relative to pad
   RAttrLine fAttrBorder{this, "border"};         ///<! line attributes for border
   RAttrFill fAttrFill{this, "fill"};             ///<! fill attributes for the frame
   RAttrAxis fAttrX{this, "x"};                   ///<! drawing attributes for X axis
   RAttrAxis fAttrY{this, "y"};                   ///<! drawing attributes for Y axis
   RAttrAxis fAttrZ{this, "z"};                   ///<! drawing attributes for Z axis
   RAttrValue<bool> fGridX{this, "gridx", false}; ///<! show grid for X axis
   RAttrValue<bool> fGridY{this, "gridy", false}; ///<! show grid for Y axis
   RAttrValue<bool> fSwapX{this, "swapx", false}; ///<! swap position of X axis
   RAttrValue<bool> fSwapY{this, "swapy", false}; ///<! swap position of Y axis
   RAttrValue<int> fTicksX{this, "ticksx", 1};    ///<! X ticks drawing:
   RAttrValue<int> fTicksY{this, "ticksy", 1};    ///<! Y ticks drawing
   std::map<unsigned, RUserRanges> fClientRanges; ///<! individual client ranges

   /// Mapping of user coordinates to normal coordinates, one entry per dimension.
   std::vector<std::unique_ptr<RPadUserAxisBase>> fUserCoord;

   RFrame(const RFrame &) = delete;
   RFrame &operator=(const RFrame &) = delete;

   // Default constructor
   RFrame() : RDrawable("frame")
   {
      GrowToDimensions(2);
   }

   /// Constructor taking user coordinate system, position and extent.
   explicit RFrame(std::vector<std::unique_ptr<RPadUserAxisBase>> &&coords);

   void SetClientRanges(unsigned connid, const RUserRanges &ranges, bool ismainconn);

protected:

   void PopulateMenu(RMenuItems &) override;

   void GetAxisRanges(unsigned ndim, const RAttrAxis &axis, RUserRanges &ranges) const;
   void AssignZoomRange(unsigned ndim, RAttrAxis &axis, const RUserRanges &ranges);

public:

   class RZoomRequest : public RDrawableRequest {
      RUserRanges ranges; // specified ranges
   public:
      RZoomRequest() = default;
      std::unique_ptr<RDrawableReply> Process() override
      {
         auto frame = dynamic_cast<RFrame *>(GetContext().GetDrawable());
         if (frame) frame->SetClientRanges(GetContext().GetConnId(), ranges, GetContext().IsMainConn());
         return nullptr;
      }
   };


   RFrame(TRootIOCtor*) : RFrame() {}

   const RAttrMargins &GetMargins() const { return fMargins; }
   RFrame &SetMargins(const RAttrMargins &margins) { fMargins = margins; return *this; }
   RAttrMargins &Margins() { return fMargins; }

   const RAttrLine &GetAttrBorder() const { return fAttrBorder; }
   RFrame &SetAttrBorder(const RAttrLine &border) { fAttrBorder = border; return *this; }
   RAttrLine &AttrBorder() { return fAttrBorder; }

   const RAttrFill &GetAttrFill() const { return fAttrFill; }
   RFrame &SetAttrFill(const RAttrFill &fill) { fAttrFill = fill; return *this; }
   RAttrFill &AttrFill() { return fAttrFill; }

   const RAttrAxis &GetAttrX() const { return fAttrX; }
   RFrame &SetAttrX(const RAttrAxis &axis) { fAttrX = axis; return *this; }
   RAttrAxis &AttrX() { return fAttrX; }

   const RAttrAxis &GetAttrY() const { return fAttrY; }
   RFrame &SetAttrY(const RAttrAxis &axis) { fAttrY = axis; return *this; }
   RAttrAxis &AttrY() { return fAttrY; }

   const RAttrAxis &GetAttrZ() const { return fAttrZ; }
   RFrame &SetAttrZ(const RAttrAxis &axis) { fAttrZ = axis; return *this; }
   RAttrAxis &AttrZ() { return fAttrZ; }

   RFrame &SetGridX(bool on = true) { fGridX = on; return *this; }
   bool GetGridX() const { return fGridX; }

   RFrame &SetGridY(bool on = true) { fGridY = on; return *this; }
   bool GetGridY() const { return fGridY; }

   RFrame &SetSwapX(bool on = true) { fSwapX = on; return *this; }
   bool GetSwapX() const { return fSwapX; }

   RFrame &SetSwapY(bool on = true) { fSwapY = on; return *this; }
   bool GetSwapY() const { return fSwapY; }

   /** Configure X ticks drawing 0 - off, 1 - as configured for axis, 2 - both sides, 3 - labels on both side */
   RFrame &SetTicksX(int v = 1) { fTicksX = v; return *this; }
   int GetTicksX() const { return fTicksX; }

   /** Configure Y ticks drawing 0 - off, 1 - as configured for axis, 2 - both sides, 3 - labels on both side */
   RFrame &SetTicksY(int v = 1) { fTicksY = v; return *this; }
   int GetTicksY() const { return fTicksY; }

   void GetClientRanges(unsigned connid, RUserRanges &ranges);

   /// Create `nDimensions` default axes for the user coordinate system.
   void GrowToDimensions(size_t nDimensions);

   /// Get the number of axes.
   size_t GetNDimensions() const { return fUserCoord.size(); }

   /// Get the current user coordinate system for a given dimension.
   RPadUserAxisBase &GetUserAxis(size_t dimension) const { return *fUserCoord[dimension]; }

   /// Set the user coordinate system.
   void SetUserAxis(std::vector<std::unique_ptr<RPadUserAxisBase>> &&axes) { fUserCoord = std::move(axes); }

   /// Convert user coordinates to normal coordinates.
   std::array<RPadLength::Normal, 2> UserToNormal(const std::array<RPadLength::User, 2> &pos) const
   {
      return {{fUserCoord[0]->ToNormal(pos[0]), fUserCoord[1]->ToNormal(pos[1])}};
   }
};


} // namespace Experimental
} // namespace ROOT

#endif
