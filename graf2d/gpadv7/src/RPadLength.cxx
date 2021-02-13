/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RPadLength.hxx"

#include "ROOT/RAttrBase.hxx" // for GPadLog()
#include "ROOT/RLogger.hxx"

#include <string>

///////////////////////////////////////////////////////////////////////
/// Converts RPadLength to string like "0.1 + 25px"
/// User coordinates not (yet) supported

std::string ROOT::Experimental::RPadLength::AsString() const
{
   std::string res;

   if (HasNormal()) {
      double v = GetNormal();
      if (v) res = std::to_string(v);
   }

   if (HasPixel()) {
      double v = GetPixel();
      if ((v > 0) && !res.empty())
         res += " + ";

      if ((v != 0) || res.empty()) {
         res += std::to_string(v);
         res += "px";
      }
   }

   if (!Empty() && res.empty())
      res = "0";

   return res;
}

///////////////////////////////////////////////////////////////////////
/// Parse string and fill RPadLength attributes
/// String can be like "0.1 + 25px"
/// User coordinates not (yet) supported

bool ROOT::Experimental::RPadLength::ParseString(const std::string &value)
{
   Clear();

   if (value.empty())
      return true;

   if (value == "0") {
      SetNormal(0);
      return true;
   }

   if (value == "0px") {
      SetPixel(0);
      return true;
   }

   std::size_t pos = 0;
   std::string val = value;
   int operand = 0;

   while (!val.empty()) {
      // skip empty spaces
      while ((pos < val.length()) && ((val[pos] == ' ') || (val[pos] == '\t')))
         ++pos;

      if (pos >= val.length())
         break;

      if ((val[pos] == '-') || (val[pos] == '+')) {
         if (operand) {
            Clear();
            R__LOG_ERROR(GPadLog()) << "Fail to parse RPadLength " << value;
            return false;
         }
         operand = (val[pos] == '-') ? -1 : 1;
         pos++;
         continue;
      }

      if (pos > 0) {
         val.erase(0, pos);
         pos = 0;
      }

      double v = 0;

      try {
         v = std::stod(val, &pos);
      } catch (...) {
         Clear();
         return false;
      }

      val.erase(0, pos);
      pos = 0;
      if (!operand) operand = 1;

      if ((val.length() > 0) && (val[0] == '%')) {
         val.erase(0, 1);
         SetNormal(GetNormal() + operand*0.01*v);
      } else if ((val.length() > 1) && (val[0] == 'p') && (val[1] == 'x')) {
         val.erase(0, 2);
         SetPixel(GetPixel() + operand*v);
      } else {
         SetNormal(GetNormal() + operand*v);
      }

      operand = 0;
   }

   return true;
}
