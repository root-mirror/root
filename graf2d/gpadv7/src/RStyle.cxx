/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RStyle.hxx>

#include "ROOT/RAttrBase.hxx" // for GPadLog()
#include <ROOT/RDrawable.hxx>
#include <ROOT/RLogger.hxx>

using namespace std::string_literals;

///////////////////////////////////////////////////////////////////////////////
/// Evaluate attribute value for provided RDrawable

const ROOT::Experimental::RAttrMap::Value_t *ROOT::Experimental::RStyle::Eval(const std::string &field, const RDrawable &drawable) const
{
   for (const auto &block : fBlocks) {
      if (drawable.MatchSelector(block.selector)) {
         auto res = block.map.Find(field);
         if (res)
            return res;
      }
   }

   return nullptr;
}

///////////////////////////////////////////////////////////////////////////////
/// Evaluate attribute value for provided selector - exact match is expected

const ROOT::Experimental::RAttrMap::Value_t *ROOT::Experimental::RStyle::Eval(const std::string &field, const std::string &selector) const
{
   for (const auto &block : fBlocks) {
      if (block.selector == selector) {
         auto res = block.map.Find(field);
         if (res)
            return res;
      }
   }

   return nullptr;
}


///////////////////////////////////////////////////////////////////////////////
/// Parse string with CSS code inside

void ROOT::Experimental::RStyle::Clear()
{
   fBlocks.clear();
}

///////////////////////////////////////////////////////////////////////////////
/// Parse string with CSS code inside
/// All data will be append to existing style records

bool ROOT::Experimental::RStyle::ParseString(const std::string &css_code)
{
   if (css_code.empty())
      return true;

   struct RParser {
      int pos{0};
      int nline{1};
      int linebeg{0};
      int len{0};
      const std::string &css_code;

      RParser(const std::string &_code) : css_code(_code)
      {
         len = css_code.length();
      }

      bool more_data() const { return pos < len; }

      char current() const { return css_code[pos]; }

      void shift() { ++pos; }

      bool check_symbol(bool isfirst = false)
      {
         auto symbol = current();
         if (((symbol >= 'a') && (symbol <= 'z')) ||
             ((symbol >= 'A') && (symbol <= 'Z')) || (symbol == '_')) return true;
         return (!isfirst && (symbol>='0') && (symbol<='9'));
      }

      std::string error_position() const
      {
         std::string res = "\nLine "s + std::to_string(nline) + ": "s;

         int p = linebeg;
         while ((p<len) && (p < linebeg+100) && (css_code[p] != '\n')) ++p;

         return res + css_code.substr(linebeg, p-linebeg);
      }

      bool skip_empty()
      {
         bool skip_until_newline = false, skip_until_endblock = false;

         while (pos < len) {
            if (current() == '\n') {
               skip_until_newline = false;
               linebeg = ++pos;
               ++nline;
               continue;
            }

            if (skip_until_endblock && (current() == '*') && (pos+1 < len) && (css_code[pos+1] == '/')) {
               pos+=2;
               skip_until_endblock = false;
               continue;
            }

            if (skip_until_newline || skip_until_endblock || (current() == ' ') || (current() == '\t')) {
               shift();
               continue;
            }

            if ((current() == '/') && (pos+1 < len)) {
               if (css_code[pos+1] == '/') {
                  pos+=2;
                  skip_until_newline = true;
                  continue;
               } else if (css_code[pos+1] == '*')  {
                  pos+=2;
                  skip_until_endblock = true;
                  continue;
               }
            }

            return true;
         }

         return false;
      }

      std::string scan_identifier(bool selector = false)
      {
         if (pos >= len) return ""s;

         int pos0 = pos;

         // start symbols of selector
         if (selector && ((current() == '.') || (current() == '#'))) shift();

         bool is_first = true;

         while ((pos < len) && check_symbol(is_first)) { shift(); is_first = false; }

         return css_code.substr(pos0, pos-pos0);
      }

      std::string scan_value()
      {
          if (pos >= len) return ""s;

          int pos0 = pos;

          while ((pos < len) && (current() != ';') && current() != '\n') shift();

          if (pos >= len)
             return ""s;

          shift();

          return css_code.substr(pos0, pos - pos0 - 1);
      }

   };

   RParser parser(css_code);

   RStyle newstyle;

   while (parser.more_data()) {

      if (!parser.skip_empty())
         return false;

      auto sel = parser.scan_identifier(true);
      if (sel.empty()) {
         R__LOG_ERROR(GPadLog()) << "Fail to find selector" << parser.error_position();
         return false;
      }

      if (!parser.skip_empty())
         return false;

      if (parser.current() != '{') {
         R__LOG_ERROR(GPadLog()) << "Fail to find starting {" << parser.error_position();
         return false;
      }

      parser.shift();

      if (!parser.skip_empty())
         return false;

      auto &map = newstyle.AddBlock(sel);

      while (parser.current() != '}') {
         auto name = parser.scan_identifier();
         if (name.empty()) {
            R__LOG_ERROR(GPadLog()) << "not able to extract identifier" << parser.error_position();
            return false;
         }

         if (!parser.skip_empty())
            return false;

         if (parser.current() != ':') {
            R__LOG_ERROR(GPadLog()) << "not able to find separator :" << parser.error_position();
            return false;
         }

         parser.shift();

         if (!parser.skip_empty())
            return false;

         if (parser.current() == ';') {
            parser.shift();
            map.AddNoValue(name);
         } else {
            auto value = parser.scan_value();
            if (value.empty()) {
               R__LOG_ERROR(GPadLog()) << "not able to find value" << parser.error_position();
               return false;
            }

            map.AddBestMatch(name, value);
         }

         if (!parser.skip_empty())
            return false;
      }

      parser.shift();

      parser.skip_empty(); // after closing } end of file is possible
   }

   // finally move all read blocks to this
   fBlocks.splice(fBlocks.end(), newstyle.fBlocks);

   return true;
}


///////////////////////////////////////////////////////////////////////////////
/// Parse CSS code and returns std::shared_ptr<RStyle> when successful

std::shared_ptr<ROOT::Experimental::RStyle> ROOT::Experimental::RStyle::Parse(const std::string &css_code)
{
   auto style = std::make_shared<RStyle>();
   if (!style->ParseString(css_code)) return nullptr;
   return style;
}
