// mathtext - A TeX/LaTeX compatible rendering library. Copyright (C)
// 2008-2012 Yue Shi Lai <ylai@users.sourceforge.net>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1 of
// the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA

#include "../inc/mathtext.h"
#include <climits>

/////////////////////////////////////////////////////////////////////

namespace mathtext {

   /////////////////////////////////////////////////////////////////

   void math_text_t::field_t::transform_script(void)
   {
      const unsigned long size = _math_list.size();

      if(size < 2)
         return;

      std::vector<item_t>::reverse_iterator last =
      _math_list.rbegin();
      std::vector<item_t>::reverse_iterator second_last =
      last + 1;

      if(last->_type == item_t::TYPE_ATOM &&
         second_last->_type == item_t::TYPE_ATOM &&
         second_last->_atom._type == atom_t::TYPE_ACC &&
         !(last->_atom._superscript.empty() &&
           last->_atom._subscript.empty())) {
            // Rule 12
            atom_t atom = field_t();

            atom._nucleus._math_list.push_back(*second_last);
            atom._nucleus._math_list.push_back(*last);
            atom._superscript =
            atom._nucleus._math_list.back()._atom._superscript;
            atom._subscript =
            atom._nucleus._math_list.back()._atom._subscript;
            atom._nucleus._math_list.back()._atom._superscript =
            field_t();
            atom._nucleus._math_list.back()._atom._subscript =
            field_t();
            _math_list.pop_back();
            _math_list.pop_back();
            _math_list.push_back(item_t(atom));
         }
   }

   void math_text_t::field_t::append(const item_t &item)
   {
      _math_list.push_back(item);
   }

   // FIXME: Check for malformed "..._a^b_c" (instead of just
   // overwriting)
   void math_text_t::field_t::
   append(const field_t &field, const bool superscript,
          const bool subscript)
   {
      if((superscript || subscript) && _math_list.empty())
         _math_list.push_back(item_t(field_t()));

      if(superscript) {
         _math_list.back()._atom._superscript = field;
         transform_script();
      }
      else if(subscript) {
         _math_list.back()._atom._subscript = field;
         transform_script();
      }
      else
         append(item_t(item_t::TYPE_ATOM, atom_t(field)));
   }

   void math_text_t::field_t::
   prepend(const unsigned int type,
           const math_symbol_t &math_symbol)
   {
      _math_list.insert(_math_list.begin(),
                        item_t(type,
                               atom_t(field_t(math_symbol))));
   }

   void math_text_t::field_t::
   append(const unsigned int type, const math_symbol_t &math_symbol,
          const bool superscript, const bool subscript)
   {
      if((superscript || subscript) && _math_list.empty())
         _math_list.push_back(item_t(field_t()));

      if(superscript) {
         _math_list.back()._atom._superscript =
         field_t(math_symbol);
         transform_script();
      }
      else if(subscript) {
         _math_list.back()._atom._subscript =
         field_t(math_symbol);
         transform_script();
      }
      else
         append(item_t(type, atom_t(field_t(math_symbol))));
   }

   math_text_t::field_t::
   field_t(const std::vector<std::string> &str_split,
           const unsigned int default_family)
   : _type(TYPE_MATH_LIST)
   {
      parse_math_list(str_split, default_family);
   }

   math_text_t::field_t::
   field_t(const std::string &str_delimiter_left,
           const std::vector<std::string> &str_split,
           const std::string &str_delimiter_right,
           const unsigned int default_family)
   : _type(TYPE_MATH_LIST)
   {
      parse_math_list(str_split, default_family);

      const math_symbol_t
      symbol_left(str_delimiter_left, default_family);

      prepend(item_t::TYPE_BOUNDARY, symbol_left);

      const math_symbol_t
      symbol_right(str_delimiter_right, default_family);

      append(item_t::TYPE_BOUNDARY, symbol_right, false, false);
   }

   bool math_text_t::field_t::generalized_fraction(void) const
   {
      if(_type == TYPE_MATH_LIST)
         for (std::vector<item_t>::const_iterator iterator = _math_list.begin(); iterator != _math_list.end();
              ++iterator)
            if(iterator->_type ==
               item_t::TYPE_GENERALIZED_FRACTION)
               return true;
      return false;
   }

   void math_text_t::atom_t::classify(void)
   {
      // Only nucleus affects the atom type (Knuth, The TeXbook,
      // 1986, p. 171)
      if(_nucleus._type == field_t::TYPE_MATH_SYMBOL)
         _type = _nucleus._math_symbol._type;
      else if(_nucleus.generalized_fraction())
         _type = atom_t::TYPE_INNER;
      else
         // FIXME: Does TeX flatten compound expressions before
         // classify them?
         _type = TYPE_ORD;
   }

   bool math_text_t::atom_t::is_combining_diacritical(void) const
   {
      return _nucleus._type == field_t::TYPE_MATH_SYMBOL &&
      _nucleus._math_symbol.is_combining_diacritical();
   }

   unsigned int math_text_t::atom_t::
   spacing(const unsigned int left_type,
           const unsigned int right_type, const bool script)
   {
#include "table/mathspacing.h"
      // Since we only handle atom types upto Inner, the upper bound
      // is type <= TYPE_INNER and not type < NTYPE.
      if(left_type == TYPE_UNKNOWN || left_type > TYPE_INNER ||
         right_type == TYPE_UNKNOWN || right_type > TYPE_INNER) {
         // Invalid
         return 0;
      }

      int index = ((left_type - TYPE_ORD) << 3) +
      (right_type - TYPE_ORD);
      int space = spacing_table[index];

      if(space == nvr) {
         // Invalid
         return 0;
      }
      // Interpret the \nonscript sign, which denotes spaces that
      // should be ignored within the script (and scriptscript)
      // style.
      if(space < 0)
         space = script ? 0 : -space;

      return space;
   }

   bool math_text_t::item_t::operator==(const item_t &item) const
   {
      switch(_type) {
         case TYPE_GENERALIZED_FRACTION:
            return item._type == TYPE_GENERALIZED_FRACTION;
            break;
         default:
            return false;
      }
   }

   std::wstring math_text_t::bad_cast(const std::string string)
   {
      std::wstring wstring;

      for (std::string::const_iterator iterator = string.begin(); iterator != string.end(); ++iterator) {
         wstring.push_back(*iterator);
      }

      return wstring;
   }

   std::wstring math_text_t::utf8_cast(const std::string string)
   {
      std::wstring wstring;

      for(std::string::const_iterator iterator = string.begin();
          iterator != string.end();) {
         wchar_t c;

         // Skip over byte ordering marks
         if((unsigned char)(*iterator) == 0xef) {
            ++iterator;
            if((unsigned char)(*iterator) == 0xbb) {
               ++iterator;
               if((unsigned char)(*iterator) == 0xbf) {
                  ++iterator;
               }
            }
         }
         if((*iterator & 0xf0) == 0xf0) {
            c = (*iterator & 0x7) << 18;
            ++iterator;
            if((*iterator & 0xc0) != 0x80) {
               continue;
            }
            c |= (*iterator & 0x3f) << 12;
            ++iterator;
            if((*iterator & 0xc0) != 0x80) {
               continue;
            }
            c |= (*iterator & 0x3f) << 6;
            ++iterator;
            if((*iterator & 0xc0) != 0x80) {
               continue;
            }
            c |= (*iterator & 0x3f);
            ++iterator;
         }
         else if((*iterator & 0xe0) == 0xe0) {
            c = (*iterator & 0xf) << 12;
            ++iterator;
            if((*iterator & 0xc0) != 0x80) {
               continue;
            }
            c |= (*iterator & 0x3f) << 6;
            ++iterator;
            if((*iterator & 0xc0) != 0x80) {
               continue;
            }
            c |= (*iterator & 0x3f);
            ++iterator;
         }
         else if((*iterator & 0xc0) == 0xc0) {
            c = (*iterator & 0x1f) << 6;
            ++iterator;
            if((*iterator & 0xc0) != 0x80) {
               continue;
            }
            c |= (*iterator & 0x3f);
            ++iterator;
         }
         else if((*iterator & 0x80) == 0x0) {
            c = (*iterator & 0x7f);
            ++iterator;
         }
         else {
            ++iterator;
            continue;
         }
         wstring.push_back(c);
      }

      return wstring;
   }

#if 0
   std::wstring math_text_t::gb_18030_cast(
                                           const std::string string)
   {
   }

   std::wstring math_text_t::shift_jis_x_0213_cast(
                                                   const std::string string)
   {
   }

   std::wstring math_text_t::euc_jis_x_0213_cast(
                                                 const std::string string)
   {
   }

   std::wstring math_text_t::ks_x_2901_cast(
                                            const std::string string)
   {
   }
#endif

   // Apply JIS X 4051:2004

   bool math_text_t::well_formed(void) const
   {
      if(_math_list._type != field_t::TYPE_MATH_LIST)
         return false;
      return true;
   }

#if 0
   std::string tex_form(const double x)
   {
      std::string retval;

      switch(std::fpclassify(x)) {
         default:
            return retval;
      }
   }
#endif

}
