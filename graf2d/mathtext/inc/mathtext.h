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

#ifndef MATHTEXT_H_
#define MATHTEXT_H_

#include <vector>
#include <string>
#include <iostream>

namespace mathtext {

   class surface_t;

   /**
    * Hierarchical representation of mathematical formulae
    *
    * The class math_text_t is a hierarchical representation of
    * mathematical formulae similar to TeX's math list.
    *
    * Limitations:
    *
    * - Only (spacing) symbols and combining diacritical marks from
    * TeX, AMS-TeX, LaTeX, AMS-LaTeX, and MathTime Professional II
    * that are representable within ISO/IEC 10646:2003/Amd.2:2006
    * Universal Character Set can be accessed.
    *
    * - Illegal TeX syntax may result in not well defined behaviors.
    * Most notably, `a^b^c' or `a_b_c' result in `a^c' and `a_c',
    * i.e. the overwriting of the previous superscripts and
    * subscripts, and `a \atop b \atop c' or `a \over b \over c'
    * result in `a \atop {b \atop c}' and `a \over {b \over c}', i.e.
    * a right associative interpretation.
    *
    * @see ISO/IEC JTC1/SC2/WG2, ISO/IEC 10646:2003/Amd.2:2006 (ISO,
    * Geneva, 2006).
    * @see D. E. Knuth, The TeXbook (Addision-Wesley, Cambridge, MA,
    * 1986).
    * @see D. E. Knuth, The METAFONTbook (Addision-Wesley, Cambridge,
    * MA, 1986).
    * @see B. Beeton, A. Freytag, M. Sargent III, Unicode support for
    * mathematics, Unicode Technical Report #25
    * @author Yue Shi Lai <ylai@phys.columbia.edu>
    * @version 1.0
    */
   class math_text_t {
   private:
      /**
       * Mathematical symbol
       *
       * The class math_symbol_t represents all (spacing) symbols
       * and combining diacritical marks representable within
       * ISO/IEC 10646:2003/Amd.2:2006 Universal Character Set from
       * TeX, AMS-TeX, LaTeX, AMS-LaTeX, and MathTime Professional
       * II.
       *
       * @author Yue Shi Lai <ylai@phys.columbia.edu>
       * @version 1.0
       */
      class math_symbol_t {
      public:
         enum {
            FAMILY_PLAIN = 0,
            // Math fonts start
            FAMILY_REGULAR,
            FAMILY_ITALIC,
            FAMILY_BOLD,
            FAMILY_BOLD_ITALIC,
            FAMILY_STIX_REGULAR,
            FAMILY_STIX_ITALIC,
            FAMILY_STIX_BOLD,
            FAMILY_STIX_BOLD_ITALIC,
            FAMILY_STIX_SIZE_1_REGULAR,
            FAMILY_STIX_SIZE_1_BOLD,
            FAMILY_STIX_SIZE_2_REGULAR,
            FAMILY_STIX_SIZE_2_BOLD,
            FAMILY_STIX_SIZE_3_REGULAR,
            FAMILY_STIX_SIZE_3_BOLD,
            FAMILY_STIX_SIZE_4_REGULAR,
            FAMILY_STIX_SIZE_4_BOLD,
            FAMILY_STIX_SIZE_5_REGULAR,
            // Below are virtual fonts for bookkeeping and do not
            // have to correspond to a physical font available in
            // surface_t.
            //
            // TeX combined styles
            FAMILY_MATH_ITALIC,
            FAMILY_MATH_BOLD_ITALIC,
            // Mathematical Alphanumerical Symbols, mostly in the
            // Unicode range U+1D500 - U+1D5FF
            FAMILY_MATH_SCRIPT_ITALIC,
            FAMILY_MATH_SCRIPT_BOLD_ITALIC,
            FAMILY_MATH_FRAKTUR_REGULAR,
            FAMILY_MATH_FRAKTUR_BOLD,
            FAMILY_MATH_BLACKBOARD_BOLD,
            FAMILY_MATH_SANS_SERIF_REGULAR,
            FAMILY_MATH_SANS_SERIF_ITALIC,
            FAMILY_MATH_SANS_SERIF_BOLD,
            FAMILY_MATH_SANS_SERIF_BOLD_ITALIC,
            FAMILY_MATH_MONOSPACE,
            NFAMILY
         };
      private:
         /**
          * Sets any family that is mathematical italic into the
          * upright variant of the same weight, thus preventing
          * them from converted into the italic form in the final
          * step.
          *
          * @see void math_italic_is_italic(void)
          */
         void math_italic_is_upright(void);
         /**
          * Sets any family that is mathematical italic into the
          * italic variant of the same weight, usually used in the
          * final step to convert all remaining character in
          * mathematical italic into the italic form.
          *
          * @see void math_italic_is_upright(void)
          */
         void math_italic_is_italic(void);
         void encode_character(void);
         void encode_control_sequence(void);
         void encode_math_blackboard_bold(void);
         void encode_math_script_italic(void);
         void encode_math_script_bold_italic(void);
         void encode_math_fraktur_regular(void);
         void encode_math_fraktur_bold(void);
         void encode_math_sans_serif_regular(void);
         void encode_math_sans_serif_italic(void);
         void encode_math_sans_serif_bold(void);
         void encode_math_sans_serif_bold_italic(void);
         void encode_math_alpha(void);
         void encode(void);
      public:
         std::string _code;
         unsigned int _family;
         wchar_t _glyph;
         unsigned int _type;
         inline math_symbol_t(void)
            : _family(FAMILY_PLAIN), _glyph(L'\0'),
              _type(atom_t::TYPE_UNKNOWN)
         {
         }
         inline math_symbol_t(std::string code, wchar_t glyph,
                         const unsigned int family)
            : _code(code), _family(family), _glyph(glyph),
              _type(atom_t::TYPE_UNKNOWN)
         {
         }
         inline math_symbol_t(std::string code,
                         const unsigned int family)
            : _code(code), _family(family), _glyph(0),
              _type(atom_t::TYPE_UNKNOWN)
         {
            encode();
         }
         inline bool is_combining_diacritical(void) const
         {
            // ISO/IEC JTC1/SC2/WG2 (2003), Annex B, but without
            // the additional characters from the list
            return
               (_glyph >= L'\u0300' && _glyph <= L'\u036f') ||
               (_glyph >= L'\u20d0' && _glyph <= L'\u20ff') ||
               (_glyph >= L'\ufe20' && _glyph <= L'\ufe2f');
         }
         bool bold(void) const;
      };
      /**
       * (Horizontal) Box
       *
       * @author Yue Shi Lai <ylai@phys.columbia.edu>
       * @version 1.0
       */
      class box_t {
      public:
         bool _vertical;
         std::wstring _string;
         box_t(void)
            : _vertical(false), _string(L"")
         {
         }
         box_t(std::wstring string)
            : _vertical(false), _string(string)
         {
         }
      };
      class atom_t;
      class item_t;
      /**
       * Math field
       *
       * @author Yue Shi Lai <ylai@phys.columbia.edu>
       * @version 1.0
       */
      class field_t {
      private:
         void transform_script(void);
         void append(const item_t &item);
         void append(const field_t &field, const bool superscript,
                  const bool subscript);
         void prepend(const unsigned int type,
                   const math_symbol_t &math_symbol);
         void append(const unsigned int type,
                  const math_symbol_t &math_symbol,
                  const bool superscript,
                  const bool subscript);
         void parse_math_list(
            const std::vector<std::string> &str_split,
            const unsigned int default_family);
      public:
         enum {
            // The explicit TYPE_EMPTY in TeX is represented here
            // by an empty math list
            TYPE_UNKNOWN = 0,
            TYPE_MATH_SYMBOL,
            TYPE_BOX,         // unused
            TYPE_MATH_LIST,
            NTYPE
         };
         unsigned int _type;
         math_symbol_t _math_symbol;
         box_t _box;
         std::vector<item_t> _math_list;
         inline field_t(void)
            : _type(TYPE_MATH_LIST)
         {
            // Empty math list == implicit TYPE_EMPTY
         }
         inline field_t(const math_symbol_t &math_symbol)
            : _type(TYPE_MATH_SYMBOL), _math_symbol(math_symbol)
         {
         }
         inline field_t(const box_t &box)
            : _type(TYPE_BOX), _box(box)
         {
         }
         inline field_t(const std::vector<item_t> &math_list)
            : _type(TYPE_MATH_LIST), _math_list(math_list)
         {
         }
         field_t(const std::vector<std::string> &str_split,
               const unsigned int default_family);
         field_t(const std::string &str_delimiter_left,
               const std::vector<std::string> &str_split,
               const std::string &str_delimiter_right,
               const unsigned int default_family);
         inline bool empty(void) const
         {
            return _type == TYPE_MATH_LIST && _math_list.empty();
         }
         bool generalized_fraction(void) const;
      };
      /**
       * Math atom
       *
       * @author Yue Shi Lai <ylai@phys.columbia.edu>
       * @version 1.0
       */
      class atom_t {
      private:
         void classify(void);
      public:
         // TeX's 13 types of math atom
         enum {
            TYPE_UNKNOWN = 0,
            TYPE_ORD,
            TYPE_OP,
            TYPE_BIN,
            TYPE_REL,
            TYPE_OPEN,
            TYPE_CLOSE,
            TYPE_PUNCT,
            TYPE_INNER,
            TYPE_OVER,   // not implemented
            TYPE_UNDER,   // not implemented
            TYPE_ACC,
            TYPE_RAD,
            TYPE_VCENT,   // unused
            NTYPE
         };
         enum {
            LIMITS_UNKNOWN = 0,
            LIMITS_LIMITS,
            LIMITS_NOLIMITS,
            LIMITS_DISPLAYLIMITS,
            NLIMITS
         };
         unsigned int _type;
         // Fields
         field_t _nucleus;
         field_t _superscript;
         field_t _subscript;
         unsigned int _limits;
         // Index the root with the radicand _nucleus, if _type ==
         // TYPE_RAD (while TeX uses absolute positioning for the
         // index, it is not possible to represent it in an font
         // independent way).
         field_t _index;
         inline atom_t(const unsigned int type,
                    const field_t &nucleus)
            : _type(type), _nucleus(nucleus),
              _limits(LIMITS_UNKNOWN)
         {
         }
         inline atom_t(const unsigned int type,
                    const field_t &nucleus,
                    const unsigned int limits)
            : _type(type), _nucleus(nucleus), _limits(limits)
         {
         }
         inline atom_t(const field_t &nucleus)
            : _nucleus(nucleus), _limits(false)
         {
            classify();
         }
         inline atom_t(const field_t &nucleus,
                    const field_t &superscript,
                    const field_t &subscript)
            : _nucleus(nucleus), _superscript(superscript),
              _subscript(subscript), _limits(LIMITS_UNKNOWN)
         {
            classify();
         }
         /**
          * Returns true if the atom is of type Acc, or if it is a
          * math symbol and its character representation is a
          * combining diacritical mark within the ISO/IEC
          * 10646:2003/Amd.2:2006 Universal Character Set, and
          * false otherwise.
          *
          * @return true if the atom is of type Acc, or if it is a
          * math symbol and its character representation is a
          * combining diacritical mark within the ISO/IEC
          * 10646:2003/Amd.2:2006 Universal Character Set, and
          * false otherwise
          */
         bool is_combining_diacritical(void) const;
         /**
          * Returns the interelement spacing between the given left
          * and right atom types, and whether the present style is
          * script or scriptscript, with 0, 1, 2, and 3
          * representing no space, \thinkmuskip, \medmuskip, and
          * \thickmuskip, respectively.
          *
          * TeX defaults to \thinkmuskip to 3 mu, \medmuskip to 4
          * mu, and \thickmuskip to 5 mu, with 1 mu being 1/18
          * quad.
          *
          * @param[in] left_type left (enum) atom type
          * @param[in] right_type right (enum) atom type
          * @return interelement spacing, with 0, 1, 2, and 3
          * representing no space, \thinkmuskip, \medmuskip, and
          * \thickmuskip, respectively
          */
         static unsigned int
         spacing(const unsigned int left_type,
               const unsigned int right_type, const bool script);
      };
      /**
       * Math item
       *
       * @author Yue Shi Lai <ylai@phys.columbia.edu>
       * @version 1.0
       */
      class item_t {
      public:
         // TeX's 9 types of math item
         enum {
            TYPE_UNKNOWN = 0,
            TYPE_ATOM,               // implemented
            TYPE_HORIZONTAL,         // unused
            TYPE_VERTICAL,            // unused
            TYPE_GLOB_OR_GLUE,         // unused
            TYPE_KERN,               // implemented
            TYPE_STYLE_CHANGE,         // not implemented
            TYPE_GENERALIZED_FRACTION,   // not implemented
            TYPE_BOUNDARY,            // implemented
            TYPE_FOUR_WAY_CHOICE,      // unused
            NTYPE
         };
         // TeX's 8 math styles
         enum {
            STYLE_UNKNOWN = 0,
            STYLE_SCRIPT_SCRIPT_PRIME,
            STYLE_SCRIPT_SCRIPT,
            STYLE_SCRIPT_PRIME,
            STYLE_SCRIPT,
            STYLE_TEXT_PRIME,
            STYLE_TEXT,
            STYLE_DISPLAY_PRIME,
            STYLE_DISPLAY,
            NSTYLE
         };
         // TeX's left and right math boundary
         enum {
            BOUNDARY_UNKNOWN = 0,
            BOUNDARY_LEFT,
            BOUNDARY_RIGHT,
            NBOUNDARY
         };
         unsigned int _type;
         atom_t _atom;
         float _length;
         unsigned int _style_change;
         unsigned int _boundary;
         inline item_t(const unsigned int type,
                    const float length = 0)
            : _type(type), _atom(field_t()), _length(length),
              _style_change(STYLE_UNKNOWN),
              _boundary(BOUNDARY_UNKNOWN)
         {
         }
         inline item_t(const atom_t &atom)
            : _type(TYPE_ATOM), _atom(atom), _length(0),
              _style_change(STYLE_UNKNOWN),
              _boundary(BOUNDARY_UNKNOWN)
         {
         }
         inline item_t(const unsigned int type,
                    const atom_t &atom)
            : _type(type), _atom(atom), _length(0),
              _style_change(STYLE_UNKNOWN),
              _boundary(BOUNDARY_UNKNOWN)
         {
         }
         bool operator==(const item_t &item) const;
      };
      std::wstring _code;
      field_t _math_list;
      bool _render_structure;
      /////////////////////////////////////////////////////////////
      void tree_view_prefix(const std::vector<bool> &branch,
                       const bool final) const;
      void tree_view(const field_t &field,
                  std::vector<bool> &branch, const bool final)
         const;
      void tree_view(const item_t &item, std::vector<bool> &branch,
                  const bool final) const;
      void tree_view(const atom_t &atom, std::vector<bool> &branch,
                  const bool final) const;
      static std::wstring bad_cast(const std::string string);
      static std::wstring utf8_cast(const std::string string);
      /////////////////////////////////////////////////////////////
      static std::vector<std::string>
      tex_split(const std::string &raw_code,
              const char escape_character = '\\');
      static std::vector<std::string>
      tex_replace(const std::vector<std::string> &code);
      field_t build_math_list(const std::vector<std::string> &
                        code_split) const
      {
         return field_t(code_split,
                     math_symbol_t::FAMILY_MATH_ITALIC);
      }
   public:
      math_text_t(void)
         : _code(), _math_list(), _render_structure(false)
      {
      }
      math_text_t(const std::string &code_string)
         : _code(bad_cast(code_string)), _render_structure(false)
      {
         std::cerr << __FILE__ << ':' << __LINE__ << ": " << std::endl;
         std::vector<std::string> code_split = tex_split(code_string);
         _math_list = build_math_list(code_split);
      }
      math_text_t(const char code_string[])
         : _code(bad_cast(code_string)), _render_structure(false)
      {
         std::vector<std::string> code_split = tex_split(code_string);
         _math_list = build_math_list(code_split);
      }
      inline std::wstring code(void) const
      {
         return _code;
      }
      inline bool render_structure(void) const
      {
         return _render_structure;
      }
      inline bool &render_structure(void)
      {
         return _render_structure;
      }
      bool well_formed(void) const;
      inline bool empty(void) const
      {
         return _math_list.empty();
      }
      inline void tree_view(void) const
      {
         std::vector<bool> branch;
         tree_view(_math_list, branch, true);
      }
      friend class math_text_renderer_t;
   };

   /**
    * Returns the TeX-formatted scientific representation of a real
    * number
    */
   extern std::string tex_form(const double x);

}

#endif // MATHTEXT_H_
