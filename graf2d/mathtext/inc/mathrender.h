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

#ifndef MATHRENDER_H_
#define MATHRENDER_H_

#if defined(_MSC_VER) && (_MSC_VER < 1800)
// Visual C++ 2008 doesn't have stdint.h
typedef __int8 int8_t;
typedef __int16 int16_t;
typedef __int32 int32_t;
typedef __int64 int64_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif
#include <math.h>
#include <string>
#include <iostream>
#include <algorithm>
#include "mathtext.h"

namespace mathtext {

   /**
    * 2D point (and vector)
    */
   class point_t {
   private:
      float _x[2];
   public:
      inline point_t(void) : _x()
      {
      }
      inline point_t(const point_t &point)
      {
         _x[0] = point._x[0];
         _x[1] = point._x[1];
      }
      inline point_t(const float x0, const float y0)
      {
         _x[0] = x0;
         _x[1] = y0;
      }
      inline point_t& operator=(const point_t &point)
      {
         _x[0] = point._x[0];
         _x[1] = point._x[1];
         return *this;
      }
      inline const float *x(void) const
      {
         return _x;
      }
      inline float *x(void)
      {
         return _x;
      }
      inline float operator[](const int n) const
      {
         return _x[n];
      }
      inline float &operator[](const int n)
      {
         return _x[n];
      }
      inline point_t operator+(const point_t &point) const
      {
         return point_t(_x[0] + point._x[0],
                        _x[1] + point._x[1]);
      }
      inline point_t operator-(const point_t &point) const
      {
         return point_t(_x[0] - point._x[0],
                        _x[1] - point._x[1]);
      }
      inline point_t operator+=(const point_t &point)
      {
         _x[0] += point._x[0];
         _x[1] += point._x[1];

         return *this;
      }
      inline point_t operator-=(const point_t &point)
      {
         _x[0] -= point._x[0];
         _x[1] -= point._x[1];

         return *this;
      }
      inline point_t operator*(const float scale) const
      {
         return point_t(_x[0] * scale, _x[1] * scale);
      }
      inline point_t operator/(const float scale) const
      {
         return point_t(_x[0] / scale, _x[1] / scale);
      }
      inline point_t operator*=(const float scale)
      {
         _x[0] *= scale;
         _x[1] *= scale;

         return *this;
      }
      inline float dot(const point_t &point) const
      {
         return _x[0] * point._x[0] + _x[1] * point._x[1];
      }
      inline float cross(const point_t &point) const
      {
         return _x[0] * point._x[1] - _x[1] * point._x[0];
      }
      inline float norm_square(void) const
      {
         return _x[0] * _x[0] + _x[1] * _x[1];
      }
      inline float norm(void) const
      {
         return sqrtf(norm_square());
      }
      inline point_t unit_vector(void) const
      {
         return *this / norm();
      }
      inline point_t rotate_cw(void) const
      {
         return point_t(_x[1], -_x[0]);
      }
      inline point_t rotate_ccw(void) const
      {
         return point_t(-_x[1], _x[0]);
      }
      inline bool operator==(const point_t &point) const
      {
         return _x[0] == point._x[0] && _x[1] == point._x[1];
      }
      inline bool operator!=(const point_t &point) const
      {
         return _x[0] != point._x[0] || _x[1] != point._x[1];
      }
      friend point_t operator*(const float scale,
                               const point_t &point);
      operator std::string(void) const;
   };

   inline point_t operator*(const float scale, const point_t &point)
   {
      return point_t(scale * point._x[0], scale * point._x[1]);
   }

   /**
    * General 2D affine transform of the Adobe Imaging Model
    *
    * The n-D affine transform is generally represented by
    *
    * ( xt )   ( a0 a1 a2 ) ( x )
    * ( yt ) = ( a3 a4 a5 ) ( y )
    * (  1 )   (  0  0  1 ) ( 1 )
    *
    * with (a0, a1, a3, a4) representing the linear component of the
    * transform, i.e. rotation and scaling, while (a2, a5) is the
    * translation vector
    *
    * @see Adobe Systems, Inc., PostScript Language Reference Manual
    * (Addison-Wesley, Reading, MA, 1999), section 4.3.3, pp.
    * 187-189.
    * @see Adobe Systems, Inc., PDF Reference, 6th Edition, Version
    * 1.7 (Adobe Systems, Inc., San Jose, CA, 2006), section
    * 4.2.2-4.2.3, pp. 204-209.
    */
   class affine_transform_t {
   private:
      float _a[6];
   public:
      static const affine_transform_t identity;
      static const affine_transform_t flip_y;
      static affine_transform_t
      translate(const float tx, const float ty);
      static affine_transform_t
      scale(const float sx, const float sy);
      static affine_transform_t rotate(const float angle);
      inline affine_transform_t(
                                const float a0, const float b, const float c,
                                const float d, const float tx, const float ty)
      {
         _a[0] = a0;
         _a[1] = b;
         _a[2] = c;
         _a[3] = d;
         _a[4] = tx;
         _a[5] = ty;
      }
      inline const float *a(void) const
      {
         return _a;
      }
      inline float *a(void)
      {
         return _a;
      }
      inline float operator[](const int n) const
      {
         return _a[n];
      }
      inline float &operator[](const int n)
      {
         return _a[n];
      }
      inline affine_transform_t linear(void) const
      {
         return affine_transform_t(
                                   _a[0], _a[1], _a[2], _a[3], 0.0F, 0.0F);
      }
      inline affine_transform_t translate(void) const
      {
         return affine_transform_t(
                                   1.0F, 0.0F, 0.0F, 1.0F, _a[4], _a[5]);
      }
      inline affine_transform_t operator*(const float s) const
      {
         return affine_transform_t(
                                   _a[0] * s, _a[1] * s, _a[2] * s, _a[3] * s,
                                   _a[4] * s, _a[5] * s);
      }
      inline affine_transform_t operator/(const float s) const
      {
         return affine_transform_t(
                                   _a[0] / s, _a[1] / s, _a[2] / s, _a[3] / s,
                                   _a[4] / s, _a[5] / s);
      }
      inline point_t operator*(const point_t x) const
      {
         return point_t(
                        _a[0] * x[0] + _a[2] * x[1] + _a[4],
                        _a[1] * x[0] + _a[3] * x[1] + _a[5]);
      }
      /**
       * Returns the affine transform of b, i.e. the matrix-vector
       * product (*this) (b^T, 1)^T, of the vector b
       *
       * @returns the affine transform of the vector b
       */
      inline affine_transform_t
      operator*(const affine_transform_t b) const
      {
         return affine_transform_t(
                                   _a[0] * b._a[0] + _a[1] * b._a[2],
                                   _a[0] * b._a[1] + _a[1] * b._a[3],
                                   _a[2] * b._a[0] + _a[3] * b._a[2],
                                   _a[2] * b._a[1] + _a[3] * b._a[3],
                                   _a[4] * b._a[0] + _a[5] * b._a[2] + b._a[4],
                                   _a[4] * b._a[1] + _a[5] * b._a[3] + b._a[5]);
      }
      /**
       * Returns the determinant det(*this) of the affine transform
       *
       * @returns the determinant det(*this) of the affine transform
       */
      inline float determinant(void) const
      {
         return _a[0] * _a[3] - _a[1] * _a[2];
      }
      /**
       * Returns the determinant (*this)^(-1) of the affine
       * transform
       *
       * @returns the determinant (*this)^(-1) of the affine
       * transform
       */
      inline affine_transform_t inverse(void) const
      {
         return affine_transform_t(
                                   _a[3], -_a[1], -_a[2], _a[0],
                                   _a[2] * _a[5] - _a[3] * _a[4],
                                   _a[1] * _a[4] - _a[0] * _a[5]) *
         (1.0F / determinant());
      }
      operator std::string(void) const;
   };

   /**
    * General TeX bounding box
    */
   // FIXME: The skewchar mechanism is missing
   class bounding_box_t {
   private:
      point_t _lower_left;
      point_t _upper_right;
      float _advance;
      float _italic_correction;
   public:
      inline bounding_box_t(void) : _lower_left(), _upper_right(), _advance(), _italic_correction()
      {
      }
      inline bounding_box_t(
                            const point_t lower_left_0, const point_t upper_right_0,
                            const float advance_0, const float italic_correction_0)
      : _lower_left(lower_left_0), _upper_right(upper_right_0),
      _advance(advance_0),
      _italic_correction(italic_correction_0)
      {
      }
      inline bounding_box_t(
                            const float left_0, const float bottom_0, const float right_0,
                            const float top_0, const float advance_0,
                            const float italic_correction_0)
      : _advance(advance_0),
      _italic_correction(italic_correction_0)
      {
         _lower_left[0] = left_0;
         _lower_left[1] = bottom_0;
         _upper_right[0] = right_0;
         _upper_right[1] = top_0;
      }
      inline point_t lower_left(void) const
      {
         return _lower_left;
      }
      inline point_t &lower_left(void)
      {
         return _lower_left;
      }
      inline point_t upper_right(void) const
      {
         return _upper_right;
      }
      inline point_t &upper_right(void)
      {
         return _upper_right;
      }
      inline float left(void) const
      {
         return _lower_left[0];
      }
      inline float &left(void)
      {
         return _lower_left[0];
      }
      inline float top(void) const
      {
         return _upper_right[1];
      }
      inline float &top(void)
      {
         return _upper_right[1];
      }
      inline float right(void) const
      {
         return _upper_right[0];
      }
      inline float &right(void)
      {
         return _upper_right[0];
      }
      inline float bottom(void) const
      {
         return _lower_left[1];
      }
      inline float &bottom(void)
      {
         return _lower_left[1];
      }
      inline float advance(void) const
      {
         return _advance;
      }
      inline float &advance(void)
      {
         return _advance;
      }
      inline float italic_correction(void) const
      {
         return _italic_correction;
      }
      inline float &italic_correction(void)
      {
         return _italic_correction;
      }
      inline float width(void) const
      {
         return _upper_right[0] - _lower_left[0];
      }
      inline float height(void) const
      {
         return _upper_right[1] - _lower_left[1];
      }
      inline float horizontal_center(void) const
      {
         return 0.5F * (_lower_left[0] + _upper_right[0]);
      }
      inline float vertical_center(void) const
      {
         return 0.5F * (_lower_left[1] + _upper_right[1]);
      }
      inline float ascent(void) const
      {
         return _upper_right[1];
      }
      inline float descent(void) const
      {
         return -_lower_left[1];
      }
      inline bounding_box_t
      merge(const bounding_box_t &bounding_box) const
      {
         bounding_box_t ret;

         ret._lower_left[0] =
         std::min(_lower_left[0],
                  bounding_box._lower_left[0]);
         ret._lower_left[1] =
         std::min(_lower_left[1],
                  bounding_box._lower_left[1]);
         if(bounding_box._upper_right[0] > _upper_right[0]) {
            ret._upper_right[0] = bounding_box._upper_right[0];
            ret._italic_correction =
            bounding_box._italic_correction;
         }
         else {
            ret._upper_right[0] = _upper_right[0];
            ret._italic_correction = _italic_correction;
         }
         ret._upper_right[1] =
         std::max(_upper_right[1],
                  bounding_box._upper_right[1]);
         ret._advance =
         std::max(_upper_right[0] + _advance,
                  bounding_box._upper_right[0] +
                  bounding_box._advance) -
         ret._upper_right[0];

         return ret;
      }
      inline bounding_box_t operator+(const point_t &point) const
      {
         return bounding_box_t(
                               _lower_left + point, _upper_right + point,
                               _advance + point[0], _italic_correction);
      }
      inline bounding_box_t operator-(const point_t &point) const
      {
         return bounding_box_t(
                               _lower_left - point, _upper_right - point,
                               _advance - point[0], _italic_correction);
      }
      inline bounding_box_t operator+=(const point_t &point)
      {
         _lower_left += point;
         _upper_right += point;
         _advance += point[0];

         return *this;
      }
      inline bounding_box_t operator-=(const point_t &point)
      {
         _lower_left -= point;
         _upper_right -= point;
         _advance -= point[0];

         return *this;
      }
      inline bounding_box_t operator*(const float scale) const
      {
         return bounding_box_t(
                               _lower_left * scale, _upper_right * scale,
                               _advance * scale, _italic_correction * scale);
      }
      inline bounding_box_t operator*=(const float scale)
      {
         _lower_left *= scale;
         _upper_right *= scale;
         _advance *= scale;
         _italic_correction *= scale;

         return *this;
      }
      friend bounding_box_t
      operator+(const point_t &, const bounding_box_t &);
      friend bounding_box_t
      operator-(const point_t &, const bounding_box_t &);
      friend bounding_box_t
      operator*(const affine_transform_t &,
                const bounding_box_t &);
   };

   inline bounding_box_t
   operator+(const point_t &point,
             const bounding_box_t &bounding_box)
   {
      return bounding_box_t(
                            point + bounding_box._lower_left,
                            point + bounding_box._upper_right,
                            point[0] + bounding_box._advance,
                            bounding_box._italic_correction);
   }

   inline bounding_box_t
   operator-(const point_t &point,
             const bounding_box_t &bounding_box)
   {
      return bounding_box_t(
                            point - bounding_box._lower_left,
                            point - bounding_box._upper_right,
                            point[0] + bounding_box._advance,
                            bounding_box._italic_correction);
   }

   inline bounding_box_t
   operator*(const affine_transform_t &transform,
             const bounding_box_t &bounding_box)
   {
      return bounding_box_t(
                            transform * bounding_box._lower_left,
                            transform * bounding_box._upper_right,
                            (transform * point_t(bounding_box._advance, 0))[0],
                            (transform * point_t(
                                                 bounding_box._italic_correction, 0))[0]);
   }

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER
   /**
    * Mathematical layout engine for formulae represented by
    * math_text_t
    *
    * The class math_text_renderer_t is a layout engine based on
    * TeX's conversion algorithm from a math list to a horizontal
    * list.
    *
    * @see ISO/IEC JTC1/SC2/WG2, ISO/IEC 10646:2003/Amd.2:2006 (ISO,
    * Geneva, 2006).
    * @see D. E. Knuth, The TeXbook (Addision-Wesley, Cambridge, MA,
    * 1986).
    * @see D. E. Knuth, The METAFONTbook (Addision-Wesley, Cambridge,
    * MA, 1986).
    * @see W. Schmidt, The macro package lucimatx (2005),
    * unpublished.
    * @see W. Schmidt, Using the MathTime Professional II fonts with
    * LaTeX (2006), unpublished.
    * @see B. Beeton, A. Freytag, M. Sargent III, Unicode support for
    * mathematics, Unicode Technical Report #25
    * @author Yue Shi Lai <ylai@phys.columbia.edu>
    * @version 1.0
    */
   class math_text_renderer_t {
   public:
      enum {
         DIRECTION_LEFT_TO_RIGHT = 0,
         DIRECTION_RIGHT_TO_LEFT,
         DIRECTION_TOP_TO_BOTTOM
      };
      enum {
         MATH_STYLE_LATIN = 0,
         MATH_STYLE_MAGHREB
      };
      enum {
         FAMILY_PLAIN = 0,
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
         NFAMILY
      };
   private:
      /////////////////////////////////////////////////////////////
      // Font parameter
      static const float script_ratio;
      static const float script_script_ratio;
      static const float thin_mu_skip;
      static const float med_mu_skip;
      static const float thick_mu_skip;
      static const float delimiter_factor;
      static const float delimiter_shortfall;
      /////////////////////////////////////////////////////////////
      static const float num_1;
      static const float num_2;
      static const float num_3;
      static const float denom_1;
      static const float denom_2;
      static const float sup_1;
      static const float sup_2;
      static const float sup_3;
      static const float sub_1;
      static const float sub_2;
      static const float sup_drop;
      static const float sub_drop;
      static const float delim_1;
      static const float delim_2;
      static const float axis_height;
      static const float default_rule_thickness;
      static const float big_op_spacing_1;
      static const float big_op_spacing_2;
      static const float big_op_spacing_3;
      static const float big_op_spacing_4;
      static const float big_op_spacing_5;
      /////////////////////////////////////////////////////////////
      static const float radical_rule_thickness;
      static const float large_operator_display_scale;
      /////////////////////////////////////////////////////////////
      static const float baselineskip_factor;
      /////////////////////////////////////////////////////////////
      // Token
      class math_token_t {
      public:
         point_t _offset;
         bounding_box_t _bounding_box;
         struct extension_t {
            wchar_t _glyph;
            unsigned int _family;
            float _size;
         };
         union {
            unsigned int _style;
            extension_t _extensible;
         };
         float _delimiter_height;
         inline math_token_t(
                             const bounding_box_t bounding_box,
                             const unsigned int style,
                             const float delimiter_height = 0.0F)
         : _offset(0, 0), _bounding_box(bounding_box),
         _style(style), _delimiter_height(delimiter_height)
         {
         }
         inline math_token_t(
                             const point_t offset,
                             const bounding_box_t bounding_box,
                             const unsigned int style,
                             const float delimiter_height = 0.0F)
         : _offset(offset), _bounding_box(bounding_box),
         _style(style), _delimiter_height(delimiter_height)
         {
         }
         inline math_token_t(
                             const bounding_box_t bounding_box,
                             const wchar_t glyph, const unsigned int family,
                             const float size)
         : _offset(0, 0), _bounding_box(bounding_box),
         _delimiter_height(0.0F)
         {
            _extensible._glyph = glyph;
            _extensible._family = family;
            _extensible._size = size;
         }
         inline math_token_t(
                             const point_t offset,
                             const bounding_box_t bounding_box,
                             const wchar_t glyph, const unsigned int family,
                             const float size)
         : _offset(offset), _bounding_box(bounding_box),
         _delimiter_height(0.0F)
         {
            _extensible._glyph = glyph;
            _extensible._family = family;
            _extensible._size = size;
         }
      };
      /////////////////////////////////////////////////////////////
      // Style test and change
      float style_size(const unsigned int style) const;
      bool is_display_style(const unsigned int style) const;
      bool is_script_style(const unsigned int style) const;
      unsigned int prime_style(const unsigned int style) const;
      bool is_prime_style(const unsigned int style) const;
      template<typename value_t>
      inline value_t
      if_else_display(const unsigned int style,
                      const value_t display_value,
                      const value_t otherwise_value) const
      {
         switch(style) {
            case math_text_t::item_t::STYLE_DISPLAY:
            case math_text_t::item_t::STYLE_DISPLAY_PRIME:
               return display_value;
            default:
               return otherwise_value;
         }
      }
      unsigned int next_superscript_style(const unsigned int style)
      const;
      unsigned int next_subscript_style(const unsigned int style)
      const;
      unsigned int next_numerator_style(const unsigned int style)
      const;
      unsigned int next_denominator_style(const unsigned int style)
      const;
      /////////////////////////////////////////////////////////////
      float x_height(const unsigned int style);
      float quad(const unsigned int style) const;
      unsigned int
      math_family(const math_text_t::math_symbol_t &math_symbol)
      const;
      void post_process_atom_type_initial(unsigned int &atom_type)
      const;
      void post_process_atom_type_interior(
                                           unsigned int &previous_atom_type,
                                           unsigned int &atom_type)
      const;
      bool valid_accent(
                        bool &vertical_alignment,
                        const std::vector<math_text_t::item_t>::const_iterator &
                        iterator,
                        const std::vector<math_text_t::item_t>::const_iterator &
                        math_list_end) const;
      float kerning_mu(float amount) const;
      float math_spacing(
                         unsigned int left_type, unsigned int right_type,
                         unsigned int style) const;
   protected:
      inline virtual affine_transform_t
      transform_logical_to_pixel(void) const = 0;
      virtual affine_transform_t
      transform_pixel_to_logical(void) const = 0;
      /////////////////////////////////////////////////////////////
      // Box rendering
      bounding_box_t math_bounding_box(
                                       const math_text_t::box_t &box, const unsigned int style);
      void math_text(
                     const point_t origin, const math_text_t::box_t &box,
                     const unsigned int style, const bool render_structure);
      /////////////////////////////////////////////////////////////
      // Symbol rendering
      static bool is_wgl_4(const wchar_t c);
      static bool is_left_to_right(const wchar_t c);
      static bool is_right_to_left(const wchar_t c);
      static bool is_top_to_bottom(const wchar_t c);
      static bool is_cyrillic(const wchar_t c);
      static bool is_cjk(const wchar_t c);
      static bool is_cjk_punctuation_open(const wchar_t c);
      static bool is_cjk_punctuation_closed(const wchar_t c);
      bounding_box_t math_bounding_box(
                                       const wchar_t &glyph, const unsigned int family,
                                       const float size);
      void math_text(
                     const point_t origin, const wchar_t &glyph,
                     const unsigned int family, const float size,
                     const bool render_structure);
      bounding_box_t math_bounding_box(
                                       const math_text_t::math_symbol_t &math_symbol,
                                       const unsigned int style);
      void math_text(
                     const point_t origin,
                     const math_text_t::math_symbol_t &math_symbol,
                     const unsigned int style, const bool render_structure);
      /////////////////////////////////////////////////////////////
      // Extensible glyph rendering
      void large_family(
                        unsigned long &nfamily, const unsigned int *&family,
                        const math_text_t::math_symbol_t &math_symbol) const;
      void extensible_glyph(
                            wchar_t glyph[4], unsigned long &nrepeat,
                            const math_text_t::math_symbol_t &math_symbol,
                            const unsigned int style, const float height);
      std::vector<math_token_t> math_tokenize(
                                              const math_text_t::math_symbol_t &math_symbol,
                                              const unsigned int style, const float height);
      bounding_box_t math_bounding_box(
                                       const math_text_t::math_symbol_t &math_symbol,
                                       const unsigned int style, const float height);
      void math_text(
                     const point_t origin,
                     const math_text_t::math_symbol_t &math_symbol,
                     const unsigned int style, const float height,
                     const bool render_structure);
      /////////////////////////////////////////////////////////////
      // Math list rendering
      std::vector<math_token_t> math_tokenize(
                                              const std::vector<math_text_t::item_t>::const_iterator &
                                              math_list_begin,
                                              const std::vector<math_text_t::item_t>::const_iterator &
                                              math_list_end,
                                              unsigned int style);
      bounding_box_t math_bounding_box(
                                       const std::vector<math_text_t::item_t>::const_iterator &
                                       math_list_begin,
                                       const std::vector<math_text_t::item_t>::const_iterator &
                                       math_list_end,
                                       unsigned int style);
      void math_text(
                     const point_t origin,
                     const std::vector<math_text_t::item_t>::const_iterator &
                     math_list_begin,
                     const std::vector<math_text_t::item_t>::const_iterator &
                     math_list_end,
                     const unsigned int style, const bool render_structure);
      /////////////////////////////////////////////////////////////
      // Field rendering
      bounding_box_t math_bounding_box(
                                       const math_text_t::field_t &field,
                                       const unsigned int style);
      void math_text(
                     const point_t origin, const math_text_t::field_t &field,
                     const unsigned int style, const bool render_structure);
      /////////////////////////////////////////////////////////////
      // Atom rendering
      std::vector<math_token_t> math_tokenize(
                                              const math_text_t::atom_t &atom, unsigned int style);
      bounding_box_t math_bounding_box(
                                       const math_text_t::atom_t &atom,
                                       const unsigned int style);
      void math_text(
                     const point_t origin, const math_text_t::atom_t &atom,
                     const unsigned int style, const bool render_structure);
      /////////////////////////////////////////////////////////////
   public:
      /////////////////////////////////////////////////////////////
      // Constructor and destructor
      inline math_text_renderer_t(void)
      {
      }
      inline virtual ~math_text_renderer_t(void)
      {
      }
      /////////////////////////////////////////////////////////////
      // Virtual functions
      virtual float font_size(
                              const unsigned int family = FAMILY_PLAIN) const = 0;
      virtual void set_font_size(
                                 const float size, const unsigned int family) = 0;
      inline virtual void set_font_size(const float size) = 0;
      virtual void reset_font_size(
                                   const unsigned int family) = 0;
      virtual void point(const float x, const float y) = 0;
      virtual void filled_rectangle(
                                    const bounding_box_t &bounding_box) = 0;
      virtual void rectangle(
                             const bounding_box_t &bounding_box) = 0;
      virtual bounding_box_t bounding_box(
                                          const std::wstring string,
                                          const unsigned int family = FAMILY_PLAIN) = 0;
      virtual void text_raw(
                            const float x, const float y, const std::wstring string,
                            const unsigned int family = FAMILY_PLAIN) = 0;
      virtual void text_with_bounding_box(
                                          const float x, const float y, const std::wstring string,
                                          const unsigned int family = FAMILY_PLAIN) = 0;
      /////////////////////////////////////////////////////////////
      // Interface
      bounding_box_t bounding_box(
                                  const math_text_t &math_text,
                                  const bool display_style = false);
      void text(
                const float x, const float y,
                const math_text_t &math_text,
                const bool display_style = false);
      /////////////////////////////////////////////////////////////
      inline float default_axis_height(
                                       const bool /*display_style = false*/) const
      {
         return axis_height * style_size(
                                         math_text_t::item_t::STYLE_TEXT);
      }
      /////////////////////////////////////////////////////////////
   };
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

}

#endif // MATHRENDER_H_
