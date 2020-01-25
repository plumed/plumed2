/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "h36.h"
#include <vector>
#include "Exception.h"

namespace PLMD {

/// Tiny namespace for hybrid36 format.
/// This namespace includes freely available tools for h36 format.
namespace h36 {


/*! C port of the hy36encode() and hy36decode() functions in the
    hybrid_36.py Python prototype/reference implementation.
    See the Python script for more information.

    This file has no external dependencies, NOT even standard C headers.
    Optionally, use hybrid_36_c.h, or simply copy the declarations
    into your code.

    This file is unrestricted Open Source (cctbx.sf.net).
    Please send corrections and enhancements to cctbx@cci.lbl.gov .

    See also: http://cci.lbl.gov/hybrid_36/

    Ralf W. Grosse-Kunstleve, Feb 2007.
 */

/* The following #include may be commented out.
   It is here only to enforce consistency of the declarations
   and the definitions.
 */
// #include <iotbx/pdb/hybrid_36_c.h>

/* All static functions below are implementation details
   (and not accessible from other translation units).
 */

static
const char*
digits_upper() { return "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"; }

static
const char*
digits_lower() { return "0123456789abcdefghijklmnopqrstuvwxyz"; }

static
const char*
value_out_of_range() { return "value out of range."; }

static
const char* invalid_number_literal() { return "invalid number literal."; }

static
const char* unsupported_width() { return "unsupported width."; }

static
void
fill_with_stars(unsigned width, char* result)
{
  while (width) {
    *result++ = '*';
    width--;
  }
  *result = '\0';
}

static
void
encode_pure(
  const char* digits,
  unsigned digits_size,
  unsigned width,
  int value,
  char* result)
{
  char buf[16];
  int rest;
  unsigned i, j;
  i = 0;
  j = 0;
  if (value < 0) {
    j = 1;
    value = -value;
  }
  while (1) {
    rest = value / digits_size;
    buf[i++] = digits[value - rest * digits_size];
    if (rest == 0) break;
    value = rest;
  }
  if (j) buf[i++] = '-';
  for(j=i; j<width; j++) *result++ = ' ';
  while (i != 0) *result++ = buf[--i];
  *result = '\0';
}

static
const char*
decode_pure(
  const int* digits_values,
  unsigned digits_size,
  const char* s,
  unsigned s_size,
  int* result)
{
  int si, dv;
  int have_minus = 0;
  int have_non_blank = 0;
  int value = 0;
  unsigned i = 0;
  for(; i<s_size; i++) {
    si = s[i];
    if (si < 0 || si > 127) {
      *result = 0;
      return invalid_number_literal();
    }
    if (si == ' ') {
      if (!have_non_blank) continue;
      value *= digits_size;
    }
    else if (si == '-') {
      if (have_non_blank) {
        *result = 0;
        return invalid_number_literal();
      }
      have_non_blank = 1;
      have_minus = 1;
      continue;
    }
    else {
      have_non_blank = 1;
      dv = digits_values[si];
      if (dv < 0 || dv >= digits_size) {
        *result = 0;
        return invalid_number_literal();
      }
      value *= digits_size;
      value += dv;
    }
  }
  if (have_minus) value = -value;
  *result = value;
  return 0;
}

/*! hybrid-36 encoder: converts integer value to string result

      width: must be 4 (e.g. for residue sequence numbers)
                  or 5 (e.g. for atom serial numbers)

      value: integer value to be converted

      result: pointer to char array of size width+1 or greater
              on return result is null-terminated

      return value: pointer to error message, if any,
                    or 0 on success

    Example usage (from C++):
      char result[4+1];
      const char* errmsg = hy36encode(4, 12345, result);
      if (errmsg) throw std::runtime_error(errmsg);
 */
const char*
hy36encode(unsigned width, int value, char* result)
{
  int i = value;
  if (width == 4U) {
    if (i >= -999) {
      if (i < 10000) {
        encode_pure(digits_upper(), 10U, 4U, i, result);
        return 0;
      }
      i -= 10000;
      if (i < 1213056 /* 26*36**3 */) {
        i += 466560 /* 10*36**3 */;
        encode_pure(digits_upper(), 36U, 0U, i, result);
        return 0;
      }
      i -= 1213056;
      if (i < 1213056) {
        i += 466560;
        encode_pure(digits_lower(), 36U, 0U, i, result);
        return 0;
      }
    }
  }
  else if (width == 5U) {
    if (i >= -9999) {
      if (i < 100000) {
        encode_pure(digits_upper(), 10U, 5U, i, result);
        return 0;
      }
      i -= 100000;
      if (i < 43670016 /* 26*36**4 */) {
        i += 16796160 /* 10*36**4 */;
        encode_pure(digits_upper(), 36U, 0U, i, result);
        return 0;
      }
      i -= 43670016;
      if (i < 43670016) {
        i += 16796160;
        encode_pure(digits_lower(), 36U, 0U, i, result);
        return 0;
      }
    }
  }
  else {
    fill_with_stars(width, result);
    return unsupported_width();
  }
  fill_with_stars(width, result);
  return value_out_of_range();
}

/*! hybrid-36 decoder: converts string s to integer result

      width: must be 4 (e.g. for residue sequence numbers)
                  or 5 (e.g. for atom serial numbers)

      s: string to be converted
         does not have to be null-terminated

      s_size: size of s
              must be equal to width, or an error message is
              returned otherwise

      result: integer holding the conversion result

      return value: pointer to error message, if any,
                    or 0 on success

    Example usage (from C++):
      int result;
      const char* errmsg = hy36decode(width, "A1T5", 4, &result);
      if (errmsg) throw std::runtime_error(errmsg);
 */
const char*
hy36decode(unsigned width, const char* s, unsigned s_size, int* result)
{
  static const std::vector<int> digits_values_upper_vector([]() {
    std::vector<int> ret(128U,-1);
    for(unsigned i=0; i<36U; i++) {
      int di = digits_upper()[i];
      if (di < 0 || di > 127) {
        plumed_error()<<"internal error hy36decode: integer value out of range";
      }
      ret[di] = i;
    }
    return ret;
  }());
  static const int* digits_values_upper=digits_values_upper_vector.data();
  static const std::vector<int> digits_values_lower_vector([]() {
    std::vector<int> ret(128U,-1);
    for(unsigned i=0; i<36U; i++) {
      int di = digits_lower()[i];
      if (di < 0 || di > 127) {
        plumed_error()<<"internal error hy36decode: integer value out of range";
      }
      ret[di] = i;
    }
    return ret;
  }());
  static const int* digits_values_lower=digits_values_lower_vector.data();
  int di;
  const char* errmsg;
  if (s_size == width) {
    di = s[0];
    if (di >= 0 && di <= 127) {
      if (digits_values_upper[di] >= 10) {
        errmsg = decode_pure(digits_values_upper, 36U, s, s_size, result);
        if (errmsg == 0) {
          /* result - 10*36**(width-1) + 10**width */
          if      (width == 4U) (*result) -= 456560;
          else if (width == 5U) (*result) -= 16696160;
          else {
            *result = 0;
            return unsupported_width();
          }
          return 0;
        }
      }
      else if (digits_values_lower[di] >= 10) {
        errmsg = decode_pure(digits_values_lower, 36U, s, s_size, result);
        if (errmsg == 0) {
          /* result + 16*36**(width-1) + 10**width */
          if      (width == 4U) (*result) += 756496;
          else if (width == 5U) (*result) += 26973856;
          else {
            *result = 0;
            return unsupported_width();
          }
          return 0;
        }
      }
      else {
        errmsg = decode_pure(digits_values_upper, 10U, s, s_size, result);
        if (errmsg) return errmsg;
        if (!(width == 4U || width == 5U)) {
          *result = 0;
          return unsupported_width();
        }
        return 0;
      }
    }
  }
  *result = 0;
  return invalid_number_literal();
}

}

}

