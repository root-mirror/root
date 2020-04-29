/*
 * RooNaNPacker.h
 *
 *  Created on: 28.04.2020
 *      Author: shageboeck
 */

#ifndef ROOFIT_ROOFITCORE_INC_ROONANPACKER_H_
#define ROOFIT_ROOFITCORE_INC_ROONANPACKER_H_

#include <TError.h>

#include <limits>
#include <cstdint>
#include <cmath>
#include <cassert>
#include <numeric>
#include <cstring>

/// Union that can pack a float into the unused bits of the mantissa of a
/// NaN double. This can be used to transport information about violation
/// of function definition ranges in RooFit.
/// To separate NaNs with packed floats from regular NaNs, a tag is
/// written into the upper bits of the mantissa.
struct RooNaNPacker {
  double _payload;

  // In double-valued NaNs, on can abuse the lowest 51 bit for payloads.
  // We use this to pack a float into the lowest 32 bits, which leaves
  // 19-bit to include a magic tag to tell NaNs with payload from ordinary
  // NaNs:
  static constexpr uint64_t magicTagMask = 0x3ffff00000000;
  static constexpr uint64_t magicTag     = 0x321ab00000000;

  constexpr RooNaNPacker() :
    _payload(0.) { }

  /// Create NaN with a packed floating point number.
  explicit RooNaNPacker(float value) :
    _payload(packFloatIntoNaN(value)) {  }

  /// Write into the packed floating point number.
  void setPayload(float payload) {
    _payload = packFloatIntoNaN(payload);

    if (!std::isnan(_payload)) {
      // Little endian machine or similar. Just return NaN.
      warn();
      _payload = std::numeric_limits<double>::quiet_NaN();
    }
  }

  /// Accumulate packed float from a NaN.
  void accumulate(double val) {
    *this += unpackNaN(val);
  }

  /// Unpack floats from NaN doubles, and sum packed values.
  template<class It_t>
  static double accumulatePayloads(It_t begin, It_t end) {
    double sum = std::accumulate(begin, end, 0.f, [](float acc, double val) {
      return acc += unpackNaN(val);
    });

    return packFloatIntoNaN(sum);
  }

  /// Add to the packed float.
  RooNaNPacker& operator+=(float val) {
    setPayload(getPayload() + val);
    return *this;
  }

  /// Multiply the packed float.
  RooNaNPacker& operator*=(float val) {
    setPayload(getPayload() * val);
    return *this;
  }

  /// Retrieve packed float. Returns zero if number is not NaN
  /// or if float wasn't packed by this class.
  float getPayload() const {
    return isNaNWithPayload(_payload) ? unpackNaN(_payload) : 0.;
  }

  /// Test if NaN has a float packed into its mantissa.
  bool isNaNWithPayload() const {
    return isNaNWithPayload(_payload);
  }

  /// Test if NaN has a float packed into its mantissa.
  static bool isNaNWithPayload(double val) {
    uint64_t tmp;
    std::memcpy(&tmp, &val, sizeof(uint64_t));
    return std::isnan(val) && (tmp & magicTagMask) == magicTag;
  }

  /// Test if NaN has a float packed into its mantissa.
  template<typename T>
  static T reinterpret(double theDouble) {
    T tmp;
    std::memcpy(&tmp, &theDouble, sizeof(T));
    return tmp;
  }

  /// Pack float into mantissa of a NaN. Adds a tag to the
  /// upper bits of the mantissa.
  static double packFloatIntoNaN(float payload) {
    double result = std::numeric_limits<double>::quiet_NaN();
    uint64_t tmp;
    std::memcpy(&tmp, &result, sizeof(uint64_t));
    tmp |= magicTag;
    std::memcpy(&tmp, &payload, sizeof(float));
    std::memcpy(&result, &tmp, sizeof(uint64_t));
    return result;
  }

  /// If `val` is NaN and a tag is found in the upper bits of the mantissa,
  /// unpack float from the mantissa.
  static float unpackNaN(double val) {
    float tmp;
    std::memcpy(&tmp, &val, sizeof(float));
    return isNaNWithPayload(val) ? tmp : 0.;
  }

  /// Warn that packing doesn't work on little-endian machines.
  static void warn() {
    static bool haveWarned = false;
    if (!haveWarned)
      Warning("RooNaNPacker", "Fast recovery from NaN not yet implemented on little-endian machines."
          " If necessary, request this on https://root.cern");
    haveWarned = true;
  }
};


#endif /* ROOFIT_ROOFITCORE_INC_ROONANPACKER_H_ */
