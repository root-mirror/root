/*
 * ULP difference for use in test comparisons of floats/doubles
 * Based on https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */

#ifndef ROOT_ULPDIFF_H
#define ROOT_ULPDIFF_H

#include <stdint.h> // For int32_t, etc.

union Floaty_t
{
  Floaty_t(float num = 0.0f) : f(num) {}
  // Portable extraction of components.
  bool Negative() const { return i < 0; }

  int32_t i;
  float f;
};

union Doubly_t
{
  Doubly_t(double num = 0.0) : f(num) {}
  // Portable extraction of components.
  bool Negative() const { return i < 0; }

  int64_t i;
  double f;
};


int ulp_diff(float A, float B) {
  Floaty_t uA(A);
  Floaty_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative())
  {
    // Check for equality to make sure +0==-0
    if (A == B)
      return true;
    return false;
  }

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);

  return ulpsDiff;
}

long int ulp_diff(double A, double B) {
  Doubly_t uA(A);
  Doubly_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative())
  {
    // Check for equality to make sure +0==-0
    if (A == B)
      return true;
    return false;
  }

  // Find the difference in ULPs.
  long int ulpsDiff = abs(uA.i - uB.i);

  return ulpsDiff;
}


template <typename F, typename I>
bool AlmostEqualUlps(F A, F B, I maxUlpsDiff) {
  I ulpsDiff = ulp_diff(A, B);
  if (ulpsDiff <= maxUlpsDiff)
    return true;

  return false;
}

#endif //ROOT_ULPDIFF_H
