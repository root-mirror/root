/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2021, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */

#include <TestStatistics/kahan_sum.h>

namespace RooFit {

std::tuple<double, double> kahan_add(double sum, double additive, double carry)
{
   double y = additive - carry;
   double t = sum + y;
   carry = (t - sum) - y;
   sum = t;

   return {sum, carry};
}

}