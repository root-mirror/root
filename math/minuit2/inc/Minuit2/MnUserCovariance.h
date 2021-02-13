// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnUserCovariance
#define ROOT_Minuit2_MnUserCovariance

#include "Minuit2/MnConfig.h"
#include <vector>
#include <cassert>

namespace ROOT {

namespace Minuit2 {

/**
   Class containing the covariance matrix data represented as a vector of
   size n*(n+1)/2
   Used to hide internal matrix representation to user
 */
class MnUserCovariance {

public:
   MnUserCovariance() : fData(std::vector<double>()), fNRow(0) {}

   // safe constructor using std::vector
   MnUserCovariance(const std::vector<double> &data, unsigned int nrow) : fData(data), fNRow(nrow)
   {
      assert(data.size() == nrow * (nrow + 1) / 2);
   }

   // unsafe constructor using just a pointer
   MnUserCovariance(const double *data, unsigned int nrow)
      : fData(std::vector<double>(data, data + nrow * (nrow + 1) / 2)), fNRow(nrow)
   {
   }

   MnUserCovariance(unsigned int n) : fData(std::vector<double>(n * (n + 1) / 2, 0.)), fNRow(n) {}

   ~MnUserCovariance() {}

   MnUserCovariance(const MnUserCovariance &cov) : fData(cov.fData), fNRow(cov.fNRow) {}

   MnUserCovariance &operator=(const MnUserCovariance &cov)
   {
      if (this != &cov) {
         fData = cov.fData;
         fNRow = cov.fNRow;
      }
      return *this;
   }

   double operator()(unsigned int row, unsigned int col) const
   {
      assert(row < fNRow && col < fNRow);
      if (row > col)
         return fData[col + row * (row + 1) / 2];
      else
         return fData[row + col * (col + 1) / 2];
   }

   double &operator()(unsigned int row, unsigned int col)
   {
      assert(row < fNRow && col < fNRow);
      if (row > col)
         return fData[col + row * (row + 1) / 2];
      else
         return fData[row + col * (col + 1) / 2];
   }

   void Scale(double f)
   {
      for (unsigned int i = 0; i < fData.size(); i++)
         fData[i] *= f;
   }

   const std::vector<double> &Data() const { return fData; }

   unsigned int Nrow() const { return fNRow; }

   // VC 7.1 warning: conversion from size_t to unsigned int
   unsigned int size() const { return static_cast<unsigned int>(fData.size()); }

private:
   std::vector<double> fData;
   unsigned int fNRow;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnUserCovariance
