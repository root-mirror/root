// @(#)root/mathcore:$Id: FitData.cxx 45049 2012-07-13 12:31:59Z mborinsk $
// Author: M. Borinsky


/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/


#include "Fit/FitData.h"

// Implementation file for class FitData

namespace ROOT {

   namespace Fit {
      FitData::FitData(unsigned int maxpoints, unsigned int dim) :
         fWrapped(false),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(dim),
         fpTmpCoordVector(nullptr)
      {
         assert(fDim >= 1);
         InitCoordsVector();
      }

      /// construct passing options and default data range
      FitData::FitData(const DataOptions &opt, unsigned int maxpoints, unsigned int dim) :
         fWrapped(false),
         fOptions(opt),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(dim),
         fpTmpCoordVector(nullptr)
      {
         assert(fDim >= 1);
         InitCoordsVector();
      }


      /// construct passing range and default options
      FitData::FitData(const DataRange &range, unsigned int maxpoints, unsigned int dim) :
         fWrapped(false),
         fRange(range),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(dim),
         fpTmpCoordVector(nullptr)
      {
         assert(fDim >= 1);
         InitCoordsVector();
      }

      /// construct passing options and data range
      FitData::FitData(const DataOptions &opt, const DataRange &range,
                       unsigned int maxpoints, unsigned int dim) :
         fWrapped(false),
         fOptions(opt),
         fRange(range),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(dim),
         fCoords(fDim),
         fCoordsPtr(fDim),
         fpTmpCoordVector(nullptr)
      {
         assert(fDim >= 1);
         InitCoordsVector();
      }

      /// constructor from external data for 1D data
      FitData::FitData(unsigned int n, const double *dataX) :
         fWrapped(true),
         fMaxPoints(n),
         fNPoints(n),
         fDim(1),
         fCoordsPtr(fDim),
         fpTmpCoordVector(nullptr)
      {
         assert(dataX);
         fCoordsPtr[0] = dataX;

         if (fpTmpCoordVector) {
            delete[] fpTmpCoordVector;
            fpTmpCoordVector = nullptr;
         }

         fpTmpCoordVector = new double [fDim];
      }

      /// constructor from external data for 2D data
      FitData::FitData(unsigned int n, const double *dataX, const double *dataY) :
         fWrapped(true),
         fMaxPoints(n),
         fNPoints(n),
         fDim(2),
         fCoordsPtr(fDim),
         fpTmpCoordVector(nullptr)
      {
         assert(dataX && dataY);
         fCoordsPtr[0] = dataX;
         fCoordsPtr[1] = dataY;

         if (fpTmpCoordVector) {
            delete[] fpTmpCoordVector;
            fpTmpCoordVector = nullptr;
         }

         fpTmpCoordVector = new double [fDim];
      }

      /// constructor from external data for 3D data
      FitData::FitData(unsigned int n, const double *dataX, const double *dataY,
                       const double *dataZ) :
         fWrapped(true),
         fMaxPoints(n),
         fNPoints(fMaxPoints),
         fDim(3),
         fCoordsPtr(fDim),
         fpTmpCoordVector(nullptr)
      {
         assert(dataX && dataY && dataZ);
         fCoordsPtr[0] = dataX;
         fCoordsPtr[1] = dataY;
         fCoordsPtr[2] = dataZ;

         if (fpTmpCoordVector) {
            delete[] fpTmpCoordVector;
            fpTmpCoordVector = nullptr;
         }

         fpTmpCoordVector = new double [fDim];
      }

      /**
        constructor for multi-dim external data and a range (data are copied inside according to the range)
        Uses as argument an iterator of a list (or vector) containing the const double * of the data
        An example could be the std::vector<const double *>::begin
      */
      FitData::FitData(const DataRange &range, unsigned int maxpoints, const double *dataX) :
         fWrapped(false),
         fRange(range),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(1),
         fpTmpCoordVector(nullptr)
      {
         InitCoordsVector();

         const double *ptrList[] = { dataX };

         InitFromRange(ptrList);
      }

      /**
        constructor for multi-dim external data and a range (data are copied inside according to the range)
        Uses as argument an iterator of a list (or vector) containing the const double * of the data
        An example could be the std::vector<const double *>::begin
      */
      FitData::FitData(const DataRange &range, unsigned int maxpoints, const double *dataX, const double *dataY) :
         fWrapped(false),
         fRange(range),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(2),
         fpTmpCoordVector(nullptr)
      {
         InitCoordsVector();

         const double *ptrList[] = { dataX, dataY };

         InitFromRange(ptrList);
      }

      /**
        constructor for multi-dim external data and a range (data are copied inside according to the range)
        Uses as argument an iterator of a list (or vector) containing the const double * of the data
        An example could be the std::vector<const double *>::begin
      */
      FitData::FitData(const DataRange &range, unsigned int maxpoints, const double *dataX, const double *dataY,
                       const double *dataZ) :
         fWrapped(false),
         fRange(range),
         fMaxPoints(maxpoints),
         fNPoints(0),
         fDim(3),
         fpTmpCoordVector(nullptr)
      {
         InitCoordsVector();

         const double *ptrList[] = { dataX, dataY, dataZ };

         InitFromRange(ptrList);
      }

      /// dummy virtual destructor
      FitData::~FitData()
      {
         assert(fWrapped == fCoords.empty());
         for (unsigned int i = 0; i < fDim; i++) {
            assert(fWrapped || fCoords[i].empty() || &fCoords[i].front() == fCoordsPtr[i]);
         }
         if (fpTmpCoordVector)  delete[] fpTmpCoordVector;

      }

      FitData::FitData(const FitData &rhs)
         : fWrapped(false), fMaxPoints(0), fNPoints(0), fDim(0),
           fpTmpCoordVector(nullptr)
      {
         *this = rhs;
      }

      FitData &FitData::operator= (const FitData &rhs)
      {
         fWrapped = rhs.fWrapped;
         fOptions = rhs.fOptions;
         fRange = rhs.fRange;
         fDim = rhs.fDim;
         fMaxPoints = rhs.fMaxPoints;

         if (fWrapped) {
            fCoords.clear();

            fCoordsPtr = rhs.fCoordsPtr;
         } else {
            fCoords = rhs.fCoords;

            fCoordsPtr.resize(fDim);

            for (unsigned int i = 0; i < fDim; i++) {
               fCoordsPtr[i] = fCoords[i].empty() ? nullptr : &fCoords[i].front();
            }
         }

         if (fpTmpCoordVector) {
            delete[] fpTmpCoordVector;
            fpTmpCoordVector = nullptr;
         }

         fpTmpCoordVector = new double [fDim];

         return *this;
      }

      void FitData::Append(unsigned int newPoints, unsigned int dim)
      {
         assert(!fWrapped);

         fMaxPoints = fMaxPoints + newPoints;
         fDim = dim;

         InitCoordsVector();
      }

   } // end namespace Fit

} // end namespace ROOT
