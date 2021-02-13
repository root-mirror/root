// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:04:59 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class DataOptions

#ifndef ROOT_Fit_DataOptions
#define ROOT_Fit_DataOptions


namespace ROOT {

   namespace Fit {


//___________________________________________________________________________________
/**
   DataOptions : simple structure holding the options on how the data are filled

   @ingroup FitData
*/
struct DataOptions {


   /**
      Default constructor: use the default options
   */
   DataOptions () :
      fIntegral(false),
      fBinVolume(false),
      fNormBinVolume(false),
      fUseEmpty(false),
      fUseRange(false),
      fErrors1(false),
      fExpErrors(false),
      fCoordErrors(true),
      fAsymErrors(true)
   {}


   bool fIntegral;      /// use integral of bin content instead of bin center (default is false)
   bool fBinVolume;     /// normalize data by the bin volume (it is used in the Poisson likelihood fits)
   bool fNormBinVolume; /// normalize data by a normalized the bin volume (bin volume divided by a reference value)
   bool fUseEmpty;      /// use empty bins (default is false) with a fixed error of 1
   bool fUseRange;      /// use the function range when creating the fit data (default is false)
   bool fErrors1;       /// use all errors equal to 1, i.e. fit without errors (default is false)
   bool fExpErrors;     /// use expected errors from the function and not from the data
   bool fCoordErrors;   /// use errors on the x coordinates when available (default is true)
   bool fAsymErrors;    /// use asymmetric errors in the value when available, selecting them according to the on sign of residual (default is true)


};

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_DataOptions */
