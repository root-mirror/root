// @(#)root/roostats:$Id$
// Author: George Lewis, Kyle Cranmer
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef HISTFACTORY_ASIMOV_H
#define HISTFACTORY_ASIMOV_H

#include <string>
#include <map>

#include "RooWorkspace.h"

namespace RooStats{
  namespace HistFactory {

    class Asimov {

    public:
      
      Asimov() {;}
      Asimov(std::string Name) : fName(Name) {;}

      void ConfigureWorkspace( RooWorkspace* );

      std::string GetName() { return fName; }
      void SetName(const std::string& name) { fName = name; }

      void SetFixedParam(const std::string& param, bool constant=true) { fParamsToFix[param] = constant; }
      void SetParamValue(const std::string& param, double value) { fParamValsToSet[param] = value; }
      
      std::map< std::string, bool >& GetParamsToFix() { return fParamsToFix; }
      std::map< std::string, double >& GetParamsToSet() { return fParamValsToSet; }

    protected:

      std::string fName;

      std::map<std::string, bool> fParamsToFix;
      std::map< std::string, double > fParamValsToSet;

    };


  } // namespace HistFactory
} // namespace RooStats

#endif
