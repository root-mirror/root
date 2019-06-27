/// \file ROOT/RFieldVisitor.hxx
/// \ingroup NTuple ROOT7
/// \author Simon Leisibach <simon.satoshi.rene.leisibach@cern.ch>
/// \date 2019-06-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RFieldVisitor
#define ROOT7_RFieldVisitor

#include <iostream>
#include <sstream>
#include <string>


namespace ROOT {
namespace Experimental {
    class RFieldRoot;
namespace Detail {
    class RFieldBase;
    }
    
    
class RNTupleVisitor {
public:
    virtual void visitField(const ROOT::Experimental::Detail::RFieldBase& fField) = 0;
};
    
class RPrintVisitor: public RNTupleVisitor {
private:
    /// holds ostream which is used for output
    std::ostream &fOutput;
    int fWidth;
    int maxNoFields;
public:
    RPrintVisitor(std::ostream& out = std::cout, int width = 69, int fNoFields = 1000): fOutput{out}, fWidth{width}, maxNoFields{fNoFields} { }
    void visitField(const ROOT::Experimental::Detail::RFieldBase& fField);
};
    
int NumDigits(int x);
int FieldDistance(unsigned int fNoFields);
std::string CutIfNecessary(const std::string &fToCut, unsigned int maxAvailableSpace);
} // namespace Experimental
} // namespace ROOT

#endif
