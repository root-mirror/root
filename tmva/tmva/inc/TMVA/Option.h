// @(#)root/tmva $Id$   
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Option                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Option container                                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

#ifndef ROOT_TMVA_Option
#define ROOT_TMVA_Option

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Option                                                               //
//                                                                      //
// Class for TMVA-option handling                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <vector>
#include <string>

#include "TObject.h"
#include "TString.h"
#include "TMVA/MsgLogger.h"

namespace TMVA {

   class Configurable;

   class OptionBase : public TObject {

   public:

      friend class Configurable;

      OptionBase( const TString& name, const TString& desc );
      virtual ~OptionBase() {}

      virtual const char* GetName() const { return fNameAllLower.Data(); }
      virtual const char* TheName() const { return fName.Data(); }
      virtual TString     GetValue(Int_t i=-1) const = 0;

      Bool_t IsSet() const { return fIsSet; }
      virtual Bool_t IsArrayOpt() const = 0;
      const TString& Description() const { return fDescription; }
      virtual Bool_t IsPreDefinedVal(const TString&) const = 0;
      virtual Bool_t HasPreDefinedVal() const = 0;
      virtual Int_t  GetArraySize() const = 0;
      virtual Bool_t SetValue( const TString& vs, Int_t i=-1 );

      using TObject::Print;
      virtual void Print( std::ostream&, Int_t levelofdetail=0 ) const = 0;

   private:

      virtual void SetValueLocal(const TString& vs, Int_t i=-1) = 0;         

      const TString fName;         // name of variable 
      TString fNameAllLower;       // name of variable 
      const TString fDescription;  // its description
      Bool_t        fIsSet;        // set by user ?

   protected:

      static MsgLogger& Log();
   protected:
       
       ClassDef(OptionBase,1);
   };
      
   // ---------------------------------------------------------------------------

   template <class T>

      class Option : public OptionBase {
     
   public:

   Option( T& ref, const TString& name, const TString& desc ) : 
      OptionBase(name, desc), fRefPtr(&ref) {}
      virtual ~Option() {}

      // getters
      virtual TString  GetValue( Int_t i=-1 ) const;
      virtual const T& Value   ( Int_t i=-1 ) const;
      virtual Bool_t HasPreDefinedVal() const { return (fPreDefs.size()!=0); }
      virtual Bool_t IsPreDefinedVal( const TString& ) const; 
      virtual Bool_t IsArrayOpt()   const { return kFALSE; }
      virtual Int_t  GetArraySize() const { return 0; }

      // setters
      virtual void AddPreDefVal(const T&);
      using OptionBase::Print;
      virtual void Print       ( std::ostream&, Int_t levelofdetail=0 ) const;
      virtual void PrintPreDefs( std::ostream&, Int_t levelofdetail=0 ) const;

   protected:

      T& Value(Int_t=-1);
  
      virtual void   SetValueLocal( const TString& val, Int_t i=-1 );
      virtual Bool_t IsPreDefinedValLocal( const T& ) const;

      T* fRefPtr;
      std::vector<T> fPreDefs;  // templated vector
   };      

   template<typename T>
      class Option<T*> : public Option<T> {

   public:

   Option( T*& ref, Int_t size, const TString& name, const TString& desc ) : 
      Option<T>(*ref,name, desc), fVRefPtr(&ref), fSize(size) {}
      virtual ~Option() {}

      TString GetValue( Int_t i ) const {
         std::stringstream str;
         str << std::scientific << Value(i);
         return str.str();
      }
      const T& Value( Int_t i ) const { return (*fVRefPtr)[i]; }
      virtual Bool_t IsArrayOpt()   const { return kTRUE; }
      virtual Int_t  GetArraySize() const { return fSize; }
   
      using Option<T>::Print;
      virtual void Print( std::ostream&, Int_t levelofdetail=0 ) const;

      virtual Bool_t SetValue( const TString& val, Int_t i=0 );

      T& Value(Int_t i) { return (*fVRefPtr)[i]; }
      T ** fVRefPtr;
      Int_t fSize;

   };

} // namespace

namespace TMVA {

   //______________________________________________________________________
   template<class T>
      inline const T& TMVA::Option<T>::Value( Int_t ) const {
      return *fRefPtr;
   }

   template<class T>
      inline T& TMVA::Option<T>::Value( Int_t ) {
      return *fRefPtr;
   }

   template<class T>
      inline TString TMVA::Option<T>::GetValue( Int_t ) const {
      std::stringstream str;      
      str << std::scientific << this->Value();
      return str.str();
   }

   template<>
      inline TString TMVA::Option<Bool_t>::GetValue( Int_t ) const {
      return Value() ? "True" : "False";
   }

   template<>
      inline TString TMVA::Option<Bool_t*>::GetValue( Int_t i ) const {
      return Value(i) ? "True" : "False";
   }

   template<class T>
      inline Bool_t TMVA::Option<T>::IsPreDefinedVal( const TString& val ) const 
      {
         // template 
         T tmpVal;
         std::stringstream str(val.Data());
         str >> tmpVal;
         return IsPreDefinedValLocal(tmpVal);
      }

   template<class T>
      inline Bool_t TMVA::Option<T>::IsPreDefinedValLocal(const T& val) const 
      {
         // template
         if (fPreDefs.size()==0) return kTRUE; // if nothing pre-defined then allow everything

         typename std::vector<T>::const_iterator predefIt;
         predefIt = fPreDefs.begin();
         for (;predefIt!=fPreDefs.end(); predefIt++) 
            if ( (*predefIt)==val ) return kTRUE;
      
         return kFALSE;
      }

   template<>
      inline Bool_t TMVA::Option<TString>::IsPreDefinedValLocal( const TString& val ) const 
      {
         // template specialization for Bool_t 
         TString tVal(val);
         tVal.ToLower();
         if (fPreDefs.size()==0) return kFALSE; // if nothing pre-defined then allow everything
         Bool_t foundPreDef = kFALSE;   
         std::vector<TString>::const_iterator predefIt;
         predefIt = fPreDefs.begin();
         for (;predefIt!=fPreDefs.end(); predefIt++) {
            TString s(*predefIt);
            s.ToLower();
            if (s==tVal) { foundPreDef = kTRUE; break; }
         }
         return foundPreDef;
      }

   //______________________________________________________________________
   template<class T>
      inline void TMVA::Option<T>::AddPreDefVal( const T& val ) 
      {
         // template 
         fPreDefs.push_back(val);
      }

   template<>
      inline void TMVA::Option<Bool_t>::AddPreDefVal( const Bool_t& ) 
      {
         // template specialization for Bool_t 
         Log() << kFATAL << "<AddPreDefVal> predefined values for Option<Bool_t> don't make sense" 
               << Endl;
      }

   template<>
      inline void TMVA::Option<Float_t>::AddPreDefVal( const Float_t& ) 
      {
         // template specialization for Float_t 
         Log() << kFATAL << "<AddPreDefVal> predefined values for Option<Float_t> don't make sense" 
               << Endl;
      }

   template<class T>
      inline void TMVA::Option<T>::Print( std::ostream& os, Int_t levelofdetail ) const 
      {
         // template specialization for TString printing
         os << TheName() << ": " << "\"" << GetValue() << "\"" << " [" << Description() << "]";
         this->PrintPreDefs(os,levelofdetail);
      }

   template<class T>
      inline void TMVA::Option<T*>::Print( std::ostream& os, Int_t levelofdetail ) const 
      {
         // template specialization for TString printing
         for (Int_t i=0; i<fSize; i++) {
            if (i==0)
               os << this->TheName() << "[" << i << "]: " << "\"" << this->GetValue(i) << "\"" << " [" << this->Description() << "]";
            else
               os << "    " << this->TheName() << "[" << i << "]: " << "\"" << this->GetValue(i) << "\"";
            if (i!=fSize-1) os << std::endl;
         }
         this->PrintPreDefs(os,levelofdetail);
      }

   //______________________________________________________________________
   template<class T>
      inline void TMVA::Option<T>::PrintPreDefs( std::ostream& os, Int_t levelofdetail ) const 
      {
         // template specialization for TString printing
         if (HasPreDefinedVal() && levelofdetail>0) {
            os << std::endl << "PreDefined - possible values are:" << std::endl;
            typename std::vector<T>::const_iterator predefIt;
            predefIt = fPreDefs.begin();
            for (;predefIt!=fPreDefs.end(); predefIt++) {
               os << "                       ";
               os << "  - " << (*predefIt) << std::endl;
            }
         }
      }
   
   //______________________________________________________________________
   template<class T>
      inline Bool_t TMVA::Option<T*>::SetValue( const TString& val, Int_t ind ) 
      {
         // template 
         if (ind >= fSize) return kFALSE;
         std::stringstream str(val.Data());
         if (ind < 0) {
            str >> Value(0);
            for (Int_t i=1; i<fSize; i++) Value(i) = Value(0);      
         } 
         else {
            str >> Value(ind);
         }
         return kTRUE;
      }

   template<class T>
      inline void TMVA::Option<T>::SetValueLocal( const TString& val, Int_t i ) 
      {
         // template 
         std::stringstream str(val.Data());
         str >> Value(i);
      }

   template<>
      inline void TMVA::Option<TString>::SetValueLocal( const TString& val, Int_t ) 
      {
         // set TString value
         TString valToSet(val);
         if (fPreDefs.size()!=0) {
            TString tVal(val);
            tVal.ToLower();
            std::vector<TString>::const_iterator predefIt;
            predefIt = fPreDefs.begin();
            for (;predefIt!=fPreDefs.end(); predefIt++) {
               TString s(*predefIt);
               s.ToLower();
               if (s==tVal) { valToSet = *predefIt; break; }
            }
         }

         std::stringstream str(valToSet.Data());
         str >> Value(-1);
      }

   template<>
      inline void TMVA::Option<Bool_t>::SetValueLocal( const TString& val, Int_t ) 
      {
         // set Bool_t value
         TString valToSet(val);
         valToSet.ToLower();
         if (valToSet=="1" || valToSet=="true" || valToSet=="ktrue" || valToSet=="t") {
            this->Value() = true;
         }
         else if (valToSet=="0" || valToSet=="false" || valToSet=="kfalse" || valToSet=="f") {
            this->Value() = false;
         }
         else {
            Log() << kFATAL << "<SetValueLocal> value \'" << val 
                  << "\' can not be interpreted as boolean" << Endl;
         }
      }
}
#endif
