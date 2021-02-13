// @(#)root/tmva:$Id$
// Author: Omar Zapata   2016

/*************************************************************************
 * Copyright (C) 2016, Omar Andres Zapata Mesa                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TMVA_OptionMap
#define ROOT_TMVA_OptionMap

#include <sstream>
#include <map>
#include <string>

#include "TMVA/MsgLogger.h"

#include "TObjString.h"

#include "TObjArray.h"


namespace TMVA {

       /**
        * \class TMVA::OptionMap
        * \ingroup TMVA
        *  class to storage options for the differents methods
        */

       class OptionMap
       {
       protected:
           TString fName;
           std::map<TString,TString> fOptMap; //
           TMVA::MsgLogger fLogger;                 //!
           class Binding
           {
           private:
               std::map<TString,TString> &fInternalMap;
               TString fInternalKey;
           public:
               Binding(std::map<TString,TString>  &fmap,TString key):fInternalMap(fmap),fInternalKey(key){}
               Binding(const Binding &obj):fInternalMap(obj.fInternalMap)
               {
                   fInternalKey  = obj.fInternalKey;
               }
               ~Binding(){}
               void SetKey(TString key){fInternalKey=key;}
               TString GetKey(){return fInternalKey;}
               Binding &operator=(const Binding &obj)
               {
                   fInternalMap  = obj.fInternalMap;
                   fInternalKey  = obj.fInternalKey;
                   return *this;
               }

               template<class T> Binding& operator=(const T &value)
               {
                   ParseValue(fInternalMap[fInternalKey],*const_cast<T*>(&value));
                   return *this;
               }

               template<class T> operator T()
               {
                   return GetValue<T>();
               }
               template<class T> T GetValue()
               {
                   T result;
                   ParseValue(fInternalMap[fInternalKey],result,kFALSE);
                   return result;
               }

               template<class T> void  ParseValue(TString &str,T &value,Bool_t input=kTRUE)
               {
                   std::stringstream fStringStream;
                   if(input)
                   {
                       fStringStream<<value;
                       str=fStringStream.str();
                   }else{
                       fStringStream<<str.Data();
                       fStringStream>>value;
                   }

               }


           };
           Binding fBinder;     //!
       public:
           OptionMap(const TString options="",const TString name="Option"):fName(name),fLogger(name.Data()),fBinder(fOptMap,""){
               ParseOption(options);
           }

           OptionMap(const Char_t *options,const TString name="Option"):fName(name),fLogger(name.Data()),fBinder(fOptMap,""){
               ParseOption(options);
           }

           virtual ~OptionMap(){}

           Bool_t IsEmpty(){return fOptMap.empty();}

           Bool_t HasKey(TString key)
           {
               return fOptMap.count( key )==1;
           }

           Binding& operator[](TString key)
           {
               fBinder.SetKey(key);
               return fBinder;
           }

           OptionMap& operator=(TString options)
           {
               ParseOption(options);
               return *this;
           }

           void Print() const
           {
               MsgLogger Log(fLogger);
               for(auto &item:fOptMap)
               {
                   Log<<kINFO<<item.first.Data()<<": "<<item.second.Data()<<Endl;
               }
           }

           template<class T> T GetValue(const TString & key)
           {
               T result;
               fBinder.ParseValue(fOptMap[key],result,kFALSE);
               return result;
           }


           template<class T> T GetValue(const TString & key) const
           {
               T result;
               std::stringstream oss;
               oss<<fOptMap.at(key);
               oss>>result;
               return result;
           }
           void ParseOption(TString options)
           {
               options.ReplaceAll(" ","");
               auto opts=options.Tokenize(":");
               for(auto opt:*opts)
               {
                   TObjString *objstr=(TObjString*)opt;

                   if(objstr->GetString().Contains("="))
                   {
                      auto pair=objstr->String().Tokenize("=");
                      TObjString *key   = (TObjString *)pair->At(0);
                      TObjString *value = (TObjString *)pair->At(1);

                      fOptMap[key->GetString()] = value->GetString();
                   }else{
                      if(objstr->GetString().BeginsWith("!"))
                      {
                          objstr->String().ReplaceAll("!","");
                          fOptMap[objstr->GetString()]=TString("0");
                      }else{
                          fOptMap[objstr->GetString()]=TString("1");
                      }
                   }
               }

           }
           ClassDef(OptionMap,1);
       };

}

#endif
