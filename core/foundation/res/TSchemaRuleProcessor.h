// @(#)root/core:$Id$
// author: Lukasz Janyst <ljanyst@cern.ch>

#ifndef ROOT_TSchemaRuleProcessor
#define ROOT_TSchemaRuleProcessor

// NOTE: #included by libCore and libCling. All symbols must be inline.

#include <cstring>
#include <string>
#include <list>
#include <utility>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "RtypesCore.h"

#include "TSchemaType.h"

namespace ROOT {
namespace Internal {
   class TSchemaRuleProcessor
   {
      public:
         //---------------------------------------------------------------------
         static void SplitList( const std::string& source,
                                std::list<std::string>& result,
                                char delimiter=',')
         {
            // Split the string producing a list of substrings

            std::string::size_type curr;
            std::string::size_type last = 0;
            std::string::size_type size;
            std::string            elem;

            result.clear();

            while( last != source.size() ) {
               curr = source.find( delimiter, last );

               if( curr == std::string::npos ) {
                  curr = source.size()-1;
                  size = curr-last+1;
               }
               else size = curr-last;

               elem = Trim( source.substr( last, size ) );
               if( !elem.empty() )
                  result.push_back( elem );

               last = curr+1;
            }
         }

         static void SplitDeclaration( const std::string& source,
                                       std::list<std::pair<ROOT::Internal::TSchemaType,std::string> >& result)
         {
            // Split a declaration string producing a list of substrings
            // Typically we have:
            //    int mem; SomeType mem2; SomeTmp<const key, const value> mem3;

            std::string::size_type curr;
            std::string::size_type last = 0;
            std::string::size_type size;
            std::string            elem;
            std::string            type;
            std::string            dims;

            result.clear();

            while( last != source.size() ) {
               // Split on semi-colons.
               curr = source.find( ';', last );

               if( curr == std::string::npos ) {
                  curr = source.size()-1;
                  size = curr-last+1;
               }
               else size = curr-last;

               // Extra spaces.
               elem = Trim( source.substr( last, size ) );
               if( !elem.empty() ) {
                  unsigned int level = 0;

                  // Split between the typename and the membername
                  // Take in consideration template names.
                  for(std::string::size_type j=elem.size(); j>0; --j) {
                     std::string::size_type i = j-1;
                     if (elem[i]=='<') { ++level; }
                     else if (elem[i]=='>') { if (level==0) { continue; } ; --level; }
                     else if (level == 0 && isspace(elem[i])) {
                        type = elem.substr( 0, i );
                        // At the first iteration we know we have a space.
                        while( elem[i]=='*' || elem[i]=='&' || isspace(elem[i]) ) {
                           ++i;
                           if (strcmp("const",elem.c_str()+i)==0 && (i+5)>elem.size()
                               && ( elem[i+5]=='*' || elem[i+5]=='&' || isspace(elem[i+5])) ) {
                              i += 5;
                              type += "const ";
                           } else if (elem[i]=='*' || elem[i]=='&') {
                              type += elem[i];
                           }
                        }
                        std::string::size_type endvar = i;
                        while( endvar!=elem.size() && elem[endvar] != '[' ) {
                           ++endvar;
                        }
                        if (endvar != elem.size() ) {
                           dims = Trim( elem.substr(endvar, elem.size()-endvar) );
                        }
                        elem = Trim( elem.substr(i, endvar-i) );
                        break;
                     }
                  }
                  result.push_back( make_pair(ROOT::Internal::TSchemaType(type,dims),elem) );
               }
               last = curr+1;
            }
         }

         //---------------------------------------------------------------------
         static std::string Trim( const std::string& source )
         {
            // Trim the whitespaces at the beginning and at the end of
            // given source string

            std::string::size_type start, end;
            for( start = 0; start < source.size(); ++start) {
               if ( isspace(source[start]) ) {
                  continue;
               } else if ( source[start] == '\\' && (start+1)<source.size() && (source[start+1]=='\n' || source[start+1]=='\r') ) {
                  ++start;
                  continue;
               } else {
                  // Not a white space.
                  break;
               }
            }
            if( start == source.size() )
               return "";
            for( end = source.size()-1; end > start; --end ) {
               if ( (source[end]=='\n' || source[end]=='\r') && end > (start+1) && source[end-1] == '\\' ) {
                  --end;
                  continue;
               } else if ( isspace(source[end]) ) {
                  continue;
               } else {
                  // Not a white space.
                  break;
               }
            }
            return source.substr( start, end-start+1 );
         }

         //---------------------------------------------------------------------
         static bool ProcessVersion( const std::string& source,
                                     std::pair<Int_t, Int_t>& result )
         {
            // Check if a version is specified correctly
            // The result is set the following way:
            //   x  :  first = x   second = x
            //  -x  :  first = -10 second = x
            // x-y  :  first = x   second = y
            // x-   :  first = x   second = 50000
            // if the given string is invalid (false is returned)
            // then the state of the result is undefined

            std::string::size_type hyphenI;
            std::string            first;
            std::string            second;

            std::string version = Trim( source );

            if( version.empty() )
               return false;

            //------------------------------------------------------------------
            // Do we have a star?
            //------------------------------------------------------------------
            if( version == "*" ) {
               result.first  = -10;
               result.second = 50000;
               return true;
            }

            //------------------------------------------------------------------
            // Check if we have a minus somewhere, if not then single version
            // number was specified
            //------------------------------------------------------------------
            hyphenI = version.find( '-' );
            if( hyphenI == std::string::npos && IsANumber( version ) ) {
               result.first = result.second = atoi( version.c_str() );
               return true;
            }

            //------------------------------------------------------------------
            // We start with the hyphen
            //------------------------------------------------------------------
            if( hyphenI == 0 ) {
               second = Trim( version.substr( 1 ) );
               if( IsANumber( second ) ) {
                  result.first  = -10;
                  result.second = atoi( second.c_str() );
                  return true;
               }
            }

            //------------------------------------------------------------------
            // We end with the hyphen
            //------------------------------------------------------------------
            if( hyphenI == version.size()-1 ) {
               first = Trim( version.substr( 0, version.size()-1 ) );
               if( IsANumber( first ) ) {
                  result.first  = atoi( first.c_str() );
                  result.second = 50000;
                  return true;
               }
            }

            //------------------------------------------------------------------
            // We have the hyphen somewhere in the middle
            //------------------------------------------------------------------
            first  = Trim( version.substr( 0, hyphenI ) );
            second = Trim( version.substr( hyphenI+1, version.size()-hyphenI-1 ) );
            if( IsANumber( first ) && IsANumber( second ) ) {
               result.first  = atoi( first.c_str() );
               result.second = atoi( second.c_str() );
               return true;
            }

            return false;
         }

         //---------------------------------------------------------------------
         /// Check if given string consists of digits.
         static bool IsANumber( const std::string& source, bool acceptHex = false )
         {
            if( source.empty() )
               return false;

            if( acceptHex && source.size() > 2 && source[0] == '0' && source[1] == 'x'
                   && std::all_of(source.begin()+2, source.end(), [](unsigned char c){return std::isxdigit(c);}) )
               return true;

            return std::all_of(source.begin(), source.end(), [](unsigned char c){return std::isdigit(c);});
         }
   };
} // namespace Internal
} // namespace ROOT

#endif // ROOT_TSchemaRuleProcessor
