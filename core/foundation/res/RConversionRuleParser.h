// @(#)root/core:$Id$
// author: Lukasz Janyst <ljanyst@cern.ch>

#ifndef R__CONVERSION_RULE_PARSER_H
#define R__CONVERSION_RULE_PARSER_H

#include <list>
#include <map>
#include <string>
#include <ostream>

#include "RConfigure.h"

#include "TSchemaType.h"
#include "DllImport.h"

namespace ROOT
{
   //---------------------------------------------------------------------------
   // Global variables
   //---------------------------------------------------------------------------
   typedef std::map<std::string, std::string> SchemaRuleMap_t;
   typedef std::map<std::string, std::list<SchemaRuleMap_t> > SchemaRuleClassMap_t;
   R__EXTERN SchemaRuleClassMap_t gReadRules;
   R__EXTERN SchemaRuleClassMap_t gReadRawRules;

   typedef std::map<std::string, ROOT::Internal::TSchemaType> MembersTypeMap_t;

   //---------------------------------------------------------------------------
   // Create the data member name-type map
   //---------------------------------------------------------------------------
   //   void CreateNameTypeMap( const clang::CXXRecordDecl &cl, MembersTypeMap_t& members );

   //---------------------------------------------------------------------------
   // Check if given rule contains references to valid data members
   //---------------------------------------------------------------------------
   bool HasValidDataMembers( SchemaRuleMap_t& rule, MembersTypeMap_t& members,
                             std::string& error_string);

   //---------------------------------------------------------------------------
   // Write the conversion function for Read rule
   //---------------------------------------------------------------------------
   void WriteReadRuleFunc( SchemaRuleMap_t& rule, int index,
                           std::string& mappedName,
                           MembersTypeMap_t& members, std::ostream& output );


   //---------------------------------------------------------------------------
   // Write the conversion function for ReadRaw rule
   //---------------------------------------------------------------------------
   void WriteReadRawRuleFunc( SchemaRuleMap_t& rule, int index,
                              std::string& mappedName,
                              MembersTypeMap_t& members, std::ostream& output );

   //---------------------------------------------------------------------------
   // Write schema rules
   //---------------------------------------------------------------------------
   void WriteSchemaList( std::list<SchemaRuleMap_t>& rules,
                         const std::string& listName, std::ostream& output );

   //---------------------------------------------------------------------------
   // Get the list of includes defined in schema rules
   //---------------------------------------------------------------------------
   void GetRuleIncludes( std::list<std::string> &result );

   //---------------------------------------------------------------------------
   // Parse read pragma
   //---------------------------------------------------------------------------
   bool ParseRule(std::string rule, ROOT::Internal::MembersMap_t &result, std::string &error_string );

   //---------------------------------------------------------------------------
   // Parse read pragma
   //---------------------------------------------------------------------------
   void ProcessReadPragma( const char* args, std::string& error_string );

   //---------------------------------------------------------------------------
   // Parse readraw pragma
   //---------------------------------------------------------------------------
   void ProcessReadRawPragma( const char* args, std::string& error_string );
}

#endif // R__CONVERSION_RULE_PARSER_H

