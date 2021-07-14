// Author: Enrico Guiraud, Danilo Piparo CERN  02/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RStringView.hxx>
#include <ROOT/TSeq.hxx>
#include <RtypesCore.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TClass.h>
#include <TClassEdit.h>
#include <TFriendElement.h>
#include <TInterpreter.h>
#include <TObject.h>
#include <TPRegexp.h>
#include <TString.h>
#include <TTree.h>

// pragma to disable warnings on Rcpp which have
// so many noise compiling
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "lexertk.hpp"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <typeinfo>

namespace ROOT {
namespace Detail {
namespace RDF {
class RDefineBase;
class RFilterBase;
class RLoopManager;
class RRangeBase;
} // namespace RDF
} // namespace Detail

namespace RDF {
class RDataSource;
} // namespace RDF

} // namespace ROOT

namespace {
using ROOT::Detail::RDF::ColumnNames_t;

/// A string expression such as those passed to Filter and Define, digested to a standardized form
struct ParsedExpression {
   /// The string expression with the dummy variable names in fVarNames in place of the original column names
   std::string fExpr;
   /// The list of valid column names that were used in the original string expression.
   /// Duplicates are removed and column aliases (created with Alias calls) are resolved.
   ColumnNames_t fUsedCols;
   /// The list of variable names used in fExpr, with same ordering and size as fUsedCols
   ColumnNames_t fVarNames;
};

static bool IsStrInVec(const std::string &str, const std::vector<std::string> &vec)
{
   return std::find(vec.cbegin(), vec.cend(), str) != vec.cend();
}

// look at expression `expr` and return a list of column names used, including aliases
static ColumnNames_t FindUsedColumns(const std::string &expr, const ColumnNames_t &treeBranchNames,
                                     const ColumnNames_t &customColNames, const ColumnNames_t &dataSourceColNames,
                                     const std::map<std::string, std::string> &aliasMap)
{
   ColumnNames_t usedCols;

   lexertk::generator tokens;
   const auto tokensOk = tokens.process(expr);
   if (!tokensOk) {
      const auto msg = "Failed to tokenize expression:\n" + expr + "\n\nMake sure it is valid C++.";
      throw std::runtime_error(msg);
   }

   // iterate over tokens in expression and fill usedCols, varNames and exprWithVars
   const auto nTokens = tokens.size();
   const auto kSymbol = lexertk::token::e_symbol;
   for (auto i = 0u; i < nTokens; ++i) {
      const auto &tok = tokens[i];
      // lexertk classifies '&' as e_symbol for some reason
      if (tok.type != kSymbol || tok.value == "&" || tok.value == "|") {
         // token is not a potential variable name, skip it
         continue;
      }

      ColumnNames_t potentialColNames({tok.value});

      // if token is the start of a dot chain (a.b.c...), a.b, a.b.c etc. are also potential column names
      auto dotChainKeepsGoing = [&](unsigned int _i) {
         return _i + 2 <= nTokens && tokens[_i + 1].value == "." && tokens[_i + 2].type == kSymbol;
      };
      while (dotChainKeepsGoing(i)) {
         potentialColNames.emplace_back(potentialColNames.back() + "." + tokens[i + 2].value);
         i += 2; // consume the tokens we looked at
      }

      // find the longest potential column name that is an actual column name
      // if it's a new match, also add it to usedCols and update varNames
      // potential columns are sorted by length, so we search from the end
      auto isRDFColumn = [&](const std::string &columnOrAlias) {
         const auto &col = ROOT::Internal::RDF::ResolveAlias(columnOrAlias, aliasMap);
         if (IsStrInVec(col, customColNames) || IsStrInVec(col, treeBranchNames) || IsStrInVec(col, dataSourceColNames))
            return true;
         return false;
      };
      const auto longestRDFColMatch = std::find_if(potentialColNames.crbegin(), potentialColNames.crend(), isRDFColumn);

      if (longestRDFColMatch != potentialColNames.crend() && !IsStrInVec(*longestRDFColMatch, usedCols)) {
         // found a new RDF column in the expression (potentially an alias)
         usedCols.emplace_back(*longestRDFColMatch);
      }
   }

   return usedCols;
}

static ParsedExpression ParseRDFExpression(std::string_view expr, const ColumnNames_t &treeBranchNames,
                                           const ColumnNames_t &customColNames, const ColumnNames_t &dataSourceColNames,
                                           const std::map<std::string, std::string> &aliasMap)
{
   // transform `#var` into `__rdf_sizeof_var`
   TString preProcessedExpr(expr);
   // match #varname at beginning of the sentence or after not-a-word, but exclude preprocessor directives like #ifdef
   TPRegexp colSizeReplacer(
      "(^|\\W)#(?!(ifdef|ifndef|if|else|elif|endif|pragma|define|undef|include|line))([a-zA-Z_][a-zA-Z0-9_]*)");
   colSizeReplacer.Substitute(preProcessedExpr, "$1__rdf_sizeof_$3", "g");

   const auto usedColsAndAliases =
      FindUsedColumns(std::string(preProcessedExpr), treeBranchNames, customColNames, dataSourceColNames, aliasMap);

   auto escapeDots = [](const std::string &s) {
      TString ss(s);
      TPRegexp dot("\\.");
      dot.Substitute(ss, "\\.", "g");
      return std::string(std::move(ss));
   };

   ColumnNames_t varNames;
   ColumnNames_t usedCols;
   // when we are done, exprWithVars willl be the same as preProcessedExpr but column names will be substituted with
   // the dummy variable names in varNames
   TString exprWithVars(preProcessedExpr);
   for (const auto &colOrAlias : usedColsAndAliases) {
      const auto col = ROOT::Internal::RDF::ResolveAlias(colOrAlias, aliasMap);
      unsigned int varIdx; // index of the variable in varName corresponding to col
      if (!IsStrInVec(col, usedCols)) {
         usedCols.emplace_back(col);
         varIdx = varNames.size();
         varNames.emplace_back("var" + std::to_string(varIdx));
      } else {
         // colOrAlias must be an alias that resolves to a column we have already seen.
         // Find back the corresponding varName
         varIdx = std::distance(usedCols.begin(), std::find(usedCols.begin(), usedCols.end(), col));
      }
      TPRegexp replacer("\\b" + escapeDots(colOrAlias) + "\\b"); // watch out: need to replace colOrAlias, not col
      replacer.Substitute(exprWithVars, varNames[varIdx], "g");
   }

   return ParsedExpression{std::string(std::move(exprWithVars)), std::move(usedCols), std::move(varNames)};
}

/// Return the static global map of Filter/Define lambda expressions that have been jitted.
/// It's used to check whether a given expression has already been jitted, and
/// to look up its associated variable name if it is.
/// Keys in the map are the body of the expression, values are the name of the
/// jitted variable that corresponds to that expression. For example, for:
///     auto lambda1 = [] { return 42; };
/// key would be "[] { return 42; }" and value would be "lambda1".
static std::unordered_map<std::string, std::string> &GetJittedExprs() {
   static std::unordered_map<std::string, std::string> jittedExpressions;
   return jittedExpressions;
}

static std::string
BuildLambdaString(const std::string &expr, const ColumnNames_t &vars, const ColumnNames_t &varTypes)
{
   R__DEBUG_ASSERT(vars.size() == varTypes.size());

   TPRegexp re(R"(\breturn\b)");
   const bool hasReturnStmt = re.MatchB(expr);

   static const std::vector<std::string> fundamentalTypes = {
      "int",
      "signed",
      "signed int",
      "Int_t",
      "unsigned",
      "unsigned int",
      "UInt_t",
      "double",
      "Double_t",
      "float",
      "Float_t",
      "char",
      "Char_t",
      "unsigned char",
      "UChar_t",
      "bool",
      "Bool_t",
      "short",
      "short int",
      "Short_t",
      "long",
      "long int",
      "long long int",
      "Long64_t",
      "unsigned long",
      "unsigned long int",
      "ULong64_t",
      "std::size_t",
      "size_t",
      "Ssiz_t"
   };

   std::stringstream ss;
   ss << "[](";
   for (auto i = 0u; i < vars.size(); ++i) {
      std::string fullType;
      const auto &type = varTypes[i];
      if (std::find(fundamentalTypes.begin(), fundamentalTypes.end(), type) != fundamentalTypes.end()) {
         // pass it by const value to help detect common mistakes such as if(x = 3)
         fullType = "const " + type + " ";
      } else {
         // We pass by reference to avoid expensive copies
         // It can't be const reference in general, as users might want/need to call non-const methods on the values
         fullType = type + "& ";
      }
      ss << fullType << vars[i] << ", ";
   }
   if (!vars.empty())
      ss.seekp(-2, ss.cur);

   if (hasReturnStmt)
      ss << "){";
   else
      ss << "){return ";
   ss << expr << "\n;}";

   return ss.str();
}

/// Declare a lambda expression to the interpreter in namespace __rdf, return the name of the jitted lambda.
/// If the lambda expression is already in GetJittedExprs, return the name for the lambda that has already been jitted.
static std::string DeclareLambda(const std::string &expr, const ColumnNames_t &vars, const ColumnNames_t &varTypes)
{
   R__LOCKGUARD(gROOTMutex);

   const auto lambdaExpr = BuildLambdaString(expr, vars, varTypes);
   auto &exprMap = GetJittedExprs();
   const auto exprIt = exprMap.find(lambdaExpr);
   if (exprIt != exprMap.end()) {
      // expression already there
      const auto lambdaName = exprIt->second;
      return lambdaName;
   }

   // new expression
   const auto lambdaBaseName = "lambda" + std::to_string(exprMap.size());
   const auto lambdaFullName = "__rdf::" + lambdaBaseName;

   const auto toDeclare = "namespace __rdf {\nauto " + lambdaBaseName + " = " + lambdaExpr + ";\nusing " +
                          lambdaBaseName + "_ret_t = typename ROOT::TypeTraits::CallableTraits<decltype(" +
                          lambdaBaseName + ")>::ret_type;\n}";
   ROOT::Internal::RDF::InterpreterDeclare(toDeclare.c_str());

   // InterpreterDeclare could throw. If it doesn't, mark the lambda as already jitted
   exprMap.insert({lambdaExpr, lambdaFullName});

   return lambdaFullName;
}

/// Each jitted lambda comes with a lambda_ret_t type alias for its return type.
/// Resolve that alias and return the true type as string.
static std::string RetTypeOfLambda(const std::string &lambdaName)
{
   const auto dt = gROOT->GetType((lambdaName + "_ret_t").c_str());
   R__ASSERT(dt != nullptr);
   const auto type = dt->GetFullTypeName();
   return type;
}

static void GetTopLevelBranchNamesImpl(TTree &t, std::set<std::string> &bNamesReg, ColumnNames_t &bNames,
                                       std::set<TTree *> &analysedTrees, const std::string friendName = "")
{
   if (!analysedTrees.insert(&t).second) {
      return;
   }

   auto branches = t.GetListOfBranches();
   if (branches) {
      for (auto branchObj : *branches) {
         const auto name = branchObj->GetName();
         if (bNamesReg.insert(name).second) {
            bNames.emplace_back(name);
         } else if (!friendName.empty()) {
            // If this is a friend and the branch name has already been inserted, it might be because the friend
            // has a branch with the same name as a branch in the main tree. Let's add it as <friendname>.<branchname>.
            // If used for a Snapshot, this name will become <friendname>_<branchname> (with an underscore).
            const auto longName = friendName + "." + name;
            if (bNamesReg.insert(longName).second)
               bNames.emplace_back(longName);
         }
      }
   }

   auto friendTrees = t.GetListOfFriends();

   if (!friendTrees)
      return;

   for (auto friendTreeObj : *friendTrees) {
      auto friendElement = static_cast<TFriendElement *>(friendTreeObj);
      auto friendTree = friendElement->GetTree();
      const std::string frName(friendElement->GetName()); // this gets us the TTree name or the friend alias if any
      GetTopLevelBranchNamesImpl(*friendTree, bNamesReg, bNames, analysedTrees, frName);
   }
}

} // anonymous namespace

namespace ROOT {
namespace Internal {
namespace RDF {

/// Take a list of column names, return that list with entries starting by '#' filtered out.
/// The function throws when filtering out a column this way.
ColumnNames_t FilterArraySizeColNames(const ColumnNames_t &columnNames, const std::string &action)
{
   ColumnNames_t columnListWithoutSizeColumns;
   ColumnNames_t filteredColumns;
   std::copy_if(columnNames.begin(), columnNames.end(), std::back_inserter(columnListWithoutSizeColumns),
                [&](const std::string &name) {
                   if (name[0] == '#') {
                     filteredColumns.emplace_back(name);
                      return false;
                   } else {
                      return true;
                   }
                });

   if (!filteredColumns.empty()) {
      std::string msg = "Column name(s) {";
      for (auto &c : filteredColumns)
         msg += c + ", ";
      msg[msg.size() - 2] = '}';
      msg += "will be ignored. Please go through a valid Alias to " + action + " an array size column";
      throw std::runtime_error(msg);
   }

   return columnListWithoutSizeColumns;
}

std::string ResolveAlias(const std::string &col, const std::map<std::string, std::string> &aliasMap)
{
   const auto it = aliasMap.find(col);
   if (it != aliasMap.end())
      return it->second;

   // #var is an alias for __rdf_sizeof_var
   if (col.size() > 1 && col[0] == '#')
      return "__rdf_sizeof_" + col.substr(1);

   return col;
}

void CheckValidCppVarName(std::string_view var, const std::string &where)
{
   bool isValid = true;

   if (var.empty())
      isValid = false;
   const char firstChar = var[0];

   // first character must be either a letter or an underscore
   auto isALetter = [](char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); };
   const bool isValidFirstChar = firstChar == '_' || isALetter(firstChar);
   if (!isValidFirstChar)
      isValid = false;

   // all characters must be either a letter, an underscore or a number
   auto isANumber = [](char c) { return c >= '0' && c <= '9'; };
   auto isValidTok = [&isALetter, &isANumber](char c) { return c == '_' || isALetter(c) || isANumber(c); };
   for (const char c : var)
      if (!isValidTok(c))
         isValid = false;

   if (!isValid) {
      const auto error =
         "RDataFrame::" + where + ": cannot define column \"" + std::string(var) + "\". Not a valid C++ variable name.";
      throw std::runtime_error(error);
   }
}

///////////////////////////////////////////////////////////////////////////////
/// Get all the top-level branches names, including the ones of the friend trees
ColumnNames_t GetTopLevelBranchNames(TTree &t)
{
   std::set<std::string> bNamesSet;
   ColumnNames_t bNames;
   std::set<TTree *> analysedTrees;
   GetTopLevelBranchNamesImpl(t, bNamesSet, bNames, analysedTrees);
   return bNames;
}

// The set here is used as a registry, the real list, which keeps the order, is
// the one in the vector
class RActionBase;

std::string DemangleTypeIdName(const std::type_info &typeInfo)
{
   int dummy(0);
   char *tn = TClassEdit::DemangleTypeIdName(typeInfo, dummy);
   std::string tname(tn);
   free(tn);
   return tname;
}

ColumnNames_t
ConvertRegexToColumns(const ColumnNames_t &colNames, std::string_view columnNameRegexp, std::string_view callerName)
{
   const auto theRegexSize = columnNameRegexp.size();
   std::string theRegex(columnNameRegexp);

   const auto isEmptyRegex = 0 == theRegexSize;
   // This is to avoid cases where branches called b1, b2, b3 are all matched by expression "b"
   if (theRegexSize > 0 && theRegex[0] != '^')
      theRegex = "^" + theRegex;
   if (theRegexSize > 0 && theRegex[theRegexSize - 1] != '$')
      theRegex = theRegex + "$";

   ColumnNames_t selectedColumns;

   // Since we support gcc48 and it does not provide in its stl std::regex,
   // we need to use TPRegexp
   TPRegexp regexp(theRegex);
   for (auto &&colName : colNames) {
      if ((isEmptyRegex || regexp.MatchB(colName.c_str())) && !IsInternalColumn(colName)) {
         selectedColumns.emplace_back(colName);
      }
   }

   if (selectedColumns.empty()) {
      std::string text(callerName);
      if (columnNameRegexp.empty()) {
         text = ": there is no column available to match.";
      } else {
         text = ": regex \"" + std::string(columnNameRegexp) + "\" did not match any column.";
      }
      throw std::runtime_error(text);
   }
   return selectedColumns;
}

/// Throw if column `definedColView` is already there.
void CheckForRedefinition(const std::string &where, std::string_view definedColView, const ColumnNames_t &customCols,
                          const std::map<std::string, std::string> &aliasMap, const ColumnNames_t &treeColumns,
                          const ColumnNames_t &dataSourceColumns)
{
   const std::string definedCol(definedColView); // convert to std::string
   std::string error;

   const auto aliasColNameIt = aliasMap.find(definedCol);
   const bool isAnAlias = aliasColNameIt != aliasMap.end();
   if (isAnAlias) {
      error = "An alias with that name, pointing to column \"" + aliasColNameIt->second + "\", already exists.";
   }

   if (error.empty()) {
      if (std::find(customCols.begin(), customCols.end(), definedCol) != customCols.end())
         error = "A column with that name has already been Define'd. Use Redefine to force redefinition.";
      // else, check if definedCol is in the list of tree branches. This is a bit better than interrogating the TTree
      // directly because correct usage of GetBranch, FindBranch, GetLeaf and FindLeaf can be tricky; so let's assume we
      // got it right when we collected the list of available branches.
      else if (std::find(treeColumns.begin(), treeColumns.end(), definedCol) != treeColumns.end())
         error =
            "A branch with that name is already present in the input TTree/TChain. Use Redefine to force redefinition.";
      else if (std::find(dataSourceColumns.begin(), dataSourceColumns.end(), definedCol) != dataSourceColumns.end())
         error =
            "A column with that name is already present in the input data source. Use Redefine to force redefinition.";
   }

   if (!error.empty()) {
      error = "RDataFrame::" + where + ": cannot define column \"" + definedCol + "\". " + error;
      throw std::runtime_error(error);
   }
}

/// Throw if column `definedColView` is _not_ already there.
void CheckForDefinition(const std::string &where, std::string_view definedColView, const ColumnNames_t &customCols,
                        const std::map<std::string, std::string> &aliasMap, const ColumnNames_t &treeColumns,
                        const ColumnNames_t &dataSourceColumns)
{
   const std::string definedCol(definedColView); // convert to std::string
   std::string error;

   const auto aliasColNameIt = aliasMap.find(definedCol);
   const bool isAnAlias = aliasColNameIt != aliasMap.end();
   if (isAnAlias) {
      error = "An alias with that name, pointing to column \"" + aliasColNameIt->second +
              "\", already exists. Aliases cannot be Redefined.";
   }

   if (error.empty()) {
      const bool isAlreadyDefined = std::find(customCols.begin(), customCols.end(), definedCol) != customCols.end();
      // check if definedCol is in the list of tree branches. This is a bit better than interrogating the TTree
      // directly because correct usage of GetBranch, FindBranch, GetLeaf and FindLeaf can be tricky; so let's assume we
      // got it right when we collected the list of available branches.
      const bool isABranch = std::find(treeColumns.begin(), treeColumns.end(), definedCol) != treeColumns.end();
      const bool isADSColumn =
         std::find(dataSourceColumns.begin(), dataSourceColumns.end(), definedCol) != dataSourceColumns.end();

      if (!isAlreadyDefined && !isABranch && !isADSColumn)
         error = "No column with that name was found in the dataset. Use Define to create a new column.";
   }

   if (!error.empty()) {
      error = "RDataFrame::" + where + ": cannot redefine column \"" + definedCol + "\". " + error;
      throw std::runtime_error(error);
   }
}

void CheckTypesAndPars(unsigned int nTemplateParams, unsigned int nColumnNames)
{
   if (nTemplateParams != nColumnNames) {
      std::string err_msg = "The number of template parameters specified is ";
      err_msg += std::to_string(nTemplateParams);
      err_msg += " while ";
      err_msg += std::to_string(nColumnNames);
      err_msg += " columns have been specified.";
      throw std::runtime_error(err_msg);
   }
}

/// Choose between local column names or default column names, throw in case of errors.
const ColumnNames_t
SelectColumns(unsigned int nRequiredNames, const ColumnNames_t &names, const ColumnNames_t &defaultNames)
{
   if (names.empty()) {
      // use default column names
      if (defaultNames.size() < nRequiredNames)
         throw std::runtime_error(
            std::to_string(nRequiredNames) + " column name" + (nRequiredNames == 1 ? " is" : "s are") +
            " required but none were provided and the default list has size " + std::to_string(defaultNames.size()));
      // return first nRequiredNames default column names
      return ColumnNames_t(defaultNames.begin(), defaultNames.begin() + nRequiredNames);
   } else {
      // use column names provided by the user to this particular transformation/action
      if (names.size() != nRequiredNames) {
         auto msg = std::to_string(nRequiredNames) + " column name" + (nRequiredNames == 1 ? " is" : "s are") +
                    " required but " + std::to_string(names.size()) + (names.size() == 1 ? " was" : " were") +
                    " provided:";
         for (const auto &name : names)
            msg += " \"" + name + "\",";
         msg.back() = '.';
         throw std::runtime_error(msg);
      }
      return names;
   }
}

ColumnNames_t FindUnknownColumns(const ColumnNames_t &requiredCols, const ColumnNames_t &datasetColumns,
                                 const ColumnNames_t &definedCols, const ColumnNames_t &dataSourceColumns)
{
   ColumnNames_t unknownColumns;
   for (auto &column : requiredCols) {
      const auto isBranch = std::find(datasetColumns.begin(), datasetColumns.end(), column) != datasetColumns.end();
      if (isBranch)
         continue;
      const auto isDefine = std::find(definedCols.begin(), definedCols.end(), column) != definedCols.end();
      if (isDefine)
         continue;
      const auto isDataSourceColumn =
         std::find(dataSourceColumns.begin(), dataSourceColumns.end(), column) != dataSourceColumns.end();
      if (isDataSourceColumn)
         continue;
      unknownColumns.emplace_back(column);
   }
   return unknownColumns;
}

std::vector<std::string> GetFilterNames(const std::shared_ptr<RLoopManager> &loopManager)
{
   return loopManager->GetFiltersNames();
}

ParsedTreePath ParseTreePath(std::string_view fullTreeName)
{
   // split name into directory and treename if needed
   std::string_view dirName = "";
   std::string_view treeName = fullTreeName;
   const auto lastSlash = fullTreeName.rfind('/');
   if (std::string_view::npos != lastSlash) {
      dirName = treeName.substr(0, lastSlash);
      treeName = treeName.substr(lastSlash + 1, treeName.size());
   }
   return {std::string(treeName), std::string(dirName)};
}

std::string PrettyPrintAddr(const void *const addr)
{
   std::stringstream s;
   // Windows-friendly
   s << std::hex << std::showbase << reinterpret_cast<size_t>(addr);
   return s.str();
}

void BookFilterJit(const std::shared_ptr<RJittedFilter> &jittedFilter,
                   std::shared_ptr<RDFDetail::RNodeBase> *prevNodeOnHeap, std::string_view name,
                   std::string_view expression, const std::map<std::string, std::string> &aliasMap,
                   const ColumnNames_t &branches, const RBookedDefines &customCols, TTree *tree, RDataSource *ds)
{
   const auto &dsColumns = ds ? ds->GetColumnNames() : ColumnNames_t{};

   const auto parsedExpr =
      ParseRDFExpression(expression, branches, customCols.GetNames(), dsColumns, aliasMap);
   const auto exprVarTypes =
      GetValidatedArgTypes(parsedExpr.fUsedCols, customCols, tree, ds, "Filter", /*vector2rvec=*/true);
   const auto lambdaName = DeclareLambda(parsedExpr.fExpr, parsedExpr.fVarNames, exprVarTypes);
   const auto type = RetTypeOfLambda(lambdaName);
   if (type != "bool")
      std::runtime_error("Filter: the following expression does not evaluate to bool:\n" + std::string(expression));

   // definesOnHeap is deleted by the jitted call to JitFilterHelper
   ROOT::Internal::RDF::RBookedDefines *definesOnHeap = new ROOT::Internal::RDF::RBookedDefines(customCols);
   const auto definesOnHeapAddr = PrettyPrintAddr(definesOnHeap);
   const auto prevNodeAddr = PrettyPrintAddr(prevNodeOnHeap);

   // Produce code snippet that creates the filter and registers it with the corresponding RJittedFilter
   // Windows requires std::hex << std::showbase << (size_t)pointer to produce notation "0x1234"
   std::stringstream filterInvocation;
   filterInvocation << "ROOT::Internal::RDF::JitFilterHelper(" << lambdaName << ", new const char*["
                    << parsedExpr.fUsedCols.size() << "]{";
   for (const auto &col : parsedExpr.fUsedCols)
      filterInvocation << "\"" << col << "\", ";
   if (!parsedExpr.fUsedCols.empty())
      filterInvocation.seekp(-2, filterInvocation.cur); // remove the last ",
   // lifetime of pointees:
   // - jittedFilter: heap-allocated weak_ptr to the actual jittedFilter that will be deleted by JitFilterHelper
   // - prevNodeOnHeap: heap-allocated shared_ptr to the actual previous node that will be deleted by JitFilterHelper
   // - definesOnHeap: heap-allocated, will be deleted by JitFilterHelper
   filterInvocation << "}, " << parsedExpr.fUsedCols.size() << ", \"" << name << "\", "
                    << "reinterpret_cast<std::weak_ptr<ROOT::Detail::RDF::RJittedFilter>*>("
                    << PrettyPrintAddr(MakeWeakOnHeap(jittedFilter)) << "), "
                    << "reinterpret_cast<std::shared_ptr<ROOT::Detail::RDF::RNodeBase>*>(" << prevNodeAddr << "),"
                    << "reinterpret_cast<ROOT::Internal::RDF::RBookedDefines*>(" << definesOnHeapAddr << ")"
                    << ");\n";

   auto lm = jittedFilter->GetLoopManagerUnchecked();
   lm->ToJitExec(filterInvocation.str());
}

// Jit a Define call
std::shared_ptr<RJittedDefine> BookDefineJit(std::string_view name, std::string_view expression, RLoopManager &lm,
                                             RDataSource *ds, const RBookedDefines &customCols,
                                             const ColumnNames_t &branches,
                                             std::shared_ptr<RNodeBase> *upcastNodeOnHeap)
{
   const auto &aliasMap = lm.GetAliasMap();
   auto *const tree = lm.GetTree();
   const auto &dsColumns = ds ? ds->GetColumnNames() : ColumnNames_t{};

   const auto parsedExpr =
      ParseRDFExpression(expression, branches, customCols.GetNames(), dsColumns, aliasMap);
   const auto exprVarTypes =
      GetValidatedArgTypes(parsedExpr.fUsedCols, customCols, tree, ds, "Define", /*vector2rvec=*/true);
   const auto lambdaName = DeclareLambda(parsedExpr.fExpr, parsedExpr.fVarNames, exprVarTypes);
   const auto type = RetTypeOfLambda(lambdaName);

   auto definesCopy = new RBookedDefines(customCols);
   auto definesAddr = PrettyPrintAddr(definesCopy);
   auto jittedDefine = std::make_shared<RDFDetail::RJittedDefine>(name, type, lm.GetNSlots(), lm.GetDSValuePtrs());

   std::stringstream defineInvocation;
   defineInvocation << "ROOT::Internal::RDF::JitDefineHelper(" << lambdaName << ", new const char*["
                    << parsedExpr.fUsedCols.size() << "]{";
   for (const auto &col : parsedExpr.fUsedCols) {
      defineInvocation << "\"" << col << "\", ";
   }
   if (!parsedExpr.fUsedCols.empty())
      defineInvocation.seekp(-2, defineInvocation.cur); // remove the last ",
   // lifetime of pointees:
   // - lm is the loop manager, and if that goes out of scope jitting does not happen at all (i.e. will always be valid)
   // - jittedDefine: heap-allocated weak_ptr that will be deleted by JitDefineHelper after usage
   // - definesAddr: heap-allocated, will be deleted by JitDefineHelper after usage
   defineInvocation << "}, " << parsedExpr.fUsedCols.size() << ", \"" << name
                    << "\", reinterpret_cast<ROOT::Detail::RDF::RLoopManager*>(" << PrettyPrintAddr(&lm)
                    << "), reinterpret_cast<std::weak_ptr<ROOT::Detail::RDF::RJittedDefine>*>("
                    << PrettyPrintAddr(MakeWeakOnHeap(jittedDefine))
                    << "), reinterpret_cast<ROOT::Internal::RDF::RBookedDefines*>(" << definesAddr
                    << "), reinterpret_cast<std::shared_ptr<ROOT::Detail::RDF::RNodeBase>*>("
                    << PrettyPrintAddr(upcastNodeOnHeap) << "));\n";

   lm.ToJitExec(defineInvocation.str());
   return jittedDefine;
}

// Jit and call something equivalent to "this->BuildAndBook<ColTypes...>(params...)"
// (see comments in the body for actual jitted code)
std::string JitBuildAction(const ColumnNames_t &cols, std::shared_ptr<RDFDetail::RNodeBase> *prevNode,
                           const std::type_info &helperArgType, const std::type_info &at, void *helperArgOnHeap,
                           TTree *tree, const unsigned int nSlots, const RBookedDefines &customCols, RDataSource *ds,
                           std::weak_ptr<RJittedAction> *jittedActionOnHeap)
{
   // retrieve type of result of the action as a string
   auto helperArgClass = TClass::GetClass(helperArgType);
   if (!helperArgClass) {
      std::string exceptionText = "An error occurred while inferring the result type of an operation.";
      throw std::runtime_error(exceptionText.c_str());
   }
   const auto helperArgClassName = helperArgClass->GetName();

   // retrieve type of action as a string
   auto actionTypeClass = TClass::GetClass(at);
   if (!actionTypeClass) {
      std::string exceptionText = "An error occurred while inferring the action type of the operation.";
      throw std::runtime_error(exceptionText.c_str());
   }
   const std::string actionTypeName = actionTypeClass->GetName();
   const std::string actionTypeNameBase = actionTypeName.substr(actionTypeName.rfind(':') + 1);

   auto definesCopy = new RBookedDefines(customCols); // deleted in jitted CallBuildAction
   auto definesAddr = PrettyPrintAddr(definesCopy);

   // Build a call to CallBuildAction with the appropriate argument. When run through the interpreter, this code will
   // just-in-time create an RAction object and it will assign it to its corresponding RJittedAction.
   std::stringstream createAction_str;
   createAction_str << "ROOT::Internal::RDF::CallBuildAction<" << actionTypeName;
   const auto columnTypeNames =
      GetValidatedArgTypes(cols, customCols, tree, ds, actionTypeNameBase, /*vector2rvec=*/true);
   for (auto &colType : columnTypeNames)
      createAction_str << ", " << colType;
   // on Windows, to prefix the hexadecimal value of a pointer with '0x',
   // one need to write: std::hex << std::showbase << (size_t)pointer
   createAction_str << ">(reinterpret_cast<std::shared_ptr<ROOT::Detail::RDF::RNodeBase>*>("
                    << PrettyPrintAddr(prevNode) << "), new const char*[" << cols.size() << "]{";
   for (auto i = 0u; i < cols.size(); ++i) {
      if (i != 0u)
         createAction_str << ", ";
      createAction_str << '"' << cols[i] << '"';
   }
   createAction_str << "}, " << cols.size() << ", " << nSlots << ", reinterpret_cast<" << helperArgClassName << "*>("
                    << PrettyPrintAddr(helperArgOnHeap)
                    << "), reinterpret_cast<std::weak_ptr<ROOT::Internal::RDF::RJittedAction>*>("
                    << PrettyPrintAddr(jittedActionOnHeap)
                    << "), reinterpret_cast<ROOT::Internal::RDF::RBookedDefines*>(" << definesAddr << "));";
   return createAction_str.str();
}

bool AtLeastOneEmptyString(const std::vector<std::string_view> strings)
{
   for (const auto &s : strings) {
      if (s.empty())
         return true;
   }
   return false;
}

std::shared_ptr<RNodeBase> UpcastNode(std::shared_ptr<RNodeBase> ptr)
{
   return ptr;
}

/// Given the desired number of columns and the user-provided list of columns:
/// * fallback to using the first nColumns default columns if needed (or throw if nColumns > nDefaultColumns)
/// * check that selected column names refer to valid branches, custom columns or datasource columns (throw if not)
/// * replace column names from aliases by the actual column name
/// Return the list of selected column names.
ColumnNames_t GetValidatedColumnNames(RLoopManager &lm, const unsigned int nColumns, const ColumnNames_t &columns,
                                      const ColumnNames_t &validDefines, RDataSource *ds)
{
   auto selectedColumns = SelectColumns(nColumns, columns, lm.GetDefaultColumnNames());

   // Resolve aliases and expand `#var` to `__rdf_sizeof_var`
   const auto &aliasMap = lm.GetAliasMap();

   for (auto &col : selectedColumns) {
      col = ResolveAlias(col, aliasMap);
   }

   // Complain if there are still unknown columns at this point
   const auto unknownColumns = FindUnknownColumns(selectedColumns, lm.GetBranchNames(), validDefines,
                                                  ds ? ds->GetColumnNames() : ColumnNames_t{});

   if (!unknownColumns.empty()) {
      std::stringstream unknowns;
      std::string delim = unknownColumns.size() > 1 ? "s: " : ": "; // singular/plural
      for (auto &unknownColumn : unknownColumns) {
         unknowns << delim << unknownColumn;
         delim = ',';
      }
      throw std::runtime_error("Unknown column" + unknowns.str());
   }

   return selectedColumns;
}

std::vector<std::string> GetValidatedArgTypes(const ColumnNames_t &colNames, const RBookedDefines &defines, TTree *tree,
                                              RDataSource *ds, const std::string &context, bool vector2rvec)
{
   auto toCheckedArgType = [&](const std::string &c) {
      RDFDetail::RDefineBase *define = defines.HasName(c) ? defines.GetColumns().at(c).get() : nullptr;
      const auto colType = ColumnName2ColumnTypeName(c, tree, ds, define, vector2rvec);
      if (colType.rfind("CLING_UNKNOWN_TYPE", 0) == 0) { // the interpreter does not know this type
         const auto msg =
            "The type of custom column \"" + c + "\" (" + colType.substr(19) +
            ") is not known to the interpreter, but a just-in-time-compiled " + context +
            " call requires this column. Make sure to create and load ROOT dictionaries for this column's class.";
         throw std::runtime_error(msg);
      }
      return colType;
   };
   std::vector<std::string> colTypes;
   colTypes.reserve(colNames.size());
   std::transform(colNames.begin(), colNames.end(), std::back_inserter(colTypes), toCheckedArgType);
   return colTypes;
}

/// Return a bitset each element of which indicates whether the corresponding element in `selectedColumns` is the
/// name of a column that must be defined via datasource. All elements of the returned vector are false if no
/// data-source is present.
std::vector<bool> FindUndefinedDSColumns(const ColumnNames_t &requestedCols, const ColumnNames_t &definedCols)
{
   const auto nColumns = requestedCols.size();
   std::vector<bool> mustBeDefined(nColumns, false);
   for (auto i = 0u; i < nColumns; ++i)
      mustBeDefined[i] = std::find(definedCols.begin(), definedCols.end(), requestedCols[i]) == definedCols.end();
   return mustBeDefined;
}

void CheckForDuplicateSnapshotColumns(const ColumnNames_t &cols)
{
   std::unordered_set<std::string> uniqueCols;
   for (auto &col : cols) {
      if (!uniqueCols.insert(col).second) {
         const auto msg = "Error: column \"" + col +
                          "\" was passed to Snapshot twice. This is not supported: only one of the columns would be "
                          "readable with RDataFrame.";
         throw std::logic_error(msg);
      }
   }
}

} // namespace RDF
} // namespace Internal
} // namespace ROOT
