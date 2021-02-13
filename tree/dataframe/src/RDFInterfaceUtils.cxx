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
#include <set>
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

static const std::string &ResolveAlias(const std::string &col, const std::map<std::string, std::string> &aliasMap)
{
   const auto it = aliasMap.find(col);
   if (it != aliasMap.end())
      return it->second;
   return col;
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
         const auto &col = ResolveAlias(columnOrAlias, aliasMap);
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

static ParsedExpression ParseRDFExpression(const std::string &expr, const ColumnNames_t &treeBranchNames,
                                           const ColumnNames_t &customColNames, const ColumnNames_t &dataSourceColNames,
                                           const std::map<std::string, std::string> &aliasMap)
{
   const auto usedColsAndAliases = FindUsedColumns(expr, treeBranchNames, customColNames, dataSourceColNames, aliasMap);

   auto escapeDots = [](const std::string &s) {
      TString ss(s);
      TPRegexp dot("\\.");
      dot.Substitute(ss, "\\.", "g");
      return std::string(std::move(ss));
   };

   ColumnNames_t varNames;
   ColumnNames_t usedCols;
   TString exprWithVars(expr); // same as expr but column names will be substituted with the variable names in varNames
   for (const auto &colOrAlias : usedColsAndAliases) {
      const auto col = ResolveAlias(colOrAlias, aliasMap);
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
   R__ASSERT(vars.size() == varTypes.size());

   TPRegexp re(R"(\breturn\b)");
   const bool hasReturnStmt = re.Match(expr) == 1;

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
                                       std::set<TTree *> &analysedTrees)
{

   if (!analysedTrees.insert(&t).second) {
      return;
   }

   auto branches = t.GetListOfBranches();
   if (branches) {
      for (auto branchObj : *branches) {
         auto name = branchObj->GetName();
         if (bNamesReg.insert(name).second) {
            bNames.emplace_back(name);
         }
      }
   }

   auto friendTrees = t.GetListOfFriends();

   if (!friendTrees)
      return;

   for (auto friendTreeObj : *friendTrees) {
      auto friendTree = ((TFriendElement *)friendTreeObj)->GetTree();
      GetTopLevelBranchNamesImpl(*friendTree, bNamesReg, bNames, analysedTrees);
   }
}

static bool IsValidCppVarName(const std::string &var)
{
   if (var.empty())
      return false;
   const char firstChar = var[0];

   // first character must be either a letter or an underscore
   auto isALetter = [](char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); };
   const bool isValidFirstChar = firstChar == '_' || isALetter(firstChar);
   if (!isValidFirstChar)
      return false;

   // all characters must be either a letter, an underscore or a number
   auto isANumber = [](char c) { return c >= '0' && c <= '9'; };
   auto isValidTok = [&isALetter, &isANumber](char c) { return c == '_' || isALetter(c) || isANumber(c); };
   for (const char c : var)
      if (!isValidTok(c))
         return false;

   return true;
}

} // anonymous namespace

namespace ROOT {
namespace Internal {
namespace RDF {

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
      if ((isEmptyRegex || 0 != regexp.Match(colName.c_str())) && !RDFInternal::IsInternalColumn(colName)) {
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

void CheckDefine(std::string_view definedCol, TTree *treePtr, const ColumnNames_t &customCols,
                 const std::map<std::string, std::string> &aliasMap, const ColumnNames_t &dataSourceColumns)
{
   const std::string definedColStr(definedCol);

   if (!IsValidCppVarName(definedColStr)) {
      const auto msg = "Cannot define column \"" + definedColStr + "\": not a valid C++ variable name.";
      throw std::runtime_error(msg);
   }

   if (treePtr != nullptr) {
      // check if definedCol is already present in TTree
      const auto branch = treePtr->GetBranch(definedColStr.c_str());
      if (branch != nullptr) {
         const auto msg = "branch \"" + definedColStr + "\" already present in TTree";
         throw std::runtime_error(msg);
      }
   }
   // check if definedCol has already been `Define`d in the functional graph
   if (std::find(customCols.begin(), customCols.end(), definedCol) != customCols.end()) {
      const auto msg = "Redefinition of column \"" + definedColStr + "\"";
      throw std::runtime_error(msg);
   }

   // Check if the definedCol is an alias
   const auto aliasColNameIt = aliasMap.find(definedColStr);
   if (aliasColNameIt != aliasMap.end()) {
      const auto msg = "An alias with name " + definedColStr + " pointing to column " +
      aliasColNameIt->second + " is already existing.";
      throw std::runtime_error(msg);
   }

   // check if definedCol is already present in the DataSource (but has not yet been `Define`d)
   if (!dataSourceColumns.empty()) {
      if (std::find(dataSourceColumns.begin(), dataSourceColumns.end(), definedCol) != dataSourceColumns.end()) {
         const auto msg = "Redefinition of column \"" + definedColStr + "\" already present in the data-source";
         throw std::runtime_error(msg);
      }
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
                   const ColumnNames_t &branches, const RDFInternal::RBookedDefines &customCols, TTree *tree,
                   RDataSource *ds)
{
   const auto &dsColumns = ds ? ds->GetColumnNames() : ColumnNames_t{};

   const auto parsedExpr =
      ParseRDFExpression(std::string(expression), branches, customCols.GetNames(), dsColumns, aliasMap);
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
   filterInvocation << "ROOT::Internal::RDF::JitFilterHelper(" << lambdaName << ", {";
   for (const auto &col : parsedExpr.fUsedCols)
      filterInvocation << "\"" << col << "\", ";
   if (!parsedExpr.fUsedCols.empty())
      filterInvocation.seekp(-2, filterInvocation.cur); // remove the last ",
   // lifetime of pointees:
   // - jittedFilter: heap-allocated weak_ptr to the actual jittedFilter that will be deleted by JitFilterHelper
   // - prevNodeOnHeap: heap-allocated shared_ptr to the actual previous node that will be deleted by JitFilterHelper
   // - definesOnHeap: heap-allocated, will be deleted by JitFilterHelper
   filterInvocation << "}, \"" << name << "\", "
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
                                             RDataSource *ds, const RDFInternal::RBookedDefines &customCols,
                                             const ColumnNames_t &branches,
                                             std::shared_ptr<RNodeBase> *upcastNodeOnHeap)
{
   const auto &aliasMap = lm.GetAliasMap();
   auto *const tree = lm.GetTree();
   const auto &dsColumns = ds ? ds->GetColumnNames() : ColumnNames_t{};

   const auto parsedExpr =
      ParseRDFExpression(std::string(expression), branches, customCols.GetNames(), dsColumns, aliasMap);
   const auto exprVarTypes =
      GetValidatedArgTypes(parsedExpr.fUsedCols, customCols, tree, ds, "Define", /*vector2rvec=*/true);
   const auto lambdaName = DeclareLambda(parsedExpr.fExpr, parsedExpr.fVarNames, exprVarTypes);
   const auto type = RetTypeOfLambda(lambdaName);

   auto definesCopy = new RDFInternal::RBookedDefines(customCols);
   auto definesAddr = PrettyPrintAddr(definesCopy);
   auto jittedDefine = std::make_shared<RDFDetail::RJittedDefine>(name, type, lm.GetNSlots(), lm.GetDSValuePtrs());

   std::stringstream defineInvocation;
   defineInvocation << "ROOT::Internal::RDF::JitDefineHelper(" << lambdaName << ", {";
   for (const auto &col : parsedExpr.fUsedCols) {
      defineInvocation << "\"" << col << "\", ";
   }
   if (!parsedExpr.fUsedCols.empty())
      defineInvocation.seekp(-2, defineInvocation.cur); // remove the last ",
   // lifetime of pointees:
   // - lm is the loop manager, and if that goes out of scope jitting does not happen at all (i.e. will always be valid)
   // - jittedDefine: heap-allocated weak_ptr that will be deleted by JitDefineHelper after usage
   // - definesAddr: heap-allocated, will be deleted by JitDefineHelper after usage
   defineInvocation << "}, \"" << name << "\", reinterpret_cast<ROOT::Detail::RDF::RLoopManager*>("
                    << PrettyPrintAddr(&lm) << "), reinterpret_cast<std::weak_ptr<ROOT::Detail::RDF::RJittedDefine>*>("
                    << PrettyPrintAddr(MakeWeakOnHeap(jittedDefine))
                    << "), reinterpret_cast<ROOT::Internal::RDF::RBookedDefines*>(" << definesAddr
                    << "), reinterpret_cast<std::shared_ptr<ROOT::Detail::RDF::RNodeBase>*>("
                    << PrettyPrintAddr(upcastNodeOnHeap) << "));\n";

   lm.ToJitExec(defineInvocation.str());
   return jittedDefine;
}

// Jit and call something equivalent to "this->BuildAndBook<ColTypes...>(params...)"
// (see comments in the body for actual jitted code)
std::string JitBuildAction(const ColumnNames_t &bl, std::shared_ptr<RDFDetail::RNodeBase> *prevNode,
                           const std::type_info &helperArgType, const std::type_info &at, void *helperArgOnHeap,
                           TTree *tree, const unsigned int nSlots, const RDFInternal::RBookedDefines &customCols,
                           RDataSource *ds, std::weak_ptr<RJittedAction> *jittedActionOnHeap)
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

   auto definesCopy = new RDFInternal::RBookedDefines(customCols); // deleted in jitted CallBuildAction
   auto definesAddr = PrettyPrintAddr(definesCopy);

   // Build a call to CallBuildAction with the appropriate argument. When run through the interpreter, this code will
   // just-in-time create an RAction object and it will assign it to its corresponding RJittedAction.
   std::stringstream createAction_str;
   createAction_str << "ROOT::Internal::RDF::CallBuildAction<" << actionTypeName;
   const auto columnTypeNames =
      GetValidatedArgTypes(bl, customCols, tree, ds, actionTypeNameBase, /*vector2rvec=*/true);
   for (auto &colType : columnTypeNames)
      createAction_str << ", " << colType;
   // on Windows, to prefix the hexadecimal value of a pointer with '0x',
   // one need to write: std::hex << std::showbase << (size_t)pointer
   createAction_str << ">(reinterpret_cast<std::shared_ptr<ROOT::Detail::RDF::RNodeBase>*>("
                    << PrettyPrintAddr(prevNode) << "), {";
   for (auto i = 0u; i < bl.size(); ++i) {
      if (i != 0u)
         createAction_str << ", ";
      createAction_str << '"' << bl[i] << '"';
   }
   createAction_str << "}, " << nSlots << ", reinterpret_cast<" << helperArgClassName << "*>("
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
   const auto &defaultColumns = lm.GetDefaultColumnNames();
   auto selectedColumns = SelectColumns(nColumns, columns, defaultColumns);
   const auto &validBranchNames = lm.GetBranchNames();
   const auto unknownColumns =
      FindUnknownColumns(selectedColumns, validBranchNames, validDefines, ds ? ds->GetColumnNames() : ColumnNames_t{});

   if (!unknownColumns.empty()) {
      // throw
      std::stringstream unknowns;
      std::string delim = unknownColumns.size() > 1 ? "s: " : ": "; // singular/plural
      for (auto &unknownColumn : unknownColumns) {
         unknowns << delim << unknownColumn;
         delim = ',';
      }
      throw std::runtime_error("Unknown column" + unknowns.str());
   }

   // Now we need to check within the aliases if some of the yet unknown names can be recovered
   auto &aliasMap = lm.GetAliasMap();
   auto aliasMapEnd = aliasMap.end();

   for (auto idx : ROOT::TSeqU(selectedColumns.size())) {
      const auto &colName = selectedColumns[idx];
      const auto aliasColumnNameIt = aliasMap.find(colName);
      if (aliasMapEnd != aliasColumnNameIt) {
         selectedColumns[idx] = aliasColumnNameIt->second;
      }
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

} // namespace RDF
} // namespace Internal
} // namespace ROOT
