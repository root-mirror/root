//------------------------------------------------------------------------------
// CLING - the C++ LLVM-based InterpreterG :)
// author:  Lukasz Janyst <ljanyst@cern.ch>
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#include "cling/Interpreter/Interpreter.h"
#include "cling/Utils/Paths.h"
#include "ClingUtils.h"

#include "DynamicLookup.h"
#include "ExternalInterpreterSource.h"
#include "ForwardDeclPrinter.h"
#include "IncrementalExecutor.h"
#include "IncrementalParser.h"
#include "MultiplexInterpreterCallbacks.h"
#include "TransactionUnloader.h"

#include "cling/Interpreter/CIFactory.h"
#include "cling/Interpreter/ClangInternalState.h"
#include "cling/Interpreter/ClingCodeCompleteConsumer.h"
#include "cling/Interpreter/CompilationOptions.h"
#include "cling/Interpreter/DynamicLibraryManager.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/Interpreter/Value.h"
#include "cling/Interpreter/AutoloadCallback.h"
#include "cling/Utils/AST.h"
#include "cling/Utils/SourceNormalization.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/Basic/TargetInfo.h"
#include "clang/Basic/SourceManager.h"
#include "clang/CodeGen/ModuleBuilder.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/Utils.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Lex/HeaderSearch.h"
#include "clang/Parse/Parser.h"
#include "clang/Sema/Sema.h"
#include "clang/Sema/SemaDiagnostic.h"

#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/Path.h"
#include "llvm/Support/raw_ostream.h"

#include <sstream>
#include <string>
#include <vector>

using namespace clang;

namespace {

  static cling::Interpreter::ExecutionResult
  ConvertExecutionResult(cling::IncrementalExecutor::ExecutionResult ExeRes) {
    switch (ExeRes) {
    case cling::IncrementalExecutor::kExeSuccess:
      return cling::Interpreter::kExeSuccess;
    case cling::IncrementalExecutor::kExeFunctionNotCompiled:
      return cling::Interpreter::kExeFunctionNotCompiled;
    case cling::IncrementalExecutor::kExeUnresolvedSymbols:
      return cling::Interpreter::kExeUnresolvedSymbols;
    default: break;
    }
    return cling::Interpreter::kExeSuccess;
  }

  static bool isPracticallyEmptyModule(const llvm::Module* M) {
    return M->empty() && M->global_empty() && M->alias_empty();
  }
} // unnamed namespace

namespace cling {

  Interpreter::PushTransactionRAII::PushTransactionRAII(const Interpreter* i)
    : m_Interpreter(i) {
    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = 0;
    CO.Debug = 0;
    CO.CodeGeneration = 1;
    CO.CodeGenerationForModule = 0;

    m_Transaction = m_Interpreter->m_IncrParser->beginTransaction(CO);
  }

  Interpreter::PushTransactionRAII::~PushTransactionRAII() {
    pop();
  }

  void Interpreter::PushTransactionRAII::pop() const {
    IncrementalParser::ParseResultTransaction PRT
      = m_Interpreter->m_IncrParser->endTransaction(m_Transaction);
    if (PRT.getPointer()) {
      assert(PRT.getPointer()==m_Transaction && "Ended different transaction?");
      m_Interpreter->m_IncrParser->commitTransaction(PRT);
    }
  }

  Interpreter::StateDebuggerRAII::StateDebuggerRAII(const Interpreter* i)
    : m_Interpreter(i) {
    if (m_Interpreter->isPrintingDebug()) {
      const CompilerInstance& CI = *m_Interpreter->getCI();
      CodeGenerator* CG = i->m_IncrParser->getCodeGenerator();

      // The ClangInternalState constructor can provoke deserialization,
      // we need a transaction.
      PushTransactionRAII pushedT(i);

      m_State.reset(new ClangInternalState(CI.getASTContext(),
                                           CI.getPreprocessor(),
                                           CG ? CG->GetModule() : 0,
                                           CG,
                                           "aName"));
    }
  }

  Interpreter::StateDebuggerRAII::~StateDebuggerRAII() {
    if (m_State) {
      // The ClangInternalState destructor can provoke deserialization,
      // we need a transaction.
      PushTransactionRAII pushedT(m_Interpreter);
      m_State->compare("aName", m_Interpreter->m_Opts.Verbose());
      m_State.reset();
    }
  }

  const Parser& Interpreter::getParser() const {
    return *m_IncrParser->getParser();
  }

  Parser& Interpreter::getParser() {
    return *m_IncrParser->getParser();
  }

  clang::SourceLocation Interpreter::getNextAvailableLoc() const {
    return m_IncrParser->getLastMemoryBufferEndLoc().getLocWithOffset(1);
  }


  bool Interpreter::isInSyntaxOnlyMode() const {
    return getCI()->getFrontendOpts().ProgramAction
      == clang::frontend::ParseSyntaxOnly;
  }

  Interpreter::Interpreter(int argc, const char* const *argv,
                           const char* llvmdir /*= 0*/, bool noRuntime,
                           const Interpreter* parentInterp) :
    m_Opts(argc, argv),
    m_UniqueCounter(parentInterp ? parentInterp->m_UniqueCounter + 1 : 0),
    m_PrintDebug(false), m_DynamicLookupDeclared(false),
    m_DynamicLookupEnabled(false), m_RawInputEnabled(false) {

    m_LLVMContext.reset(new llvm::LLVMContext);
    m_DyLibManager.reset(new DynamicLibraryManager(getOptions()));
    m_IncrParser.reset(new IncrementalParser(this, llvmdir));

    Sema& SemaRef = getSema();
    Preprocessor& PP = SemaRef.getPreprocessor();
    // Enable incremental processing, which prevents the preprocessor destroying
    // the lexer on EOF token.
    PP.enableIncrementalProcessing();

    m_LookupHelper.reset(new LookupHelper(new Parser(PP, SemaRef,
                                                     /*SkipFunctionBodies*/false,
                                                     /*isTemp*/true), this));

    if (!isInSyntaxOnlyMode())
      m_Executor.reset(new IncrementalExecutor(SemaRef.Diags, *getCI()));

    // Tell the diagnostic client that we are entering file parsing mode.
    DiagnosticConsumer& DClient = getCI()->getDiagnosticClient();
    DClient.BeginSourceFile(getCI()->getLangOpts(), &PP);

    llvm::SmallVector<IncrementalParser::ParseResultTransaction, 2>
      IncrParserTransactions;
    m_IncrParser->Initialize(IncrParserTransactions, parentInterp);

    handleFrontendOptions();

    if (!noRuntime) {
      if (getCI()->getLangOpts().CPlusPlus)
        IncludeCXXRuntime();
      else
        IncludeCRuntime();
    }
    // Commit the transactions, now that gCling is set up. It is needed for
    // static initialization in these transactions through local_cxa_atexit().
    for (auto&& I: IncrParserTransactions)
      m_IncrParser->commitTransaction(I);
    // Disable suggestions for ROOT
    bool showSuggestions = !llvm::StringRef(ClingStringify(CLING_VERSION)).startswith("ROOT");

    // We need InterpreterCallbacks only if it is a parent Interpreter.
    if (!parentInterp) {
      std::unique_ptr<InterpreterCallbacks>
         AutoLoadCB(new AutoloadCallback(this, showSuggestions));
      setCallbacks(std::move(AutoLoadCB));
    }

    m_IncrParser->SetTransformers(parentInterp);
  }

  ///\brief Constructor for the child Interpreter.
  /// Passing the parent Interpreter as an argument.
  ///
  Interpreter::Interpreter(const Interpreter &parentInterpreter, int argc,
                           const char* const *argv,
                           const char* llvmdir /*= 0*/, bool noRuntime) :
    Interpreter(argc, argv, llvmdir, noRuntime, &parentInterpreter) {
    // Do the "setup" of the connection between this interpreter and
    // its parent interpreter.

    // The "bridge" between the interpreters.
    ExternalInterpreterSource *myExternalSource =
      new ExternalInterpreterSource(&parentInterpreter, this);

    llvm::IntrusiveRefCntPtr <ExternalASTSource>
      astContextExternalSource(myExternalSource);

    getCI()->getASTContext().setExternalSource(astContextExternalSource);

    // Inform the Translation Unit Decl of I2 that it has to search somewhere
    // else to find the declarations.
    getCI()->getASTContext().getTranslationUnitDecl()->setHasExternalVisibleStorage(true);

    // Give my IncrementalExecutor a pointer to the Incremental executor of the
    // parent Interpreter.
    m_Executor->setExternalIncrementalExecutor(parentInterpreter.m_Executor.get());
  }

  Interpreter::~Interpreter() {
    if (m_Executor)
      m_Executor->shuttingDown();
    for (size_t i = 0, e = m_StoredStates.size(); i != e; ++i)
      delete m_StoredStates[i];
    getCI()->getDiagnostics().getClient()->EndSourceFile();
    // LookupHelper's ~Parser needs the PP from IncrParser's CI, so do this
    // first:
    m_LookupHelper.reset();

    // We want to keep the callback alive during the shutdown of Sema, CodeGen
    // and the ASTContext. For that to happen we shut down the IncrementalParser
    // explicitly, before the implicit destruction (through the unique_ptr) of
    // the callbacks.
    m_IncrParser.reset(0);
  }

  const char* Interpreter::getVersion() const {
    return ClingStringify(CLING_VERSION);
  }

  void Interpreter::handleFrontendOptions() {
    if (m_Opts.ShowVersion) {
      llvm::errs() << getVersion() << '\n';
    }
    if (m_Opts.Help) {
      m_Opts.PrintHelp();
    }
  }

  void Interpreter::IncludeCXXRuntime() {
    // Set up common declarations which are going to be available
    // only at runtime
    // Make sure that the universe won't be included to compile time by using
    // -D __CLING__ as CompilerInstance's arguments

    std::stringstream initializer;
#ifdef _WIN32
    // We have to use the #defined __CLING__ on windows first.
    //FIXME: Find proper fix.
    initializer << "#ifdef __CLING__ \n#endif\n";
#endif

    initializer << "#include \"cling/Interpreter/RuntimeUniverse.h\"\n";

    if (!isInSyntaxOnlyMode()) {
      // Set up the gCling variable if it can be used
      initializer << "namespace cling {namespace runtime { "
        "cling::Interpreter *gCling=(cling::Interpreter*)"
        << "0x" << std::hex << (uintptr_t)this << " ;} }";
    }
    declare(initializer.str());
  }

  void Interpreter::IncludeCRuntime() {
    // Set up the gCling variable if it can be used
    std::stringstream initializer;
    initializer << "void* gCling=(void*)" << (uintptr_t)this << ';';
    declare(initializer.str());
    // declare("void setValueNoAlloc(void* vpI, void* vpSVR, void* vpQT);");
    // declare("void setValueNoAlloc(void* vpI, void* vpV, void* vpQT, float value);");
    // declare("void setValueNoAlloc(void* vpI, void* vpV, void* vpQT, double value);");
    // declare("void setValueNoAlloc(void* vpI, void* vpV, void* vpQT, long double value);");
    // declare("void setValueNoAlloc(void* vpI, void* vpV, void* vpQT, unsigned long long value);");
    // declare("void setValueNoAlloc(void* vpI, void* vpV, void* vpQT, const void* value);");
    // declare("void* setValueWithAlloc(void* vpI, void* vpV, void* vpQT);");

    declare("#include \"cling/Interpreter/CValuePrinter.h\"");
  }

  void Interpreter::AddIncludePaths(llvm::StringRef PathStr, const char* Delm) {
    CompilerInstance* CI = getCI();
    HeaderSearchOptions& HOpts = CI->getHeaderSearchOpts();

    // Save the current number of entries
    size_t Idx = HOpts.UserEntries.size();
    utils::AddIncludePaths(PathStr, HOpts, Delm);

    Preprocessor& PP = CI->getPreprocessor();
    SourceManager& SM = PP.getSourceManager();
    FileManager& FM = SM.getFileManager();
    HeaderSearch& HSearch = PP.getHeaderSearchInfo();
    const bool isFramework = false;

    // Add all the new entries into Preprocessor
    for (const size_t N = HOpts.UserEntries.size(); Idx < N; ++Idx) {
      const HeaderSearchOptions::Entry& E = HOpts.UserEntries[Idx];
      if (const clang::DirectoryEntry *DE = FM.getDirectory(E.Path)) {
        HSearch.AddSearchPath(DirectoryLookup(DE, SrcMgr::C_User, isFramework),
                              E.Group == frontend::Angled);
      }
    }
  }

  void Interpreter::DumpIncludePath(llvm::raw_ostream* S) {
    utils::DumpIncludePaths(getCI()->getHeaderSearchOpts(), S ? *S : llvm::outs(),
                            true /*withSystem*/, true /*withFlags*/);
  }

  void Interpreter::storeInterpreterState(const std::string& name) const {
    // This may induce deserialization
    PushTransactionRAII RAII(this);
    CodeGenerator* CG = m_IncrParser->getCodeGenerator();
    ClangInternalState* state
      = new ClangInternalState(getCI()->getASTContext(),
                               getCI()->getPreprocessor(),
                               getLastTransaction()->getModule(),
                               CG, name);
    m_StoredStates.push_back(state);
  }

  void Interpreter::compareInterpreterState(const std::string& name) const {
    short foundAtPos = -1;
    for (short i = 0, e = m_StoredStates.size(); i != e; ++i) {
      if (m_StoredStates[i]->getName() == name) {
        foundAtPos = i;
        break;
      }
    }
    if (foundAtPos < 0) {
      llvm::errs() << "The store point name " << name << " does not exist."
      "Unbalanced store / compare\n";
      return;
    }

    // This may induce deserialization
    PushTransactionRAII RAII(this);
    m_StoredStates[foundAtPos]->compare(name, m_Opts.Verbose());
  }

  void Interpreter::printIncludedFiles(llvm::raw_ostream& Out) const {
    ClangInternalState::printIncludedFiles(Out, getCI()->getSourceManager());
  }


  void Interpreter::GetIncludePaths(llvm::SmallVectorImpl<std::string>& incpaths,
                                   bool withSystem, bool withFlags) {
    utils::CopyIncludePaths(getCI()->getHeaderSearchOpts(), incpaths,
                            withSystem, withFlags);
  }

  CompilerInstance* Interpreter::getCI() const {
    return m_IncrParser->getCI();
  }

  Sema& Interpreter::getSema() const {
    return getCI()->getSema();
  }

  ///\brief Maybe transform the input line to implement cint command line
  /// semantics (declarations are global) and compile to produce a module.
  ///
  Interpreter::CompilationResult
  Interpreter::process(const std::string& input, Value* V /* = 0 */,
                       Transaction** T /* = 0 */) {
    std::string wrapReadySource = input;
    size_t wrapPoint = std::string::npos;
    if (!isRawInputEnabled())
      wrapPoint = utils::getWrapPoint(wrapReadySource, getCI()->getLangOpts());

    if (isRawInputEnabled() || wrapPoint == std::string::npos) {
      CompilationOptions CO;
      CO.DeclarationExtraction = 0;
      CO.ValuePrinting = 0;
      CO.ResultEvaluation = 0;
      CO.DynamicScoping = isDynamicLookupEnabled();
      CO.Debug = isPrintingDebug();
      CO.CheckPointerValidity = 1;
      return DeclareInternal(input, CO, T);
    }

    CompilationOptions CO;
    CO.DeclarationExtraction = 1;
    CO.ValuePrinting = CompilationOptions::VPAuto;
    CO.ResultEvaluation = (bool)V;
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();
    CO.CheckPointerValidity = 1;
    if (EvaluateInternal(wrapReadySource, CO, V, T, wrapPoint)
                                                     == Interpreter::kFailure) {
      return Interpreter::kFailure;
    }

    return Interpreter::kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::parse(const std::string& input, Transaction** T /*=0*/) const {
    CompilationOptions CO;
    CO.CodeGeneration = 0;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();

    return DeclareInternal(input, CO, T);
  }

  Interpreter::CompilationResult
  Interpreter::loadModuleForHeader(const std::string& headerFile) {
    Preprocessor& PP = getCI()->getPreprocessor();
    //Copied from clang's PPDirectives.cpp
    bool isAngled = false;
    // Clang doc says:
    // "LookupFrom is set when this is a \#include_next directive, it specifies
    // the file to start searching from."
    const DirectoryLookup* FromDir = 0;
    const FileEntry* FromFile = 0;
    const DirectoryLookup* CurDir = 0;

    ModuleMap::KnownHeader suggestedModule;
    // PP::LookupFile uses it to issue 'nice' diagnostic
    SourceLocation fileNameLoc;
    PP.LookupFile(fileNameLoc, headerFile, isAngled, FromDir, FromFile, CurDir,
                  /*SearchPath*/0, /*RelativePath*/ 0, &suggestedModule,
                  /*SkipCache*/false, /*OpenFile*/ false, /*CacheFail*/ false);
    if (!suggestedModule)
      return Interpreter::kFailure;

    // Copied from PPDirectives.cpp
    SmallVector<std::pair<IdentifierInfo *, SourceLocation>, 2> path;
    for (Module *mod = suggestedModule.getModule(); mod; mod = mod->Parent) {
      IdentifierInfo* II
        = &getSema().getPreprocessor().getIdentifierTable().get(mod->Name);
      path.push_back(std::make_pair(II, fileNameLoc));
    }

    std::reverse(path.begin(), path.end());

    // Pretend that the module came from an inclusion directive, so that clang
    // will create an implicit import declaration to capture it in the AST.
    bool isInclude = true;
    SourceLocation includeLoc;
    if (getCI()->loadModule(includeLoc, path, Module::AllVisible, isInclude)) {
      // After module load we need to "force" Sema to generate the code for
      // things like dynamic classes.
      getSema().ActOnEndOfTranslationUnit();
      return Interpreter::kSuccess;
    }

    return Interpreter::kFailure;
  }

  Interpreter::CompilationResult
  Interpreter::parseForModule(const std::string& input) {
    CompilationOptions CO;
    CO.CodeGeneration = 1;
    CO.CodeGenerationForModule = 1;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();

    // When doing parseForModule avoid warning about the user code
    // being loaded ... we probably might as well extend this to
    // ALL warnings ... but this will suffice for now (working
    // around a real bug in QT :().
    DiagnosticsEngine& Diag = getCI()->getDiagnostics();
    Diag.setSeverity(clang::diag::warn_field_is_uninit,
                     clang::diag::Severity::Ignored, SourceLocation());
    CompilationResult Result = DeclareInternal(input, CO);
    Diag.setSeverity(clang::diag::warn_field_is_uninit,
                     clang::diag::Severity::Warning, SourceLocation());
    return Result;
  }


  Interpreter::CompilationResult
  Interpreter::CodeCompleteInternal(const std::string& input, unsigned offset) {

    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();
    CO.CheckPointerValidity = 0;

    std::string wrapped = input;
    size_t wrapPos = utils::getWrapPoint(wrapped, getCI()->getLangOpts());
    const std::string& Src = WrapInput(wrapped, wrapped, wrapPos);

    CO.CodeCompletionOffset = offset + wrapPos;

    StateDebuggerRAII stateDebugger(this);

    // This triggers the FileEntry to be created and the completion
    // point to be set in clang.
    m_IncrParser->Compile(Src, CO);

    return kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::declare(const std::string& input, Transaction** T/*=0 */) {
    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();
    CO.CheckPointerValidity = 0;

    return DeclareInternal(input, CO, T);
  }

  Interpreter::CompilationResult
  Interpreter::evaluate(const std::string& input, Value& V) {
    // Here we might want to enforce further restrictions like: Only one
    // ExprStmt can be evaluated and etc. Such enforcement cannot happen in the
    // worker, because it is used from various places, where there is no such
    // rule
    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 1;

    return EvaluateInternal(input, CO, &V);
  }

  Interpreter::CompilationResult
  Interpreter::codeComplete(const std::string& line, size_t& cursor,
                            std::vector<std::string>& completions) const {

    const char * const argV = "cling";
    std::string resourceDir = this->getCI()->getHeaderSearchOpts().ResourceDir;
    // Remove the extra 3 directory names "/lib/clang/3.9.0"
    StringRef parentResourceDir = llvm::sys::path::parent_path(
                                  llvm::sys::path::parent_path(
                                  llvm::sys::path::parent_path(resourceDir)));
    std::string llvmDir = parentResourceDir.str();

    Interpreter childInterpreter(*this, 1, &argV, llvmDir.c_str());

    auto childCI = childInterpreter.getCI();
    clang::Sema &childSemaRef = childCI->getSema();

    // Create the CodeCompleteConsumer with InterpreterCallbacks
    // from the parent interpreter and set the consumer for the child
    // interpreter.
    ClingCodeCompleteConsumer* consumer = new ClingCodeCompleteConsumer(
                getCI()->getFrontendOpts().CodeCompleteOpts, completions);
    // Child interpreter CI will own consumer!
    childCI->setCodeCompletionConsumer(consumer);
    childSemaRef.CodeCompleter = consumer;

    // Ignore diagnostics when we tab complete.
    // This is because we get redefinition errors due to the import of the decls.
    clang::IgnoringDiagConsumer* ignoringDiagConsumer =
                                            new clang::IgnoringDiagConsumer();                      
    childSemaRef.getDiagnostics().setClient(ignoringDiagConsumer, true);
    DiagnosticsEngine& parentDiagnostics = this->getCI()->getSema().getDiagnostics();

    std::unique_ptr<DiagnosticConsumer> ownerDiagConsumer = 
                                                parentDiagnostics.takeClient();
    auto clientDiagConsumer = parentDiagnostics.getClient();
    parentDiagnostics.setClient(ignoringDiagConsumer, /*owns*/ false);

    // The child will desirialize decls from *this. We need a transaction RAII.
    PushTransactionRAII RAII(this);

    // Triger the code completion.
    childInterpreter.CodeCompleteInternal(line, cursor);

    // Restore the original diagnostics client for parent interpreter.
    parentDiagnostics.setClient(clientDiagConsumer,
                                ownerDiagConsumer.release() != nullptr);
    parentDiagnostics.Reset(/*soft=*/true);

    return kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::echo(const std::string& input, Value* V /* = 0 */) {
    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = CompilationOptions::VPEnabled;
    CO.ResultEvaluation = (bool)V;

    return EvaluateInternal(input, CO, V);
  }

  Interpreter::CompilationResult
  Interpreter::execute(const std::string& input) {
    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = 0;
    CO.Debug = isPrintingDebug();
    return EvaluateInternal(input, CO);
  }

  Interpreter::CompilationResult Interpreter::emitAllDecls(Transaction* T) {
    assert(!isInSyntaxOnlyMode() && "No CodeGenerator?");
    m_IncrParser->emitTransaction(T);
    m_IncrParser->addTransaction(T);
    T->setState(Transaction::kCollecting);
    auto PRT = m_IncrParser->endTransaction(T);
    m_IncrParser->commitTransaction(PRT);

    if ((T = PRT.getPointer()))
      if (executeTransaction(*T))
        return Interpreter::kSuccess;

    return Interpreter::kFailure;
  }

  const std::string& Interpreter::WrapInput(const std::string& Input,
                                            std::string& Output,
                                            size_t& WrapPoint) const {
    // If wrapPoint is > length of input, nothing is wrapped!
    if (WrapPoint < Input.size()) {
      std::string Header("void ");
      Header.append(createUniqueWrapper());
      Header.append("(void* vpClingValue) {\n ");

      // Suppport Input and Output begin the same string
      std::string Wrapper = Input.substr(WrapPoint);
      Wrapper.insert(0, Header);
      Wrapper.append("\n;\n}");
      Wrapper.insert(0, Input.substr(0, WrapPoint));
      Wrapper.swap(Output);
      WrapPoint += Header.size();
      return Output;
    }
    // in-case std::string::npos was passed
    WrapPoint = 0;
    return Input;
  }

  Interpreter::ExecutionResult
  Interpreter::RunFunction(const FunctionDecl* FD, Value* res /*=0*/) {
    if (getCI()->getDiagnostics().hasErrorOccurred())
      return kExeCompilationError;

    if (isInSyntaxOnlyMode()) {
      return kExeNoCodeGen;
    }

    if (!FD)
      return kExeUnkownFunction;

    std::string mangledNameIfNeeded;
    utils::Analyze::maybeMangleDeclName(FD, mangledNameIfNeeded);
    IncrementalExecutor::ExecutionResult ExeRes =
       m_Executor->executeWrapper(mangledNameIfNeeded, res);
    return ConvertExecutionResult(ExeRes);
  }

  const FunctionDecl* Interpreter::DeclareCFunction(StringRef name,
                                                    StringRef code,
                                                    bool withAccessControl) {
    /*
    In CallFunc we currently always (intentionally and somewhat necessarily)
    always fully specify member function template, however this can lead to
    an ambiguity with a class template.  For example in
    roottest/cling/functionTemplate we get:

    input_line_171:3:15: warning: lookup of 'set' in member access expression
    is ambiguous; using member of 't'
    ((t*)obj)->set<int>(*(int*)args[0]);
               ^
    roottest/cling/functionTemplate/t.h:19:9: note: lookup in the object type
    't' refers here
    void set(T targ) {
         ^
    /usr/include/c++/4.4.5/bits/stl_set.h:87:11: note: lookup from the
    current scope refers here
    class set
          ^
    This is an intention warning implemented in clang, see
    http://llvm.org/viewvc/llvm-project?view=revision&revision=105518

    which 'should have been' an error:

    C++ [basic.lookup.classref] requires this to be an error, but,
    because it's hard to work around, Clang downgrades it to a warning as
    an extension.</p>

    // C++98 [basic.lookup.classref]p1:
    // In a class member access expression (5.2.5), if the . or -> token is
    // immediately followed by an identifier followed by a <, the identifier
    // must be looked up to determine whether the < is the beginning of a
    // template argument list (14.2) or a less-than operator. The identifier
    // is first looked up in the class of the object expression. If the
    // identifier is not found, it is then looked up in the context of the
    // entire postfix-expression and shall name a class or function template. If
    // the lookup in the class of the object expression finds a template, the
    // name is also looked up in the context of the entire postfix-expression
    // and
    // -- if the name is not found, the name found in the class of the
    // object expression is used, otherwise
    // -- if the name is found in the context of the entire postfix-expression
    // and does not name a class template, the name found in the class of the
    // object expression is used, otherwise
    // -- if the name found is a class template, it must refer to the same
    // entity as the one found in the class of the object expression,
    // otherwise the program is ill-formed.

    See -Wambiguous-member-template

    An alternative to disabling the diagnostics is to use a pointer to
    member function:

    #include <set>
    using namespace std;

    extern "C" int printf(const char*,...);

    struct S {
    template <typename T>
    void set(T) {};

    virtual void virtua() { printf("S\n"); }
    };

    struct T: public S {
    void virtua() { printf("T\n"); }
    };

    int main() {
    S *s = new T();
    typedef void (S::*Func_p)(int);
    Func_p p = &S::set<int>;
    (s->*p)(12);

    typedef void (S::*Vunc_p)(void);
    Vunc_p q = &S::virtua;
    (s->*q)(); // prints "T"
    return 0;
    }
    */
    DiagnosticsEngine& Diag = getCI()->getDiagnostics();
    Diag.setSeverity(clang::diag::ext_nested_name_member_ref_lookup_ambiguous,
                     clang::diag::Severity::Ignored, SourceLocation());


    LangOptions& LO = const_cast<LangOptions&>(getCI()->getLangOpts());
    bool savedAccessControl = LO.AccessControl;
    LO.AccessControl = withAccessControl;
    cling::Transaction* T = 0;
    cling::Interpreter::CompilationResult CR = declare(code, &T);
    LO.AccessControl = savedAccessControl;

    Diag.setSeverity(clang::diag::ext_nested_name_member_ref_lookup_ambiguous,
                     clang::diag::Severity::Warning, SourceLocation());

    if (CR != cling::Interpreter::kSuccess)
      return 0;

    for (cling::Transaction::const_iterator I = T->decls_begin(),
           E = T->decls_end(); I != E; ++I) {
      if (I->m_Call != cling::Transaction::kCCIHandleTopLevelDecl)
        continue;
      if (const LinkageSpecDecl* LSD
          = dyn_cast<LinkageSpecDecl>(*I->m_DGR.begin())) {
        DeclContext::decl_iterator DeclBegin = LSD->decls_begin();
        if (DeclBegin == LSD->decls_end())
          continue;
        if (const FunctionDecl* D = dyn_cast<FunctionDecl>(*DeclBegin)) {
          const IdentifierInfo* II = D->getDeclName().getAsIdentifierInfo();
          if (II && II->getName() == name)
            return D;
        }
      }
    }
    return 0;
  }

  void*
  Interpreter::compileFunction(llvm::StringRef name, llvm::StringRef code,
                               bool ifUnique, bool withAccessControl) {
    //
    //  Compile the wrapper code.
    //

    if (isInSyntaxOnlyMode())
      return 0;

    if (ifUnique) {
      if (void* Addr = (void*)getAddressOfGlobal(name)) {
        return Addr;
      }
    }

    const FunctionDecl* FD = DeclareCFunction(name, code, withAccessControl);
    if (!FD)
      return 0;
    //
    //  Get the wrapper function pointer
    //  from the ExecutionEngine (the JIT).
    //
    if (const llvm::GlobalValue* GV
        = getLastTransaction()->getModule()->getNamedValue(name))
      return m_Executor->getPointerToGlobalFromJIT(*GV);

    return 0;
  }

  void*
  Interpreter::compileDtorCallFor(const clang::RecordDecl* RD) {
    void* &addr = m_DtorWrappers[RD];
    if (addr)
      return addr;

    std::string funcname;
    {
      llvm::raw_string_ostream namestr(funcname);
      namestr << "__cling_Destruct_" << RD;
    }

    std::string code = "extern \"C\" void ";
    clang::QualType RDQT(RD->getTypeForDecl(), 0);
    std::string typeName
      = utils::TypeName::GetFullyQualifiedName(RDQT, RD->getASTContext());
    std::string dtorName = RD->getNameAsString();
    code += funcname + "(void* obj){((" + typeName + "*)obj)->~"
      + dtorName + "();}";

    // ifUniq = false: we know it's unique, no need to check.
    addr = compileFunction(funcname, code, false /*ifUniq*/,
                           false /*withAccessControl*/);
    return addr;
  }

  void Interpreter::createUniqueName(std::string& out) {
    out += utils::Synthesize::UniquePrefix;
    llvm::raw_string_ostream(out) << m_UniqueCounter++;
  }

  bool Interpreter::isUniqueName(llvm::StringRef name) {
    return name.startswith(utils::Synthesize::UniquePrefix);
  }

  llvm::StringRef Interpreter::createUniqueWrapper() const {
    llvm::SmallString<128> out(utils::Synthesize::UniquePrefix);
    llvm::raw_svector_ostream(out) << m_UniqueCounter++;
    return (getCI()->getASTContext().Idents.getOwn(out)).getName();
  }

  bool Interpreter::isUniqueWrapper(llvm::StringRef name) {
    return name.startswith(utils::Synthesize::UniquePrefix);
  }

  Interpreter::CompilationResult
  Interpreter::DeclareInternal(const std::string& input,
                               const CompilationOptions& CO,
                               Transaction** T /* = 0 */) const {
    assert(CO.DeclarationExtraction == 0
           && CO.ValuePrinting == 0
           && CO.ResultEvaluation == 0
           && "Compilation Options not compatible with \"declare\" mode.");

    StateDebuggerRAII stateDebugger(this);

    IncrementalParser::ParseResultTransaction PRT
      = m_IncrParser->Compile(input, CO);
    if (PRT.getInt() == IncrementalParser::kFailed)
      return Interpreter::kFailure;

    if (T)
      *T = PRT.getPointer();
    return Interpreter::kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::EvaluateInternal(const std::string& input,
                                CompilationOptions CO,
                                Value* V, /* = 0 */
                                Transaction** T /* = 0 */,
                                size_t wrapPoint /* = 0*/) {
    StateDebuggerRAII stateDebugger(this);

    // Wrap the expression
    std::string WrapperBuffer;
    const std::string& Wrapper = WrapInput(input, WrapperBuffer, wrapPoint);

    // We have wrapped and need to disable warnings that are caused by
    // non-default C++ at the prompt:
    CO.IgnorePromptDiags = 1;

    IncrementalParser::ParseResultTransaction PRT
      = m_IncrParser->Compile(Wrapper, CO);
    Transaction* lastT = PRT.getPointer();
    if (lastT && lastT->getState() != Transaction::kCommitted) {
      assert((lastT->getState() == Transaction::kCommitted
              || lastT->getState() == Transaction::kRolledBack
              || lastT->getState() == Transaction::kRolledBackWithErrors)
             && "Not committed?");
      if (V)
        *V = Value();
      return kFailure;
    }

    // Might not have a Transaction
    if (PRT.getInt() == IncrementalParser::kFailed) {
      if (V)
        *V = Value();
      return kFailure;
    }

    if (!lastT) {
      // Empty transactions are good, too!
      if (V)
        *V = Value();
      return kSuccess;
    }

    Value resultV;
    if (!V)
      V = &resultV;
    if (!lastT->getWrapperFD()) // no wrapper to run
      return Interpreter::kSuccess;
    else if (RunFunction(lastT->getWrapperFD(), V) < kExeFirstError){
      if (lastT->getCompilationOpts().ValuePrinting
          != CompilationOptions::VPDisabled
          && V->isValid()
          // the !V->needsManagedAllocation() case is handled by
          // dumpIfNoStorage.
          && V->needsManagedAllocation())
        V->dump();
      return Interpreter::kSuccess;
    }
    return Interpreter::kSuccess;
  }

  std::string Interpreter::lookupFileOrLibrary(llvm::StringRef file) {
    std::string canonicalFile = DynamicLibraryManager::normalizePath(file);
    if (canonicalFile.empty())
      canonicalFile = file;
    const FileEntry* FE = 0;

    //Copied from clang's PPDirectives.cpp
    bool isAngled = false;
    // Clang doc says:
    // "LookupFrom is set when this is a \#include_next directive, it
    // specifies the file to start searching from."
    const DirectoryLookup* FromDir = 0;
    const FileEntry* FromFile = 0;
    const DirectoryLookup* CurDir = 0;
    Preprocessor& PP = getCI()->getPreprocessor();
    // PP::LookupFile uses it to issue 'nice' diagnostic
    SourceLocation fileNameLoc;
    FE = PP.LookupFile(fileNameLoc, canonicalFile, isAngled, FromDir, FromFile,
                       CurDir, /*SearchPath*/0, /*RelativePath*/ 0,
                       /*suggestedModule*/0, /*SkipCache*/false,
                       /*OpenFile*/ false, /*CacheFail*/ false);
    if (FE)
      return FE->getName();
    return getDynamicLibraryManager()->lookupLibrary(canonicalFile);
  }

  Interpreter::CompilationResult
  Interpreter::loadFile(const std::string& filename,
                        bool allowSharedLib /*=true*/,
                        Transaction** T /*= 0*/) {
    DynamicLibraryManager* DLM = getDynamicLibraryManager();
    std::string canonicalLib = DLM->lookupLibrary(filename);
    if (allowSharedLib && !canonicalLib.empty()) {
      switch (DLM->loadLibrary(canonicalLib, /*permanent*/false, /*resolved*/true)) {
      case DynamicLibraryManager::kLoadLibSuccess: // Intentional fall through
      case DynamicLibraryManager::kLoadLibAlreadyLoaded:
        return kSuccess;
      case DynamicLibraryManager::kLoadLibNotFound:
        assert(0 && "Cannot find library with existing canonical name!");
        return kFailure;
      default:
        // Not a source file (canonical name is non-empty) but can't load.
        return kFailure;
      }
    }

    std::string code;
    code += "#include \"" + filename + "\"";

    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();
    CO.CheckPointerValidity = 1;
    CompilationResult res = DeclareInternal(code, CO, T);
    return res;
  }

  void Interpreter::unload(Transaction& T) {
    if (InterpreterCallbacks* callbacks = getCallbacks())
      callbacks->TransactionUnloaded(T);
    if (m_Executor) { // we also might be in fsyntax-only mode.
      m_Executor->runAndRemoveStaticDestructors(&T);
      if (!T.getExecutor()) {
        // this transaction might be queued in the executor
        m_Executor->unloadFromJIT(T.getModule(),
                                  Transaction::ExeUnloadHandle({(void*)(size_t)-1}));
      }
    }

    // We can revert the most recent transaction or a nested transaction of a
    // transaction that is not in the middle of the transaction collection
    // (i.e. at the end or not yet added to the collection at all).
    assert(!T.getTopmostParent()->getNext() &&
           "Can not revert previous transactions");
    assert((T.getState() != Transaction::kRolledBack ||
            T.getState() != Transaction::kRolledBackWithErrors) &&
           "Transaction already rolled back.");
    if (getOptions().ErrorOut)
      return;

    if (InterpreterCallbacks* callbacks = getCallbacks())
      callbacks->TransactionRollback(T);

    TransactionUnloader U(this, &getCI()->getSema(),
                          m_IncrParser->getCodeGenerator(),
                          m_Executor.get());
    if (U.RevertTransaction(&T))
      T.setState(Transaction::kRolledBack);
    else
      T.setState(Transaction::kRolledBackWithErrors);

    m_IncrParser->deregisterTransaction(T);
  }

  void Interpreter::unload(unsigned numberOfTransactions) {
    while(true) {
      cling::Transaction* T = m_IncrParser->getLastTransaction();
      if (!T) {
        llvm::errs() << "cling: invalid last transaction; unload failed!\n";
        return;
      }
      unload(*T);
      if (!--numberOfTransactions)
        break;
    }

  }

  static void runAndRemoveStaticDestructorsImpl(IncrementalExecutor &executor,
                                std::vector<const Transaction*> &transactions,
                                         unsigned int begin, unsigned int end) {

    for(auto i = begin; i != end; --i) {
      if (transactions[i-1] != nullptr)
        executor.runAndRemoveStaticDestructors(const_cast<Transaction*>(transactions[i-1]));
    }
  }

  void Interpreter::runAndRemoveStaticDestructors(unsigned numberOfTransactions) {
    if (!m_Executor)
      return;
    auto transactions( m_IncrParser->getAllTransactions() );
    unsigned int min = 0;
    if (transactions.size() > numberOfTransactions) {
      min = transactions.size() - numberOfTransactions;
    }
    runAndRemoveStaticDestructorsImpl(*m_Executor, transactions,
                                      transactions.size(), min);
  }

  void Interpreter::runAndRemoveStaticDestructors() {
    if (!m_Executor)
      return;
    auto transactions( m_IncrParser->getAllTransactions() );
    runAndRemoveStaticDestructorsImpl(*m_Executor, transactions,
                                      transactions.size(), 0);
  }

  void Interpreter::installLazyFunctionCreator(void* (*fp)(const std::string&)) {
    if (m_Executor)
      m_Executor->installLazyFunctionCreator(fp);
  }

  Value Interpreter::Evaluate(const char* expr, DeclContext* DC,
                                       bool ValuePrinterReq) {
    Sema& TheSema = getCI()->getSema();
    // The evaluation should happen on the global scope, because of the wrapper
    // that is created.
    //
    // We can't PushDeclContext, because we don't have scope.
    Sema::ContextRAII pushDC(TheSema,
                             TheSema.getASTContext().getTranslationUnitDecl());

    Value Result;
    getCallbacks()->SetIsRuntime(true);
    if (ValuePrinterReq)
      echo(expr, &Result);
    else
      evaluate(expr, Result);
    getCallbacks()->SetIsRuntime(false);

    return Result;
  }

  void Interpreter::setCallbacks(std::unique_ptr<InterpreterCallbacks> C) {
    // We need it to enable LookupObject callback.
    if (!m_Callbacks) {
      m_Callbacks.reset(new MultiplexInterpreterCallbacks(this));
      // FIXME: Move to the InterpreterCallbacks.cpp;
      if (DynamicLibraryManager* DLM = getDynamicLibraryManager())
        DLM->setCallbacks(m_Callbacks.get());
    }

    static_cast<MultiplexInterpreterCallbacks*>(m_Callbacks.get())
      ->addCallback(std::move(C));
  }

  const Transaction* Interpreter::getFirstTransaction() const {
    return m_IncrParser->getFirstTransaction();
  }

  const Transaction* Interpreter::getLastTransaction() const {
    return m_IncrParser->getLastTransaction();
  }

  const Transaction* Interpreter::getCurrentTransaction() const {
    return m_IncrParser->getCurrentTransaction();
  }


  void Interpreter::enableDynamicLookup(bool value /*=true*/) {
    if (!m_DynamicLookupDeclared && value) {
      // No dynlookup for the dynlookup header!
      m_DynamicLookupEnabled = false;
      if (loadModuleForHeader("cling/Interpreter/DynamicLookupRuntimeUniverse.h")
          != kSuccess)
      declare("#include \"cling/Interpreter/DynamicLookupRuntimeUniverse.h\"");
    }
    m_DynamicLookupDeclared = true;

    // Enable it *after* parsing the headers.
    m_DynamicLookupEnabled = value;
  }

  Interpreter::ExecutionResult
  Interpreter::executeTransaction(Transaction& T) {
    assert(!isInSyntaxOnlyMode() && "Running on what?");
    assert(T.getState() == Transaction::kCommitted && "Must be committed");

    IncrementalExecutor::ExecutionResult ExeRes
       = IncrementalExecutor::kExeSuccess;
    if (!isPracticallyEmptyModule(T.getModule())) {
      T.setExeUnloadHandle(m_Executor.get(), m_Executor->emitToJIT());

      // Forward to IncrementalExecutor; should not be called by
      // anyone except for IncrementalParser.
      ExeRes = m_Executor->runStaticInitializersOnce(T);
    }

    return ConvertExecutionResult(ExeRes);
  }

  bool Interpreter::addSymbol(const char* symbolName,  void* symbolAddress) {
    // Forward to IncrementalExecutor;
    if (!symbolName || !symbolAddress )
      return false;

    return m_Executor->addSymbol(symbolName, symbolAddress);
  }

  void Interpreter::addModule(llvm::Module* module) {
     m_Executor->addModule(module);
  }


  void* Interpreter::getAddressOfGlobal(const GlobalDecl& GD,
                                        bool* fromJIT /*=0*/) const {
    // Return a symbol's address, and whether it was jitted.
    std::string mangledName;
    utils::Analyze::maybeMangleDeclName(GD, mangledName);
    return getAddressOfGlobal(mangledName, fromJIT);
  }

  void* Interpreter::getAddressOfGlobal(llvm::StringRef SymName,
                                        bool* fromJIT /*=0*/) const {
    // Return a symbol's address, and whether it was jitted.
    if (isInSyntaxOnlyMode())
      return 0;
    return m_Executor->getAddressOfGlobal(SymName, fromJIT);
  }

  void Interpreter::AddAtExitFunc(void (*Func) (void*), void* Arg) {
    m_Executor->AddAtExitFunc(Func, Arg);
  }

  void Interpreter::GenerateAutoloadingMap(llvm::StringRef inFile,
                                           llvm::StringRef outFile,
                                           bool enableMacros,
                                           bool enableLogs) {

    const char *const dummy="cling_fwd_declarator";
    // Create an interpreter without any runtime, producing the fwd decls.
    // FIXME: CIFactory appends extra 3 folders to the llvmdir.
    std::string llvmdir
      = getCI()->getHeaderSearchOpts().ResourceDir + "/../../../";
    cling::Interpreter fwdGen(1, &dummy, llvmdir.c_str(), true);

    // Copy the same header search options to the new instance.
    Preprocessor& fwdGenPP = fwdGen.getCI()->getPreprocessor();
    HeaderSearchOptions headerOpts = getCI()->getHeaderSearchOpts();
    clang::ApplyHeaderSearchOptions(fwdGenPP.getHeaderSearchInfo(), headerOpts,
                                    fwdGenPP.getLangOpts(),
                                    fwdGenPP.getTargetInfo().getTriple());


    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.ValuePrinting = 0;
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = 0;
    CO.Debug = isPrintingDebug();


    std::string includeFile = std::string("#include \"") + inFile.str() + "\"";
    IncrementalParser::ParseResultTransaction PRT
      = fwdGen.m_IncrParser->Compile(includeFile, CO);
    cling::Transaction* T = PRT.getPointer();

    // If this was already #included we will get a T == 0.
    if (PRT.getInt() == IncrementalParser::kFailed || !T)
      return;

    std::error_code EC;
    llvm::raw_fd_ostream out(outFile.data(), EC,
                             llvm::sys::fs::OpenFlags::F_None);
    llvm::raw_fd_ostream log((outFile + ".skipped").str().c_str(),
                             EC, llvm::sys::fs::OpenFlags::F_None);
    log << "Generated for :" << inFile << "\n";
    forwardDeclare(*T, fwdGenPP, fwdGen.getCI()->getSema().getASTContext(),
                   out, enableMacros,
                   &log);
  }

  void Interpreter::forwardDeclare(Transaction& T, Preprocessor& P,
                                   clang::ASTContext& Ctx,
                                   llvm::raw_ostream& out,
                                   bool enableMacros /*=false*/,
                                   llvm::raw_ostream* logs /*=0*/,
                                   IgnoreFilesFunc_t ignoreFiles /*= return always false*/) const {
    llvm::raw_null_ostream null;
    if (!logs)
      logs = &null;

    ForwardDeclPrinter visitor(out, *logs, P, Ctx, T, 0, false, ignoreFiles);
    visitor.printStats();

    // Avoid assertion in the ~IncrementalParser.
    T.setState(Transaction::kCommitted);
  }



} //end namespace cling
