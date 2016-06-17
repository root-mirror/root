//--------------------------------------------------------------------*- C++ -*-
// CLING - the C++ LLVM-based InterpreterG :)
// author:  Axel Naumann <axel@cern.ch>
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_INCREMENTAL_JIT_H
#define CLING_INCREMENTAL_JIT_H

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "llvm/IR/Mangler.h"
#include "llvm/IR/GlobalValue.h"
#include "llvm/ExecutionEngine/JITEventListener.h"
#include "llvm/ExecutionEngine/Orc/CompileUtils.h"
#include "llvm/ExecutionEngine/Orc/IRCompileLayer.h"
#include "llvm/ExecutionEngine/Orc/LazyEmittingLayer.h"
#include "llvm/ExecutionEngine/Orc/ObjectLinkingLayer.h"
#include "llvm/ExecutionEngine/RTDyldMemoryManager.h"
#include "llvm/Target/TargetMachine.h"

namespace llvm {
  class Module;
  class RTDyldMemoryManager;
}

namespace cling {
class Azog;
class IncrementalExecutor;

class IncrementalJIT {
  friend class Azog;

  ///\brief The IncrementalExecutor who owns us.
  IncrementalExecutor& m_Parent;
  llvm::JITEventListener* m_GDBListener; // owned by llvm::ManagedStaticBase

  struct JITSymbol: public llvm::orc::JITSymbol {
    JITSymbol(): llvm::orc::JITSymbol(nullptr) {}
    JITSymbol& operator=(const llvm::orc::JITSymbol& RHS) {
      llvm::orc::JITSymbol::operator=(RHS);
      return *this;
    }
  };
  llvm::StringMap<JITSymbol> m_SymbolMap;

  class NotifyObjectLoadedT {
  public:
    typedef std::vector<std::unique_ptr<llvm::object::OwningBinary<llvm::object::ObjectFile>>> ObjListT;
    typedef std::vector<std::unique_ptr<llvm::RuntimeDyld::LoadedObjectInfo>>
        LoadedObjInfoListT;

    NotifyObjectLoadedT(IncrementalJIT &jit) : m_JIT(jit) {}

    void operator()(llvm::orc::ObjectLinkingLayerBase::ObjSetHandleT H,
                    const ObjListT &Objects,
                    const LoadedObjInfoListT &Infos) const {
      m_JIT.m_UnfinalizedSections[H]
        = std::move(m_JIT.m_SectionsAllocatedSinceLastLoad);
      m_JIT.m_SectionsAllocatedSinceLastLoad = SectionAddrSet();
      assert(Objects.size() == Infos.size() &&
             "Incorrect number of Infos for Objects.");
      if (auto GDBListener = m_JIT.m_GDBListener) {
        for (size_t I = 0, N = Objects.size(); I < N; ++I)
          GDBListener->NotifyObjectEmitted(*Objects[I]->getBinary(),
                                           *Infos[I]);
      }

      for (const auto& Object: Objects) {
        for (const auto &Symbol: Object->getBinary()->symbols()) {
          auto Flags = Symbol.getFlags();
          if (Flags & llvm::object::BasicSymbolRef::SF_Undefined)
            continue;
          // FIXME: this should be uncommented once we serve incremental
          // modules from a TU module.
          //if (!(Flags & llvm::object::BasicSymbolRef::SF_Exported))
          //  continue;
          auto NameOrError = Symbol.getName();
          assert(NameOrError);
          auto Name = NameOrError.get();
          if (m_JIT.m_SymbolMap.find(Name) == m_JIT.m_SymbolMap.end()) {
            llvm::orc::JITSymbol Sym = m_JIT.m_CompileLayer.findSymbolIn(H, Name, true);
            if (Sym.getAddress())
              m_JIT.m_SymbolMap[Name] = Sym;
          }
        }
      }
    };

  private:
    IncrementalJIT &m_JIT;
  };

  typedef llvm::orc::ObjectLinkingLayer<NotifyObjectLoadedT> ObjectLayerT;
  typedef llvm::orc::IRCompileLayer<ObjectLayerT> CompileLayerT;
  typedef llvm::orc::LazyEmittingLayer<CompileLayerT> LazyEmitLayerT;
  typedef LazyEmitLayerT::ModuleSetHandleT ModuleSetHandleT;

  std::unique_ptr<llvm::TargetMachine> m_TM;
  llvm::DataLayout m_TMDataLayout;

  ///\brief The RTDyldMemoryManager used to communicate with the
  /// IncrementalExecutor to handle missing or special symbols.
  std::unique_ptr<llvm::RTDyldMemoryManager> m_ExeMM;

  NotifyObjectLoadedT m_NotifyObjectLoaded;

  ObjectLayerT m_ObjectLayer;
  CompileLayerT m_CompileLayer;
  LazyEmitLayerT m_LazyEmitLayer;

  // We need to store ObjLayerT::ObjSetHandles for each of the object sets
  // that have been emitted but not yet finalized so that we can forward the
  // mapSectionAddress calls appropriately.
  typedef std::set<const void *> SectionAddrSet;
  struct ObjSetHandleCompare {
    bool operator()(ObjectLayerT::ObjSetHandleT H1,
                    ObjectLayerT::ObjSetHandleT H2) const {
      return &*H1 < &*H2;
    }
  };
  SectionAddrSet m_SectionsAllocatedSinceLastLoad;
  std::map<ObjectLayerT::ObjSetHandleT, SectionAddrSet, ObjSetHandleCompare>
      m_UnfinalizedSections;

  ///\brief Vector of ModuleSetHandleT. UnloadHandles index into that
  /// vector.
  std::vector<ModuleSetHandleT> m_UnloadPoints;


  std::string Mangle(llvm::StringRef Name) {
    std::string MangledName;
    {
      llvm::raw_string_ostream MangledNameStream(MangledName);
      llvm::Mangler::getNameWithPrefix(MangledNameStream, Name,
                                       m_TMDataLayout);
    }
    return MangledName;
  }

  llvm::orc::JITSymbol getInjectedSymbols(llvm::StringRef Name) const;

public:
  IncrementalJIT(IncrementalExecutor& exe,
                 std::unique_ptr<llvm::TargetMachine> TM);

  ///\brief Get the address of a symbol from the JIT or the memory manager,
  /// mangling the name as needed. Use this to resolve symbols as coming
  /// from clang's mangler.
  /// \param Name - name to look for. This name might still get mangled
  ///   (prefixed by '_') to make IR versus symbol names.
  /// \param AlsoInProcess - Sometimes you only care about JITed symbols. If so,
  ///   pass `false` here to not resolve the symbol through dlsym().
  uint64_t getSymbolAddress(llvm::StringRef Name, bool AlsoInProcess) {
    return getSymbolAddressWithoutMangling(Mangle(Name), AlsoInProcess)
      .getAddress();
  }

  ///\brief Get the address of a symbol from the JIT or the memory manager.
  /// Use this to resolve symbols of known, target-specific names.
  llvm::orc::JITSymbol getSymbolAddressWithoutMangling(llvm::StringRef Name,
                                                       bool AlsoInProcess);

  size_t addModules(std::vector<llvm::Module*>&& modules);
  void removeModules(size_t handle);

  IncrementalExecutor& getParent() const { return m_Parent; }

  void
  RemoveUnfinalizedSection(llvm::orc::ObjectLinkingLayerBase::ObjSetHandleT H) {
    m_UnfinalizedSections.erase(H);
  }
};
} // end cling
#endif // CLING_INCREMENTAL_EXECUTOR_H
