//===-- SystemZTargetMachine.cpp - Define TargetMachine for SystemZ -------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "SystemZTargetMachine.h"
#include "llvm/CodeGen/Passes.h"
#include "llvm/Support/TargetRegistry.h"
#include "llvm/Transforms/Scalar.h"

using namespace llvm;

extern "C" void LLVMInitializeSystemZTarget() {
  // Register the target.
  RegisterTargetMachine<SystemZTargetMachine> X(TheSystemZTarget);
}

SystemZTargetMachine::SystemZTargetMachine(const Target &T, StringRef TT,
                                           StringRef CPU, StringRef FS,
                                           const TargetOptions &Options,
                                           Reloc::Model RM, CodeModel::Model CM,
                                           CodeGenOpt::Level OL)
    : LLVMTargetMachine(T, TT, CPU, FS, Options, RM, CM, OL),
      Subtarget(TT, CPU, FS, *this) {
  initAsmInfo();
}

namespace {
/// SystemZ Code Generator Pass Configuration Options.
class SystemZPassConfig : public TargetPassConfig {
public:
  SystemZPassConfig(SystemZTargetMachine *TM, PassManagerBase &PM)
    : TargetPassConfig(TM, PM) {}

  SystemZTargetMachine &getSystemZTargetMachine() const {
    return getTM<SystemZTargetMachine>();
  }

  void addIRPasses() override;
  bool addInstSelector() override;
  bool addPreSched2() override;
  bool addPreEmitPass() override;
};
} // end anonymous namespace

void SystemZPassConfig::addIRPasses() {
  TargetPassConfig::addIRPasses();
}

bool SystemZPassConfig::addInstSelector() {
  addPass(createSystemZISelDag(getSystemZTargetMachine(), getOptLevel()));
  return false;
}

bool SystemZPassConfig::addPreSched2() {
  if (getOptLevel() != CodeGenOpt::None &&
      getSystemZTargetMachine().getSubtargetImpl()->hasLoadStoreOnCond())
    addPass(&IfConverterID);
  return true;
}

bool SystemZPassConfig::addPreEmitPass() {
  // We eliminate comparisons here rather than earlier because some
  // transformations can change the set of available CC values and we
  // generally want those transformations to have priority.  This is
  // especially true in the commonest case where the result of the comparison
  // is used by a single in-range branch instruction, since we will then
  // be able to fuse the compare and the branch instead.
  //
  // For example, two-address NILF can sometimes be converted into
  // three-address RISBLG.  NILF produces a CC value that indicates whether
  // the low word is zero, but RISBLG does not modify CC at all.  On the
  // other hand, 64-bit ANDs like NILL can sometimes be converted to RISBG.
  // The CC value produced by NILL isn't useful for our purposes, but the
  // value produced by RISBG can be used for any comparison with zero
  // (not just equality).  So there are some transformations that lose
  // CC values (while still being worthwhile) and others that happen to make
  // the CC result more useful than it was originally.
  //
  // Another reason is that we only want to use BRANCH ON COUNT in cases
  // where we know that the count register is not going to be spilled.
  //
  // Doing it so late makes it more likely that a register will be reused
  // between the comparison and the branch, but it isn't clear whether
  // preventing that would be a win or not.
  if (getOptLevel() != CodeGenOpt::None)
    addPass(createSystemZElimComparePass(getSystemZTargetMachine()));
  if (getOptLevel() != CodeGenOpt::None)
    addPass(createSystemZShortenInstPass(getSystemZTargetMachine()));
  addPass(createSystemZLongBranchPass(getSystemZTargetMachine()));
  return true;
}

TargetPassConfig *SystemZTargetMachine::createPassConfig(PassManagerBase &PM) {
  return new SystemZPassConfig(this, PM);
}
