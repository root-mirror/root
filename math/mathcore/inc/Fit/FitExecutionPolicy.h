#include "ROOT/EExecutionPolicy.hxx"
#include "ROOT/RConfig.hxx"

namespace ROOT {
namespace Fit {
namespace ExecutionPolicy {

constexpr ROOT::EExecutionPolicy
   _R__DEPRECATED_626("ROOT::Fit::ExecutionPolicy is being deprecated."
                      "Use ROOT::EExecutionPolicy::kMultiThread instead.")
      kMultithread = ROOT::EExecutionPolicy::kMultiThread;

constexpr ROOT::EExecutionPolicy
   _R__DEPRECATED_626("ROOT::Fit::ExecutionPolicy is being deprecated."
                      "Use ROOT::EExecutionPolicy::kMultiProcess instead.")
      kMultiprocess = ROOT::EExecutionPolicy::kMultiProcess;

constexpr ROOT::EExecutionPolicy _R__DEPRECATED_626(
   "ROOT::Fit::ExecutionPolicy is being deprecated."
   "Use ROOT::EExecutionPolicy::kSequential instead.")
      kSerial = ROOT::EExecutionPolicy::kSequential;

} // namespace ExecutionPolicy
} // Fit ns
} // namespace ROOT
