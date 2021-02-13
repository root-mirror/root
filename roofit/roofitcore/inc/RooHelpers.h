// Author: Stephan Hageboeck, CERN  01/2019

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_ROOFITCORE_INC_ROOHELPERS_H_
#define ROOFIT_ROOFITCORE_INC_ROOHELPERS_H_

#include "RooMsgService.h"
#include "RooAbsArg.h"
#include "RooAbsReal.h"

#include <sstream>
#include <vector>
#include <string>

namespace RooHelpers {

/// Switches the message service to a different level while the instance is alive.
/// Can also temporarily activate / deactivate message topics.
/// Use as
/// ~~~{.cpp}
/// RooHelpers::LocalChangeMessageLevel changeMsgLvl(RooFit::WARNING);
/// [ statements that normally generate a lot of output ]
/// ~~~
class LocalChangeMsgLevel {
  public:
    /// Change message level (and topics) while this object is alive, reset when it goes out of scope.
    /// \param[in] lvl The desired message level. Defaults to verbose.
    /// \param[in] extraTopics Extra topics to be switched on. These will only switched on in the last stream to prevent all streams are printing.
    /// \param[in] removeTopics Message topics to be switched off
    /// \param[in] overrideExternalLevel Override the user message level.
    LocalChangeMsgLevel(RooFit::MsgLevel lvl = RooFit::DEBUG,
        unsigned int extraTopics = 0u,
        unsigned int removeTopics = 0u,
        bool overrideExternalLevel = true);

    ~LocalChangeMsgLevel();

  private:
    RooFit::MsgLevel fOldKillBelow;
    std::vector<RooMsgService::StreamConfig> fOldConf;
    int fExtraStream{-1};
};


/// Wrap an object into a TObject. Sometimes needed to avoid reinterpret_cast or enable RTTI.
template<typename T>
struct WrapIntoTObject : public TObject {
  WrapIntoTObject(T& obj) : _payload(&obj) { }
  T* _payload;
};


/// Hijacks all messages with given level and topic (and optionally object name) while alive.
/// Use this like an ostringstream afterwards. The messages can e.g. be retrieved using `str()`.
/// Useful for unit tests / debugging.
class HijackMessageStream{
  public:
    HijackMessageStream(RooFit::MsgLevel level, RooFit::MsgTopic topics, const char* objectName = nullptr);
    template<typename T>
    const HijackMessageStream& operator<<(const T& v) const {
      _str << v;
      return *this;
    }
    std::string str() { return _str.str(); }
    std::ostringstream& stream() { return _str; };
    ~HijackMessageStream();

  private:
    std::ostringstream _str;
    RooFit::MsgLevel _oldKillBelow;
    std::vector<RooMsgService::StreamConfig> _oldConf;
    Int_t _thisStream;
};


std::vector<std::string> tokenise(const std::string &str, const std::string &delims, bool returnEmptyToken = true);


/// Check if the parameters have a range, and warn if the range extends below / above the set limits.
void checkRangeOfParameters(const RooAbsReal* callingClass, std::initializer_list<const RooAbsReal*> pars,
    double min = -std::numeric_limits<double>::max(), double max = std::numeric_limits<double>::max(),
    bool limitsInAllowedRange = false, std::string extraMessage = "");


/// Disable all caches for sub-branches in an expression tree.
/// This is helpful when an expression with cached sub-branches needs to be integrated numerically.
struct DisableCachingRAII {
  /// Inhibit all dirty-state propagation, and assume every node as dirty.
  /// \param[in] oldState Restore this state when going out of scope.
  DisableCachingRAII(bool oldState):
  _oldState(oldState) {
    RooAbsArg::setDirtyInhibit(true);
  }

  ~DisableCachingRAII() {
    RooAbsArg::setDirtyInhibit(_oldState);
  }
  bool _oldState;
};



}



#endif /* ROOFIT_ROOFITCORE_INC_ROOHELPERS_H_ */
