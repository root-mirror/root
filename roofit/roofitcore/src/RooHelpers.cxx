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

#include "RooHelpers.h"
#include "RooAbsRealLValue.h"

#include "TClass.h"

namespace RooHelpers {

LocalChangeMsgLevel::LocalChangeMsgLevel(RooFit::MsgLevel lvl,
    unsigned int extraTopics, unsigned int removeTopics, bool overrideExternalLevel) {
  auto& msg = RooMsgService::instance();
  fOldKillBelow = msg.globalKillBelow();
  if (overrideExternalLevel) msg.setGlobalKillBelow(lvl);

  for (int i = 0; i < msg.numStreams(); ++i) {
    fOldConf.push_back(msg.getStream(i));
    if (overrideExternalLevel) msg.getStream(i).minLevel = lvl;
    msg.getStream(i).removeTopic(static_cast<RooFit::MsgTopic>(removeTopics));
    msg.setStreamStatus(i, true);
  }

  if (extraTopics != 0) {
    fExtraStream = msg.addStream(lvl);
    msg.getStream(fExtraStream).addTopic(static_cast<RooFit::MsgTopic>(extraTopics));
  }
}

LocalChangeMsgLevel::~LocalChangeMsgLevel() {
  auto& msg = RooMsgService::instance();
  msg.setGlobalKillBelow(fOldKillBelow);
  for (int i=0; i < msg.numStreams(); ++i) {
    if (i < static_cast<int>(fOldConf.size()))
      msg.getStream(i) = fOldConf[i];
  }

  if (fExtraStream > 0)
    msg.deleteStream(fExtraStream);
}


/// Tokenise the string by splitting at the characters in delims.
/// Consecutive delimiters are collapsed, so that no delimiters will appear in the
/// tokenised strings, and no emtpy strings are returned.
/// \param[in] str String to tokenise.
/// \param[in] delims One or more delimiters used to split the string.
/// \param[in] returnEmptyToken If the string is empty, return one empty token. Default is to return an empty vector.
std::vector<std::string> tokenise(const std::string &str, const std::string &delims, bool returnEmptyToken /*= true*/) {
  if (str.empty())
    return std::vector<std::string>(returnEmptyToken ? 1 : 0);

  std::vector<std::string> tokens;

  auto beg = str.find_first_not_of(delims, 0);
  auto end = str.find_first_of(delims, beg);
  do {
    tokens.emplace_back(str.substr(beg, end-beg));
    beg = str.find_first_not_of(delims, end);
    end = str.find_first_of(delims, beg);
  } while (beg != std::string::npos);

  return tokens;
}


/// Hijack all messages with given level and topics while this object is alive.
/// \param[in] level Minimum level to hijack. Higher levels also get captured.
/// \param[in] topics Topics to hijack. Use `|` to combine different topics, and cast to `RooFit::MsgTopic` if necessary.
/// \param[in] objectName Only hijack messages from an object with the given name. Defaults to any object.
HijackMessageStream::HijackMessageStream(RooFit::MsgLevel level, RooFit::MsgTopic topics, const char* objectName)
{
  auto& msg = RooMsgService::instance();
  _oldKillBelow = msg.globalKillBelow();
  msg.setGlobalKillBelow(level);
  for (int i = 0; i < msg.numStreams(); ++i) {
    _oldConf.push_back(msg.getStream(i));
    msg.getStream(i).removeTopic(topics);
    msg.setStreamStatus(i, true);
  }

  _thisStream = msg.addStream(level,
      RooFit::Topic(topics),
      RooFit::OutputStream(_str),
      objectName ? RooFit::ObjectName(objectName) : RooCmdArg());
}

HijackMessageStream::~HijackMessageStream() {
  auto& msg = RooMsgService::instance();
  msg.setGlobalKillBelow(_oldKillBelow);
  for (unsigned int i = 0; i < _oldConf.size(); ++i) {
    msg.getStream(i) = _oldConf[i];
  }
  msg.deleteStream(_thisStream);
}


/// \param[in] callingClass Class that's calling. Needed to include name and type name of the class in error message.
/// \param[in] pars List of all parameters to be checked.
/// \param[in] min Minimum of allowed range. `min` itself counts as disallowed.
/// \param[in] max Maximum of allowed range. `max` itself counts as disallowed.
/// \param[in] limitsInAllowedRange If true, the limits passed as parameters are part of the allowed range.
/// \param[in] extraMessage Message that should be appended to the warning.
void checkRangeOfParameters(const RooAbsReal* callingClass, std::initializer_list<const RooAbsReal*> pars,
    double min, double max, bool limitsInAllowedRange, std::string extraMessage) {
  const char openBr = limitsInAllowedRange ? '[' : '(';
  const char closeBr = limitsInAllowedRange ? ']' : ')';

  for (auto parameter : pars) {
    auto par = dynamic_cast<const RooAbsRealLValue*>(parameter);
    if (par && (
        (par->getMin() < min || par->getMax() > max)
        || (!limitsInAllowedRange && (par->getMin() == min || par->getMax() == max)) )) {
      std::stringstream rangeMsg;
      rangeMsg << openBr;
      if (min > -std::numeric_limits<double>::max())
        rangeMsg << min << ", ";
      else
        rangeMsg << "-inf, ";

      if (max < std::numeric_limits<double>::max())
        rangeMsg << max << closeBr;
      else
        rangeMsg << "inf" << closeBr;

      oocoutW(callingClass, InputArguments) << "The parameter '" << par->GetName() << "' with range [" << par->getMin("") << ", "
          << par->getMax() << "] of the " << callingClass->IsA()->GetName() << " '" << callingClass->GetName()
          << "' exceeds the safe range of " << rangeMsg.str() << ". Advise to limit its range."
          << (!extraMessage.empty() ? "\n" : "") << extraMessage << std::endl;
    }
  }
}

}
