# -*- coding:utf-8 -*-

#-----------------------------------------------------------------------------
#  Author: Danilo Piparo <Danilo.Piparo@cern.ch> CERN
#  Author: Enric Tejedor <enric.tejedor.saavedra@cern.ch> CERN
#-----------------------------------------------------------------------------

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from JupyROOT.helpers import utils
import ROOT

# Jit a wrapper for the ttabcom
_TTabComHookCode = """
std::vector<std::string> _TTabComHook(const char* pattern){
   static auto ttc = new TTabCom;
   const size_t lineBufSize = 2*1024;  // must be equal to/larger than BUF_SIZE in TTabCom.cxx
   std::unique_ptr<char[]> completed(new char[lineBufSize]);
   strncpy(completed.get(), pattern, lineBufSize);
   completed[lineBufSize-1] = '\\0';
   int pLoc = strlen(completed.get());
   std::ostringstream oss;
   Int_t firstChange = ttc->Hook(completed.get(), &pLoc, oss);
   if (firstChange == -2) { // got some completions in oss
      auto completions = oss.str();
      vector<string> completions_v;
      istringstream f(completions);
      string s;
      while (getline(f, s, '\\n')) {
         completions_v.push_back(s);
      }
      return completions_v;
   }
   if (firstChange == -1) { // found no completions
      return vector<string>();
   }
   // found exactly one completion
   return vector<string>(1, completed.get());
}
"""

class CppCompleter(object):
    '''
    Completer which interfaces to the TTabCom of ROOT. It is activated
    (deactivated) upon the load(unload) of the load of the extension.

    >>> comp = CppCompleter()
    >>> comp.activate()
    >>> for suggestion in comp._completeImpl("TTreeF"):
    ...     print(suggestion)
    TTreeFormula
    TTreeFormulaManager
    TTreeFriendLeafIter
    >>> garbage = ROOT.gInterpreter.ProcessLine("TH1F* h")
    >>> for suggestion in comp._completeImpl("h->GetA"):
    ...     print(suggestion)
    h->GetArray
    h->GetAsymmetry
    h->GetAt
    h->GetAxisColor
    >>> garbage = ROOT.gInterpreter.ProcessLine("TH1F aa")
    >>> for suggestion in comp._completeImpl("aa.Add("):
    ...     print(suggestion.replace("\\t"," "))
    <BLANKLINE>
    Bool_t Add(TF1* h1, Double_t c1 = 1, Option_t* option = "")
    Bool_t Add(const TH1* h, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1)  // *MENU*
    Bool_t Add(const TH1* h1, Double_t c1 = 1)
    >>> for suggestion in comp._completeImpl("TROOT::Is"):
    ...     print(suggestion)
    TROOT::IsA
    TROOT::IsBatch
    TROOT::IsEqual
    TROOT::IsEscaped
    TROOT::IsExecutingMacro
    TROOT::IsFolder
    TROOT::IsInterrupted
    TROOT::IsLineProcessing
    TROOT::IsModified
    TROOT::IsOnHeap
    TROOT::IsProofServ
    TROOT::IsRootFile
    TROOT::IsSortable
    TROOT::IsWebDisplay
    TROOT::IsWebDisplayBatch
    TROOT::IsWritable
    TROOT::IsZombie
    >>> comp.deactivate()
    >>> for suggestion in comp._completeImpl("TG"):
    ...     print(suggestion)
    '''

    def __init__(self):
        self.hook = None
        self.active = True
        self.firstActivation = True
        self.accessors = [".", "->", "::"]

    def activate(self):
        self.active = True
        if self.firstActivation:
            utils.declareCppCode('#include "dlfcn.h"')
            dlOpenRint = 'dlopen("libRint.so",RTLD_NOW);'
            utils.processCppCode(dlOpenRint)
            utils.declareCppCode(_TTabComHookCode)
            self.hook = ROOT._TTabComHook
            self.firstActivation = False

    def deactivate(self):
        self.active = False

    def _getSuggestions(self,line):
        if self.active:
            return self.hook(line)
        return []

    def _getLastAccessorPos(self,line):
        accessorPos = -1
        for accessor in self.accessors:
            tmpAccessorPos = line.rfind(accessor)
            if accessorPos < tmpAccessorPos:
                accessorPos = tmpAccessorPos+len(accessor)
        return accessorPos

    def _completeImpl(self, line):
        line=line.split()[-1]
        suggestions = [ str(s) for s in self._getSuggestions(line) ]
        suggestions = filter(lambda s: len(s.strip()) != 0, suggestions)
        suggestions = sorted(suggestions)
        if not suggestions: return []
        # Remove combinations of opening and closing brackets and just opening
        # brackets at the end of a line. Jupyter seems to expect functions
        # without these brackets to work properly. The brackets of'operator()'
        # must not be removed
        suggestions = [sugg[:-2] if sugg[-2:] == '()' and sugg != 'operator()' else sugg for sugg in suggestions]
        suggestions = [sugg[:-1] if sugg[-1:] == '(' else sugg for sugg in suggestions]
        # If a function signature is encountered, add an empty item to the
        # suggestions. Try to guess a function signature by an opening bracket
        # ignoring 'operator()'.
        are_signatures = "(" in "".join(filter(lambda s: s != 'operator()', suggestions))
        accessorPos = self._getLastAccessorPos(line)
        if are_signatures:
            suggestions = [" "] + suggestions
        elif accessorPos > 0:
            # Prepend variable name to suggestions. Do not prepend if the
            # suggestion already contains the variable name, this can happen if
            # e.g. there is only one valid completion
            if len(suggestions) > 1 or line[:accessorPos] != suggestions[0][:accessorPos]:
                suggestions = [line[:accessorPos]+sugg for sugg in suggestions]
        return suggestions

    def complete(self, ip, event) :
        '''
        Autocomplete interfacing to TTabCom. If an accessor of a scope is
        present in the line, the suggestions are prepended with the line.
        That's how completers work. For example:
        myGraph.Set<tab> will return "myGraph.Set+suggestion in the list of
        suggestions.
        '''
        return self._completeImpl(event.line)


_cppCompleter = CppCompleter()

def load_ipython_extension(ipython):
    _cppCompleter.activate()
    ipython.set_hook('complete_command', _cppCompleter.complete, re_key=r"[(.*)[\.,::,\->](.*)]|(.*)")

def unload_ipython_extension(ipython):
    _cppCompleter.deactivate()

