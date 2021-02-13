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

"""This preprocessor marks cell's metadata so that the appropriate
highlighter can be used in the `highlight` filter afterwards.
More precisely, the language of a cell is set to C++ in two scenarios:
- Python notebooks: cells with `%%cpp` magic extension.
- ROOT C++ notebooks: all cells.
This preprocessor relies on the metadata section of the notebook to
find out about the notebook's language.
"""

from IPython.nbconvert.preprocessors.base import Preprocessor
import re


class CppHighlighter(Preprocessor):
    """
    Detects and tags code cells that use the C++ language.
    """

    magics = [ '%%cpp' ]
    cpp = 'cpp'
    python = 'python'

    def __init__(self, config=None, **kw):
        super(CppHighlighter, self).__init__(config=config, **kw)

        # Build regular expressions to catch language extensions or switches and choose
        # an adequate pygments lexer
        any_magic = "|".join(self.magics)
        self.re_magic_language = re.compile(r"^\s*({0}).*".format(any_magic), re.DOTALL)

    def matches(self, source, reg_expr):
        m = reg_expr.match(source)
        if m:
            return True
        else:
            return False

    def _preprocess_cell_python(self, cell, resources, cell_index):
        # Mark %%cpp and %%dcl code cells as cpp
        if cell.cell_type == "code" and self.matches(cell.source, self.re_magic_language):
            cell['metadata']['magics_language'] = self.cpp

        return cell, resources

    def _preprocess_cell_cpp(self, cell, resources, cell_index):
        # Mark all code cells as cpp
        if cell.cell_type == "code":
            cell['metadata']['magics_language'] = self.cpp

        return cell, resources

    def preprocess(self, nb, resources):
        self.preprocess_cell = self._preprocess_cell_python
        try:
            if nb.metadata.kernelspec.language == "c++":
                self.preprocess_cell = self._preprocess_cell_cpp
        except:
            # if no language metadata, keep python as default
            pass
        return super(CppHighlighter, self).preprocess(nb, resources)
