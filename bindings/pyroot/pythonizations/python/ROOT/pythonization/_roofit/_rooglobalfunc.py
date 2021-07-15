# Authors:
# * Jonas Rembser 06/2021
# * Harshal Shende 06/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################


r"""
/**
\class RooGlobalFunc
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly

## PyROOT

Some member functions of RooGlobalFunc that take a RooCmdArg as argument also support keyword arguments.
So far, this applies to FitOptions, Format, Frame, MultiArg, YVar and ZVar.
\code{.py}
# Directly passing a RooCmdArg:
ROOT.RooMCStudy(model, ROOT.RooArgSet(x), ROOT.RooFit.FitOptions=(ROOT.RooFit.Save(True), ROOT.RooFit.PrintEvalErrors(0)))

# With keyword arguments:
ROOT.RooMCStudy(model, ROOT.RooArgSet(x), FitOptions=dict(Save=True, PrintEvalErrors=0))
\endcode

\htmlonly
</div>
\endhtmlonly
*/
"""

from ._utils import _kwargs_to_roocmdargs, _string_to_root_attribute, _dict_to_std_map


# Color and Style dictionary to define matplotlib conventions
_color_map = {
    "r": "kRed",
    "b": "kBlue",
    "g": "kGreen",
    "y": "kYellow",
    "w": "kWhite",
    "k": "kBlack",
    "m": "kMagenta",
    "c": "kCyan",
}
_style_map = {"-": "kSolid", "--": "kDashed", ":": "kDotted", "-.": "kDashDotted"}


def FitOptions(*args, **kwargs):
    # Redefinition of `FitOptions` for keyword arguments.
    # The keywords must correspond to the CmdArg of the `FitOptions` function.
    from cppyy.gbl import RooFit

    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._FitOptions(*args, **kwargs)


def Format(*args, **kwargs):
    # Redefinition of `Format` for keyword arguments.
    # The keywords must correspond to the CmdArg of the `Format` function.
    from cppyy.gbl import RooFit

    if "what" in kwargs:
        args = (kwargs["what"],) + args
        del kwargs["what"]
    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._Format(*args, **kwargs)


def Frame(*args, **kwargs):
    # Redefinition of `Frame` for keyword arguments.
    # The keywords must correspond to the CmdArg of the `Frame` function.
    from cppyy.gbl import RooFit

    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._Frame(*args, **kwargs)


def MultiArg(*args, **kwargs):
    # Redefinition of `MultiArg` for keyword arguments.
    # The keywords must correspond to the CmdArg of the `MultiArg` function.
    from cppyy.gbl import RooFit

    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._MultiArg(*args, **kwargs)


def YVar(*args, **kwargs):
    # Redefinition of `YVar` for keyword arguments.
    # The keywords must correspond to the CmdArg of the `YVar` function.
    from cppyy.gbl import RooFit

    if "var" in kwargs:
        args = (kwargs["var"],) + args
        del kwargs["var"]
    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._YVar(*args, **kwargs)


def ZVar(*args, **kwargs):
    # Redefinition of `ZVar` for keyword arguments.
    # The keywords must correspond to the CmdArg of the `ZVar` function.
    from cppyy.gbl import RooFit

    if "var" in kwargs:
        args = (kwargs["var"],) + args
        del kwargs["var"]
    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._ZVar(*args, **kwargs)


def Slice(*args, **kwargs):
    # Redefinition of `Slice` for keyword arguments and converting python dict to std::map.
    # The keywords must correspond to the CmdArg of the `Slice` function.
    # The instances in the dict must correspond to the template argument in std::map of the `Slice` function.
    from cppyy.gbl import RooFit

    if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], dict):
        args = list(args)
        args[0] = _dict_to_std_map(args[0])
        return RooFit._Slice(args[0])

    return RooFit._Slice(*args, **kwargs)


def Import(*args, **kwargs):
    # Redefinition of `Import` for keyword arguments and converting python dict to std::map.
    # The keywords must correspond to the CmdArg of the `Import` function.
    # The instances in the dict must correspond to the template argument in std::map of the `Import` function.
    from cppyy.gbl import RooFit

    if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], dict):
        args = list(args)
        args[0] = _dict_to_std_map(args[0])
        return RooFit._Import(args[0])

    return RooFit._Import(*args, **kwargs)


def Link(*args, **kwargs):
    # Redefinition of `Link` for keyword arguments and converting python dict to std::map.
    # The keywords must correspond to the CmdArg of the `Link` function.
    # The instances in the dict must correspond to the template argument in std::map of the `Link` function.
    from cppyy.gbl import RooFit

    if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], dict):
        args = list(args)
        args[0] = _dict_to_std_map(args[0])
        return RooFit._Import(args[0])

    return RooFit._Link(*args, **kwargs)


def LineColor(color):
    # Redefinition of `LineColor` for matplotlib conventions and string arguments.
    from cppyy.gbl import RooFit

    return RooFit._LineColor(_string_to_root_attribute(color, _color_map))


def FillColor(color):
    # Redefinition of `FillColor` for matplotlib conventions and string arguments.
    from cppyy.gbl import RooFit

    return RooFit._FillColor(_string_to_root_attribute(color, _color_map))


def MarkerColor(color):
    # Redefinition of `MarkerColor` for matplotlib conventions and string arguments.
    from cppyy.gbl import RooFit

    return RooFit._MarkerColor(_string_to_root_attribute(color, _color_map))


def LineStyle(style):
    # Redefinition of `LineStyle` for matplotlib conventions and string arguments.
    from cppyy.gbl import RooFit

    return RooFit._LineStyle(_string_to_root_attribute(style, _style_map))


def FillStyle(style):
    # Redefinition of `FillStyle` for matplotlib conventions and string arguments.
    from cppyy.gbl import RooFit

    return RooFit._FillStyle(_string_to_root_attribute(style, {}))


def MarkerStyle(style):
    # Redefinition of `MarkerStyle` for matplotlib conventions and string arguments.
    from cppyy.gbl import RooFit

    return RooFit._MarkerStyle(_string_to_root_attribute(style, {}))
