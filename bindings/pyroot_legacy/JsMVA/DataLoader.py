# -*- coding: utf-8 -*-
## @package JsMVA.DataLoader
# @author Attila Bagoly <battila93@gmail.com>
# DataLoader module with the functions to be inserted to TMVA::DataLoader class and helper functions


from ROOT import TH1F, TMVA, TBufferJSON
import sys
if sys.version_info >= (3, 0):
    from JsMVA import JPyInterface
else:
    import JPyInterface
from JsMVA.Utils import xrange
import ROOT


## Creates the input variable histogram and perform the transformations if necessary
# @param dl DataLoader object
# @param className string Signal/Background
# @param variableName string containing the variable name
# @param numBin for creating the histogram
# @param processTrfs string containing the list of transformations to be used on input variable; eg. "I;N;D;P;U;G,D"
def GetInputVariableHist(dl, className, variableName, numBin, processTrfs=""):
    dsinfo = dl.GetDefaultDataSetInfo()
    vi = 0
    ivar = 0
    for i in range(dsinfo.GetNVariables()):
        if dsinfo.GetVariableInfo(i).GetLabel()==variableName:
            vi   = dsinfo.GetVariableInfo(i)
            ivar = i
            break
    if vi==0:
        return 0

    h = TH1F(className, str(vi.GetExpression()) + " ("+className+")", numBin, vi.GetMin(), vi.GetMax())

    clsn = dsinfo.GetClassInfo(className).GetNumber()
    ds   = dsinfo.GetDataSet()

    trfsDef = processTrfs.split(';')
    trfs    = []
    for trfDef in trfsDef:
        trfs.append(TMVA.TransformationHandler(dsinfo, "DataLoader"))
        TMVA.CreateVariableTransforms( trfDef, dsinfo, trfs[-1], dl.Log())

    inputEvents = ds.GetEventCollection()
    transformed = 0
    tmp         = 0
    for trf in trfs:
        if transformed==0:
            transformed = trf.CalcTransformations(inputEvents, 1)
        else:
            tmp = trf.CalcTransformations(transformed, 1)
            transformed = tmp

    if transformed!=0:
        for event in transformed:
            if event.GetClass() != clsn:
                continue
            h.Fill(event.GetValue(ivar))
    else:
        for event in inputEvents:
            if event.GetClass() != clsn:
                continue
            h.Fill(event.GetValue(ivar))
    return (h)


## Get correlation matrix in JSON format
# This function is used by OutputTransformer
# @param dl the object pointer
# @param className Signal/Background
def GetCorrelationMatrixInJSON(className, varNames, matrix):
    m = ROOT.TMatrixD(len(matrix), len(matrix))
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix)):
            m[i][j] = matrix[i][j]
    th2 = ROOT.TH2D(m)
    th2.SetTitle("Correlation matrix ("+className+")")
    for i in xrange(len(varNames)):
        th2.GetXaxis().SetBinLabel(i+1, varNames[i])
        th2.GetYaxis().SetBinLabel(i+1, varNames[i])
    th2.Scale(100.0)
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix)):
            th2.SetBinContent(i+1, j+1, int(th2.GetBinContent(i+1, j+1)))
    th2.SetStats(0)
    th2.SetMarkerSize(1.5)
    th2.SetMarkerColor(0)
    labelSize = 0.040
    th2.GetXaxis().SetLabelSize(labelSize)
    th2.GetYaxis().SetLabelSize(labelSize)
    th2.LabelsOption("d")
    th2.SetLabelOffset(0.011)
    th2.SetMinimum(-100.0)
    th2.SetMaximum(+100.0)
    dat = TBufferJSON.ConvertToJSON(th2)
    return str(dat).replace("\n", "")

## Draw correlation matrix
# This function uses the TMVA::DataLoader::GetCorrelationMatrix function added newly to root
# @param dl the object pointer
# @param className Signal/Background
def DrawCorrelationMatrix(dl, className):
    th2 = dl.GetCorrelationMatrix(className)
    th2.SetMarkerSize(1.5)
    th2.SetMarkerColor(0)
    labelSize = 0.040
    th2.GetXaxis().SetLabelSize(labelSize)
    th2.GetYaxis().SetLabelSize(labelSize)
    th2.LabelsOption("d")
    th2.SetLabelOffset(0.011)
    JPyInterface.JsDraw.Draw(th2, 'drawTH2')

## Draw input variables
# This function uses the previously defined GetInputVariableHist function to create the histograms
# @param dl The object pointer
# @param variableName string containing the variable name
# @param numBin for creating the histogram
# @param processTrfs list of transformations to be used on input variable; eg. ["I", "N", "D", "P", "U", "G"]"
def DrawInputVariable(dl, variableName, numBin=100, processTrfs=[]):
    processTrfsSTR = ""
    if len(processTrfs)>0:
        for o in processTrfs:
            processTrfsSTR += str(o) + ";"
        processTrfsSTR = processTrfsSTR[:-1]
    sig = GetInputVariableHist(dl, "Signal",     variableName, numBin, processTrfsSTR)
    bkg = GetInputVariableHist(dl, "Background", variableName, numBin, processTrfsSTR)
    c, l = JPyInterface.JsDraw.sbPlot(sig, bkg, {"xaxis": sig.GetTitle(),
                                    "yaxis": "Number of events",
                                    "plot": "Input variable: "+sig.GetTitle()})
    JPyInterface.JsDraw.Draw(c)

## Rewrite TMVA::DataLoader::PrepareTrainingAndTestTree
# @param *args positional parameters
# @param **kwargs named parameters: this will be transformed to option string
def ChangeCallOriginalPrepareTrainingAndTestTree(*args, **kwargs):
    if len(kwargs)==0:
        originalFunction, args = JPyInterface.functions.ProcessParameters(0, *args, **kwargs)
        return originalFunction(*args)
    try:
        args, kwargs = JPyInterface.functions.ConvertSpecKwargsToArgs(["SigCut", "BkgCut"], *args, **kwargs)
    except AttributeError:
        try:
            args, kwargs = JPyInterface.functions.ConvertSpecKwargsToArgs(["Cut"], *args, **kwargs)
        except AttributeError:
            raise AttributeError
    originalFunction, args = JPyInterface.functions.ProcessParameters(3, *args, **kwargs)
    return originalFunction(*args)
