# -*- coding: utf-8 -*-
## @package JsMVA.Factory
#  @author  Attila Bagoly <battila93@gmail.com>
# Factory module with the functions to be inserted to TMVA::Factory class and helper functions and classes


import ROOT
from ROOT import TMVA
import sys
if sys.version_info >= (3, 0):
    from JsMVA import JPyInterface
else:
    import JPyInterface
from JsMVA.Utils import xrange
from xml.etree.ElementTree import ElementTree
import json
from IPython.core.display import display, HTML, clear_output
from ipywidgets import widgets
from threading import Thread
import time
from string import Template


# This class contains the necessary HTML, JavaScript, CSS codes (templates)
# for the new Factory methods. Some parts of these variables will be replaced and the new string will be the cell output.
class __HTMLJSCSSTemplates:
    # stop button
    button = """
     <script type="text/javascript">
     require(["jquery"], function(jQ){
         jQ("input.stopTrainingButton").on("click", function(){
             IPython.notebook.kernel.interrupt();
             jQ(this).css({
                 "background-color": "rgba(200, 0, 0, 0.8)",
                 "color": "#fff",
                 "box-shadow": "0 3px 5px rgba(0, 0, 0, 0.3)",
             });
         });
     });
     </script>
     <style type="text/css">
     input.stopTrainingButton {
         background-color: #fff;
         border: 1px solid #ccc;
         width: 100%;
         font-size: 16px;
         font-weight: bold;
         padding: 6px 12px;
         cursor: pointer;
         border-radius: 6px;
         color: #333;
     }
     input.stopTrainingButton:hover {
         background-color: rgba(204, 204, 204, 0.4);
     }
     </style>
     <input type="button" value="Stop" class="stopTrainingButton" />
     """
    # progress bar
    inc = Template("""
     <script type="text/javascript" id="progressBarScriptInc">
     require(["jquery"], function(jQ){
         jQ("#jsmva_bar_$id").css("width", $progress + "%");
         jQ("#jsmva_label_$id").text($progress + '%');
         jQ("#progressBarScriptInc").parent().parent().remove();
     });
     </script>
     """)
    progress_bar = Template("""
     <style>
     #jsmva_progress_$id {
         position: relative;
         float: left;
         height: 30px;
         width: 100%;
         background-color: #f5f5f5;
         border-radius: 3px;
         box-shadow: inset 0 3px 6px rgba(0, 0, 0, 0.1);
     }
     #jsmva_bar_$id {
         position: absolute;
         width: 1%;
         height: 100%;
         background-color: #337ab7;
     }
     #jsmva_label_$id {
         text-align: center;
         line-height: 30px;
         color: white;
     }
     </style>
     <div id="jsmva_progress_$id">
       <div id="jsmva_bar_$id">
         <div id="jsmva_label_$id">0%</div>
       </div>
     </div>
     """)


## Getting method object from factory
# @param fac the TMVA::Factory object
# @param datasetName selecting the dataset
# @param methodName which method we want to get
def GetMethodObject(fac, datasetName, methodName):
    method = []
    for methodMapElement in fac.fMethodsMap:
        if methodMapElement[0] != datasetName:
            continue
        methods = methodMapElement[1]
        for m in methods:
            m.GetName._threaded = True
            if m.GetName() == methodName:
                method.append( m )
                break
    if len(method) != 1:
        print("Factory.GetMethodObject: no method object found")
        return None
    return (method[0])

## Reads deep neural network weights from file and returns it in JSON format. READER FOR OLD XML STRUCTURE OF DNN.
# IT WILL BE REMOVED IN THE FUTURE.
# @param xml_file path to DNN weight file
# @param returnObj if False it will return a JSON string, if True it will return the JSON object itself
def GetDeepNetworkOld(xml_file, returnObj=False):
    tree = ElementTree()
    tree.parse(xml_file)
    roottree = tree.getroot()
    network  = {}
    network["variables"] = []
    for v in roottree.find("Variables"):
        network["variables"].append(v.get('Title'))
    layout = roottree.find("Weights").find("Layout")
    net    = []
    for layer in layout:
        net.append({"Connection": layer.get("Connection"),
        "Nodes": layer.get("Nodes"),
        "ActivationFunction": layer.get("ActivationFunction"),
        "OutputMode": layer.get("OutputMode")
         })
    network["layers"] = net
    Synapses = roottree.find("Weights").find("Synapses")
    synapses = {
        "InputSize": Synapses.get("InputSize"),
        "OutputSize": Synapses.get("OutputSize"),
        "NumberSynapses": Synapses.get("NumberSynapses"),
        "synapses": []
    }
    for i in Synapses.text.split(" "):
        tmp = i.replace("\n", "")
        if len(tmp)>1:
            synapses["synapses"].append(tmp)
    network["synapses"] = synapses
    if returnObj:
        return network
    return json.dumps(network, sort_keys = True)

## Reads deep neural network weights from file and returns it in JSON format.
# @param xml_file path to DNN weight file
# @param returnObj if False it will return a JSON string, if True it will return the JSON object itself
def GetDeepNetwork(xml_file, returnObj=False):
    tree = ElementTree()
    tree.parse(xml_file)
    roottree = tree.getroot()
    network  = {}
    network["variables"] = []
    for v in roottree.find("Variables"):
        network["variables"].append(v.get('Title'))
    layers = []
    for lyr in roottree.find("Weights"):
        weights = lyr.find("Weights")
        biases = lyr.find("Biases")
        wd = []
        for i in weights.text.split(" "):
            tmp = i.replace("\n", "")
            if len(tmp)>=1:
                wd.append(float(tmp))
        bd = []
        for i in biases.text.split(" "):
            tmp = i.replace("\n", "")
            if len(tmp)>=1:
                bd.append(float(tmp))
        layer = {
            "ActivationFunction": lyr.get("ActivationFunction"),
            "Weights": {
                "row": int(weights.get("rows")),
                "cols": int(weights.get("cols")),
                "data": wd
            },
            "Biases": {
                "row": int(biases.get("rows")),
                "cols": int(biases.get("cols")),
                "data": bd
            }
        }
        layers.append(layer)
    network["layers"] = layers
    if returnObj:
        return network
    return json.dumps(network, sort_keys = True)

## Reads neural network weights from file and returns it in JSON format
# @param xml_file path to weight file
def GetNetwork(xml_file):
    tree = ElementTree()
    tree.parse(xml_file)
    roottree = tree.getroot()
    network  = {}
    network["variables"] = []
    for v in roottree.find("Variables"):
        network["variables"].append(v.get('Title'))
    layout = roottree.find("Weights").find("Layout")

    net    = { "nlayers": layout.get("NLayers") }
    for layer in layout:
        neuron_num = int(layer.get("NNeurons"))
        neurons    = { "nneurons": neuron_num }
        i = 0
        for neuron in layer:
            label = "neuron_"+str(i)
            i    += 1
            nsynapses = int(neuron.get('NSynapses'))
            neurons[label] = {"nsynapses": nsynapses}
            if nsynapses == 0:
                break
            text = str(neuron.text)
            wall = text.replace("\n","").split(" ")
            weights = []
            for w in wall:
                if w!="":
                    weights.append(float(w))
            neurons[label]["weights"] = weights
        net["layer_"+str(layer.get('Index'))] = neurons
    network["layout"] = net
    return json.dumps(network, sort_keys = True)

## Helper class for reading decision tree from XML file
class TreeReader:

    ## Standard Constructor
    # @param self object pointer
    # @oaran fileName path to XML file
    def __init__(self, fileName):
        self.__xmltree = ElementTree()
        self.__xmltree.parse(fileName)
        self.__NTrees = int(self.__xmltree.find("Weights").get('NTrees'))

    ## Returns the number of trees
    # @param self object pointer
    def getNTrees(self):
        return (self.__NTrees)

    # Returns DOM object to selected tree
    # @param self object pointer
    # @param itree the index of tree
    def __getBinaryTree(self, itree):
        if self.__NTrees<=itree:
            print( "to big number, tree number must be less then %s"%self.__NTrees )
            return 0
        return self.__xmltree.find("Weights").find("BinaryTree["+str(itree+1)+"]")

    ## Reads the tree
    # @param self the object pointer
    # @param binaryTree the tree DOM object to be read
    # @param tree empty object, this will be filled
    # @param depth current depth
    def __readTree(self, binaryTree, tree={}, depth=0):
        nodes = binaryTree.findall("Node")
        if len(nodes)==0:
            return
        if len(nodes)==1 and nodes[0].get("pos")=="s":
            info = {
                "IVar":   nodes[0].get("IVar"),
                "Cut" :   nodes[0].get("Cut"),
                "purity": nodes[0].get("purity"),
                "pos":    0
            }
            tree["info"]     = info
            tree["children"] = []
            self.__readTree(nodes[0], tree, 1)
            return
        for node in nodes:
            info = {
                "IVar":   node.get("IVar"),
                "Cut" :   node.get("Cut"),
                "purity": node.get("purity"),
                "pos":    node.get("pos")
            }
            tree["children"].append({
               "info": info,
                "children": []
            })
            self.__readTree(node, tree["children"][-1], depth+1)

    ## Public function which returns the specified tree object
    # @param self the object pointer
    # @param itree selected tree index
    def getTree(self, itree):
        binaryTree = self.__getBinaryTree(itree)
        if binaryTree==0:
            return {}
        tree = {}
        self.__readTree(binaryTree, tree)
        return tree

    ## Returns a list with input variable names
    # @param self the object pointer
    def getVariables(self):
        varstree = self.__xmltree.find("Variables").findall("Variable")
        variables = [None]*len(varstree)
        for v in varstree:
            variables[int(v.get('VarIndex'))] = v.get('Expression')
        return variables


## Draw ROC curve
# @param fac the object pointer
# @param datasetName the dataset name
def DrawROCCurve(fac, datasetName):
    canvas = fac.GetROCCurve(datasetName)
    JPyInterface.JsDraw.Draw(canvas)

## Draw output distributions
# @param fac the object pointer
# @param datasetName the dataset name
# @param methodName we want to see the output distribution of this method
def DrawOutputDistribution(fac, datasetName, methodName):
    method = GetMethodObject(fac, datasetName, methodName)
    if method==None:
        return None
    mvaRes = method.Data().GetResults(method.GetMethodName(), TMVA.Types.kTesting, TMVA.Types.kMaxAnalysisType)
    sig    = mvaRes.GetHist("MVA_S")
    bgd    = mvaRes.GetHist("MVA_B")
    c, l = JPyInterface.JsDraw.sbPlot(sig, bgd, {"xaxis": methodName+" response",
                                    "yaxis": "(1/N) dN^{ }/^{ }dx",
                                    "plot": "TMVA response for classifier: "+methodName})
    JPyInterface.JsDraw.Draw(c)

## Draw output probability distribution
# @param fac the object pointer
# @param datasetName the dataset name
# @param methodName we want to see the output probability distribution of this method
def DrawProbabilityDistribution(fac, datasetName, methodName):
    method = GetMethodObject(fac, datasetName, methodName)
    if method==0:
        return
    mvaRes = method.Data().GetResults(method.GetMethodName(), TMVA.Types.kTesting, TMVA.Types.kMaxAnalysisType)
    sig    = mvaRes.GetHist("Prob_S")
    bgd    = mvaRes.GetHist("Prob_B") #Rar_S
    c, l   = JPyInterface.JsDraw.sbPlot(sig, bgd, {"xaxis": "Signal probability",
                                        "yaxis": "(1/N) dN^{ }/^{ }dx",
                                        "plot": "TMVA probability for classifier: "+methodName})
    JPyInterface.JsDraw.Draw(c)

## Draw cut efficiencies
# @param fac the object pointer
# @param datasetName the dataset name
# @param methodName we want to see the cut efficiencies of this method
def DrawCutEfficiencies(fac, datasetName, methodName):
    #reading histograms
    method = GetMethodObject(fac, datasetName, methodName)
    if method==0:
        return
    mvaRes = method.Data().GetResults(method.GetMethodName(), TMVA.Types.kTesting, TMVA.Types.kMaxAnalysisType)
    sigE = mvaRes.GetHist("MVA_EFF_S")
    bgdE = mvaRes.GetHist("MVA_EFF_B")

    fNSignal = 1000
    fNBackground = 1000

    f = ROOT.TFormula("sigf", "x/sqrt(x+y)")

    pname    = "purS_"         + methodName
    epname   = "effpurS_"      + methodName
    ssigname = "significance_" + methodName

    nbins = sigE.GetNbinsX()
    low   = sigE.GetBinLowEdge(1)
    high  = sigE.GetBinLowEdge(nbins+1)

    purS    = ROOT.TH1F(pname, pname, nbins, low, high)
    sSig    = ROOT.TH1F(ssigname, ssigname, nbins, low, high)
    effpurS = ROOT.TH1F(epname, epname, nbins, low, high)

    # formating the style of histograms
    #chop off useless stuff
    sigE.SetTitle( "Cut efficiencies for "+methodName+" classifier")

    TMVA.TMVAGlob.SetSignalAndBackgroundStyle( sigE, bgdE )
    TMVA.TMVAGlob.SetSignalAndBackgroundStyle( purS, bgdE )
    TMVA.TMVAGlob.SetSignalAndBackgroundStyle( effpurS, bgdE )
    sigE.SetFillStyle( 0 )
    bgdE.SetFillStyle( 0 )
    sSig.SetFillStyle( 0 )
    sigE.SetLineWidth( 3 )
    bgdE.SetLineWidth( 3 )
    sSig.SetLineWidth( 3 )

    purS.SetFillStyle( 0 )
    purS.SetLineWidth( 2 )
    purS.SetLineStyle( 5 )
    effpurS.SetFillStyle( 0 )
    effpurS.SetLineWidth( 2 )
    effpurS.SetLineStyle( 6 )
    sig = 0
    maxSigErr = 0
    for i in range(1,sigE.GetNbinsX()+1):
        eS = sigE.GetBinContent( i )
        S = eS * fNSignal
        B = bgdE.GetBinContent( i ) * fNBackground
        if (S+B)==0:
            purS.SetBinContent( i, 0)
        else:
            purS.SetBinContent( i, S/(S+B) )

        sSig.SetBinContent( i, f.Eval(S,B) )
        effpurS.SetBinContent( i, eS*purS.GetBinContent( i ) )

    maxSignificance = sSig.GetMaximum()
    maxSignificanceErr = 0
    sSig.Scale(1/maxSignificance)

    c = ROOT.TCanvas( "canvasCutEff","Cut efficiencies for "+methodName+" classifier", JPyInterface.JsDraw.jsCanvasWidth,
                      JPyInterface.JsDraw.jsCanvasHeight)

    c.SetGrid(1)
    c.SetTickx(0)
    c.SetTicky(0)

    TMVAStyle = ROOT.gROOT.GetStyle("Plain")
    TMVAStyle.SetLineStyleString( 5, "[32 22]" )
    TMVAStyle.SetLineStyleString( 6, "[12 22]" )

    c.SetTopMargin(.2)

    effpurS.SetTitle("Cut efficiencies and optimal cut value")
    if methodName.find("Cuts")!=-1:
        effpurS.GetXaxis().SetTitle( "Signal Efficiency" )
    else:
        effpurS.GetXaxis().SetTitle( "Cut value applied on " + methodName + " output" )
    effpurS.GetYaxis().SetTitle( "Efficiency (Purity)" )
    TMVA.TMVAGlob.SetFrameStyle( effpurS )

    c.SetTicks(0,0)
    c.SetRightMargin ( 2.0 )

    effpurS.SetMaximum(1.1)
    effpurS.Draw("histl")

    purS.Draw("samehistl")

    sigE.Draw("samehistl")
    bgdE.Draw("samehistl")

    signifColor = ROOT.TColor.GetColor( "#00aa00" )

    sSig.SetLineColor( signifColor )
    sSig.Draw("samehistl")

    effpurS.Draw( "sameaxis" )

    #Adding labels and other informations to plots.

    legend1 = ROOT.TLegend( c.GetLeftMargin(), 1 - c.GetTopMargin(),
                                     c.GetLeftMargin() + 0.4, 1 - c.GetTopMargin() + 0.12 )
    legend1.SetFillStyle( 1 )
    legend1.AddEntry(sigE,"Signal efficiency","L")
    legend1.AddEntry(bgdE,"Background efficiency","L")
    legend1.Draw("same")
    legend1.SetBorderSize(1)
    legend1.SetMargin( 0.3 )


    legend2 = ROOT.TLegend( c.GetLeftMargin() + 0.4, 1 - c.GetTopMargin(),
                                     1 - c.GetRightMargin(), 1 - c.GetTopMargin() + 0.12 )
    legend2.SetFillStyle( 1 )
    legend2.AddEntry(purS,"Signal purity","L")
    legend2.AddEntry(effpurS,"Signal efficiency*purity","L")
    legend2.AddEntry(sSig, "S/#sqrt{ S+B }","L")
    legend2.Draw("same")
    legend2.SetBorderSize(1)
    legend2.SetMargin( 0.3 )

    effline = ROOT.TLine( sSig.GetXaxis().GetXmin(), 1, sSig.GetXaxis().GetXmax(), 1 )
    effline.SetLineWidth( 1 )
    effline.SetLineColor( 1 )
    effline.Draw()

    c.Update()

    tl = ROOT.TLatex()
    tl.SetNDC()
    tl.SetTextSize( 0.033 )
    maxbin = sSig.GetMaximumBin()
    line1 = tl.DrawLatex( 0.15, 0.23, "For %1.0f signal and %1.0f background"%(fNSignal, fNBackground))
    tl.DrawLatex( 0.15, 0.19, "events the maximum S/#sqrt{S+B} is")

    if maxSignificanceErr > 0:
        line2 = tl.DrawLatex( 0.15, 0.15, "%5.2f +- %4.2f when cutting at %5.2f"%(
                                                      maxSignificance,
                                                      maxSignificanceErr,
                                                      sSig.GetXaxis().GetBinCenter(maxbin)) )
    else:
        line2 = tl.DrawLatex( 0.15, 0.15, "%4.2f when cutting at %5.2f"%(
                                                      maxSignificance,
                                                      sSig.GetXaxis().GetBinCenter(maxbin)) )

    if methodName.find("Cuts")!=-1:
        tl.DrawLatex( 0.13, 0.77, "Method Cuts provides a bundle of cut selections, each tuned to a")
        tl.DrawLatex(0.13, 0.74, "different signal efficiency. Shown is the purity for each cut selection.")

    wx = (sigE.GetXaxis().GetXmax()+abs(sigE.GetXaxis().GetXmin()))*0.135
    rightAxis = ROOT.TGaxis( sigE.GetXaxis().GetXmax()+wx,
                             c.GetUymin()-0.3,
                             sigE.GetXaxis().GetXmax()+wx,
                             0.7, 0, 1.1*maxSignificance,510, "+L")
    rightAxis.SetLineColor ( signifColor )
    rightAxis.SetLabelColor( signifColor )
    rightAxis.SetTitleColor( signifColor )

    rightAxis.SetTitleSize( sSig.GetXaxis().GetTitleSize() )
    rightAxis.SetTitle( "Significance" )
    rightAxis.Draw()

    c.Update()

    JPyInterface.JsDraw.Draw(c)

## Draw neural network
# @param fac the object pointer
# @param datasetName the dataset name
# @param methodName we want to see the network created by this method
def DrawNeuralNetwork(fac, datasetName, methodName):
    m = GetMethodObject(fac, datasetName, methodName)
    if m==None:
        return None
    if m.GetMethodType() == ROOT.TMVA.Types.kDNN:
        try:
            net = GetDeepNetwork(str(m.GetWeightFileName()))
        except AttributeError:
            net = GetDeepNetworkOld(str(m.GetWeightFileName()))
    else:
        net = GetNetwork(str(m.GetWeightFileName()))
    JPyInterface.JsDraw.Draw(net, "drawNeuralNetwork", True)

## Draw deep neural network
# @param fac the object pointer
# @param datasetName the dataset name
# @param methodName we want to see the deep network created by this method
def DrawDecisionTree(fac, datasetName, methodName):
    m = GetMethodObject(fac, datasetName, methodName)
    if m==None:
        return None
    tr = TreeReader(str(m.GetWeightFileName()))

    variables = tr.getVariables();

    def clicked(b):
        if treeSelector.value>tr.getNTrees():
            treeSelector.value = tr.getNTrees()
        clear_output()
        toJs = {
            "variables": variables,
            "tree": tr.getTree(treeSelector.value)
        }
        json_str = json.dumps(toJs)
        JPyInterface.JsDraw.Draw(json_str, "drawDecisionTree", True)

    mx = str(tr.getNTrees()-1)

    treeSelector = widgets.IntText(value=0, font_weight="bold")
    drawTree     = widgets.Button(description="Draw", font_weight="bold")
    label        = widgets.HTML("<div style='padding: 6px;font-weight:bold;color:#333;'>Decision Tree [0-"+mx+"]:</div>")

    drawTree.on_click(clicked)
    container = widgets.HBox([label,treeSelector, drawTree])
    display(container)

## This function puts the main thread to sleep until data points for tracking plots appear.
# @param m Method object
# @param sleep_time default sleeping time
def GotoSleepUntilTrackingReady(m, sleep_time):
    sleep_index = 1
    while m.GetCurrentIter() == 0 and not m.TrainingEnded():
        time.sleep(sleep_time * sleep_index)
        grs = m.GetInteractiveTrainingError().GetListOfGraphs()
        if grs and len(grs) > 1 and grs[0].GetN() > 1:
            break
        if sleep_index < 120:
            sleep_index += 1
    if sleep_index > 2:
        sleep_time = float(sleep_time * sleep_index / 2)
    return sleep_time

## Rewrite function for TMVA::Factory::TrainAllMethods. This function provides interactive training.
# The training will be started on separated thread. The main thread will periodically check for updates and will create
# the JS output which will update the plots and progress bars. The main thread must contain `while True`, because, if not
# it will cause crash (output will be flushed by tornado IOLoop (runs on main thread), but the output streams are
# C++ atomic types)
# @param fac the factory object pointer
def ChangeTrainAllMethods(fac):
    clear_output()
    #stop button
    button = __HTMLJSCSSTemplates.button
    #progress bar
    inc = __HTMLJSCSSTemplates.inc
    progress_bar = __HTMLJSCSSTemplates.progress_bar
    progress_bar_idx = 0
    TTypes = ROOT.TMVA.Types
    error_plot_supported = [int(TTypes.kMLP), int(TTypes.kDNN), int(TTypes.kBDT)]
    exit_button_supported = [int(TTypes.kSVM), int(TTypes.kCuts), int(TTypes.kBoost), int(TTypes.kBDT)]

    for methodMapElement in fac.fMethodsMap:
        sleep_time = 0.5
        display(HTML("<center><h1>Dataset: "+str(methodMapElement[0])+"</h1></center>"))
        for m in methodMapElement[1]:
            m.GetMethodType._threaded = True
            m.GetName._threaded = True
            method_type = int(m.GetMethodType())
            name = str(m.GetName())
            display(HTML("<h2><b>Train method: "+name+"</b></h2>"))
            m.InitIPythonInteractive()
            t = Thread(target=ROOT.TMVA.MethodBase.TrainMethod, args=[m])
            t.start()
            if method_type in error_plot_supported:
                time.sleep(sleep_time)
                sleep_time = GotoSleepUntilTrackingReady(m, sleep_time)
                display(HTML(button))
                if m.GetMaxIter() != 0:
                    display(HTML(progress_bar.substitute({"id": progress_bar_idx})))
                    display(HTML(inc.substitute({"id": progress_bar_idx, "progress": 100 * m.GetCurrentIter() / m.GetMaxIter()})))
                JPyInterface.JsDraw.Draw(m.GetInteractiveTrainingError(), "drawTrainingTestingErrors")
                try:
                    while not m.TrainingEnded():
                        JPyInterface.JsDraw.InsertData(m.GetInteractiveTrainingError())
                        if m.GetMaxIter() != 0:
                            display(HTML(inc.substitute({
                                "id": progress_bar_idx,
                                "progress": 100 * m.GetCurrentIter() / m.GetMaxIter()
                            })))
                        time.sleep(sleep_time)
                except KeyboardInterrupt:
                    m.ExitFromTraining()
            else:
                if method_type in exit_button_supported:
                    display(HTML(button))
                time.sleep(sleep_time)
                if m.GetMaxIter()!=0:
                    display(HTML(progress_bar.substitute({"id": progress_bar_idx})))
                    display(HTML(inc.substitute({"id": progress_bar_idx, "progress": 100*m.GetCurrentIter()/m.GetMaxIter()})))
                else:
                    display(HTML("<b>Training...</b>"))
                if method_type in exit_button_supported:
                    try:
                        while not m.TrainingEnded():
                            if m.GetMaxIter()!=0:
                                display(HTML(inc.substitute({
                                    "id": progress_bar_idx,
                                    "progress": 100 * m.GetCurrentIter() / m.GetMaxIter()
                                })))
                            time.sleep(sleep_time)
                    except KeyboardInterrupt:
                        m.ExitFromTraining()
                else:
                    while not m.TrainingEnded():
                        if m.GetMaxIter() != 0:
                            display(HTML(inc.substitute({
                                "id": progress_bar_idx,
                                "progress": 100 * m.GetCurrentIter() / m.GetMaxIter()
                            })))
                        time.sleep(sleep_time)
            if m.GetMaxIter() != 0:
                display(HTML(inc.substitute({
                    "id": progress_bar_idx,
                    "progress": 100 * m.GetCurrentIter() / m.GetMaxIter()
                })))
            else:
                display(HTML("<b>End</b>"))
            progress_bar_idx += 1
            t.join()
    return

## Rewrite the constructor of TMVA::Factory
# @param *args positional parameters
# @param **kwargs named parameters: this will be transformed to option string
def ChangeCallOriginal__init__(*args,  **kwargs):
    hasColor = False
    args = list(args)
    for arg_idx in xrange(len(args)):
        # basestring==(str, unicode) in Python2, which translates to str in Python3
        if sys.version_info >= (3, 0):
            is_string = isinstance(args[arg_idx], str)
        else:
            is_string = isinstance(args[arg_idx], basestring)
        if is_string and args[arg_idx].find(":")!=-1:
            if args[arg_idx].find("Color")!=-1:
                hasColor = True
                if args[arg_idx].find("!Color")==-1:
                    args[arg_idx] = args[arg_idx].replace("Color", "!Color")
            else:
                kwargs["Color"] = False
    args = tuple(args)
    if not hasColor:
        kwargs["Color"] = False
    try:
        args, kwargs = JPyInterface.functions.ConvertSpecKwargsToArgs(["JobName", "TargetFile"], *args, **kwargs)
    except AttributeError:
        try:
            args, kwargs = JPyInterface.functions.ConvertSpecKwargsToArgs(["JobName"], *args, **kwargs)
        except AttributeError:
            raise AttributeError
    originalFunction, args = JPyInterface.functions.ProcessParameters(3, *args, **kwargs)
    return originalFunction(*args)

## Rewrite TMVA::Factory::BookMethod
# @param *args positional parameters
# @param **kwargs named parameters: this will be transformed to option string
def ChangeCallOriginalBookMethod(*args,  **kwargs):
    compositeOpts = False
    composite = False
    if "Composite" in kwargs:
        composite = kwargs["Composite"]
        del kwargs["Composite"]
        if "CompositeOptions" in kwargs:
            compositeOpts = kwargs["CompositeOptions"]
            del kwargs["CompositeOptions"]
    args, kwargs = JPyInterface.functions.ConvertSpecKwargsToArgs(["DataLoader", "Method", "MethodTitle"], *args, **kwargs)
    originalFunction, args = JPyInterface.functions.ProcessParameters(4, *args, **kwargs)
    if composite!=False:
        args = list(args)
        args.append(composite)
        args = tuple(args)
    if compositeOpts!=False:
        o, compositeOptStr = JPyInterface.functions.ProcessParameters(-10, **compositeOpts)
        args = list(args)
        args.append(compositeOptStr[0])
        args = tuple(args)
    return originalFunction(*args)

## Rewrite the constructor of TMVA::Factory::EvaluateImportance
# @param *args positional parameters
# @param **kwargs named parameters: this will be transformed to option string
def ChangeCallOriginalEvaluateImportance(*args,  **kwargs):
    if len(kwargs) == 0:
        originalFunction, args = JPyInterface.functions.ProcessParameters(0, *args, **kwargs)
        return originalFunction(*args)
    args, kwargs = JPyInterface.functions.ConvertSpecKwargsToArgs(["DataLoader", "VIType", "Method", "MethodTitle"], *args, **kwargs)
    originalFunction, args = JPyInterface.functions.ProcessParameters(5, *args, **kwargs)
    hist = originalFunction(*args)
    JPyInterface.JsDraw.Draw(hist)
    return hist

## Background booking method for BookDNN
__BookDNNHelper = None

## Graphical interface for booking DNN
# @param self object pointer
# @param loader the DataLoader object
# @param title classifier title
def BookDNN(self, loader, title="DNN"):
    global __BookDNNHelper
    def __bookDNN(optString):
        self.BookMethod(loader, ROOT.TMVA.Types.kDNN, title, optString)
        return
    __BookDNNHelper = __bookDNN
    clear_output()
    JPyInterface.JsDraw.InsertCSS("NetworkDesigner.min.css")
    JPyInterface.JsDraw.Draw("", "NetworkDesigner", True)

## This function gets the classifier information and weights in JSON formats, and the selected layers and it will create
# the weight heat map.
# @param net DNN in JSON format
# @param selectedLayers the selected layers
def CreateWeightHist(net, selectedLayers):
    firstLayer = int(selectedLayers.split("->")[0])
    weights = net["layers"][firstLayer]["Weights"]
    n1 = int(weights["row"])
    n2 = int(weights["cols"])
    m = ROOT.TMatrixD(n1, n2+1)
    vec = weights["data"]
    for i in xrange(n1):
        for j in xrange(n2):
            m[i][j] = vec[j+n2*i]
    bvec = net["layers"][firstLayer]["Biases"]["data"]
    if n1!=len(bvec):
        print("Something wrong.. Number of bias weights not equal with the neuron number ("+str(n1)+"!="+str(len(bvec))+")")
        return
    for i in xrange(n1):
        m[i][n2] = bvec[i]
    th2 = ROOT.TH2D(m)
    th2.SetTitle("Weight map for DNN")
    for i in xrange(n2):
        th2.GetXaxis().SetBinLabel(i + 1, str(i))
    th2.GetXaxis().SetBinLabel(n2+1, "B")
    for i in xrange(n1):
        th2.GetYaxis().SetBinLabel(i + 1, str(i))
    th2.GetXaxis().SetTitle("Layer: "+str(firstLayer))
    th2.GetYaxis().SetTitle("Layer: "+str(firstLayer+1))
    th2.SetStats(0)
    th2.SetMarkerSize(1.5)
    th2.SetMarkerColor(0)
    labelSize = 0.040
    th2.GetXaxis().SetLabelSize(labelSize)
    th2.GetYaxis().SetLabelSize(labelSize)
    th2.LabelsOption("d")
    th2.SetLabelOffset(0.011)
    clear_output()
    JPyInterface.JsDraw.Draw(th2, 'drawDNNMap')

## Show DNN weights in a heat map. It will produce an ipywidget element, where the layers can be selected.
# @param fac object pointer
# @oaram datasetName name of current dataset
# @param methodName DNN's name
def DrawDNNWeights(fac, datasetName, methodName="DNN"):
    m = GetMethodObject(fac, datasetName, methodName)
    if m == None:
        return None
    try:
        net = GetDeepNetwork(str(m.GetWeightFileName()), True)
    except AttributeError:
        print("STANDARD architecture not supported! If you want to use this function you must use CPU or GPU architecture")
    numOfLayers = len(net["layers"])
    options = []
    vals=[]
    for layer in xrange(numOfLayers):
        options.append(str(layer)+"->"+str(layer+1))
        vals.append(layer)
    selectLayer=widgets.Dropdown(
        options=options,
        value=options[0],
        description='Layer'
    )
    def drawWrapper(e):
        CreateWeightHist(net, selectLayer.value)
        pass
    button = widgets.Button(description="Draw", font_weight="bold", font_size="16")
    button.on_click(drawWrapper)
    box = widgets.HBox([selectLayer, button])
    display(box)
