# -*- coding: utf-8 -*-
## @package JsMVA.FormatedOutput
#  @author Attila Bagoly <battila93@gmail.com>
# This class will transform the TMVA original output to HTML formated output.

from JsMVA import DataLoader
from JsMVA.Utils import xrange
import cgi
import re

## The output transformer class. This class contains all the methods which are used to transform the C++ style output
# to HTML formatted output.
class transformTMVAOutputToHTML:

    ## Constructor
    # @param self object pointer
    def __init__(self):
        self.__eventsUID = 0

    # Transform one line of C++ output string to HTML.
    # @param self object pointer
    # @param line line to transform
    def __processGroupContentLine(self, line):
        if re.match(r"^\s*:?\s*-*\s*$", line)!=None:
            return ""
        if len(self.__outputFlagClass) > 1:
            line = cgi.escape(str(line))
        self.addClassForOutputFlag(line)
        if line.find("Booking method")!=-1:
            line = "Booking method: <b>" + line.split(":")[1].replace("\033[1m", "").replace("[0m", "")+"</b>"
        if line.find("time")!=-1:
            lr = line.split(":")
            if len(lr)==2:
                line = "<b>" + lr[0] + " : <b style='color:rgb(0,0,179)'>" +  lr[1].replace("\033[1;31m", "").replace("[0m", "") +"</b></b>"
        outFlag = self.__outputFlagClass
        if line.find("\033[0;36m")!=-1:
            line = "<b style='color:#006666'>" + line.replace("\033[0;36m", "").replace("[0m", "") + "</b>"
        if line.find("\033[1m")!=-1:
            line = "<b>" + line.replace("\033[1m", "").replace("[0m", "") + "</b>"
        lre = re.match(r"(\s*using\s*input\s*file\s*):?\s*(.*)", line, re.I)
        if lre:
            line = lre.group(1)+":"+"<b>"+lre.group(2)+"</b>"
        return "<td"+outFlag+">"+str(line)+"</td>"

    ## Checks if a line is empty.
    # @param self object pointer
    # @param line line to check
    def __isEmpty(self, line):
        return re.match(r"^\s*:?\s*-*$", line)!=None

    ## Transforms all data set specific content.
    # @param self object pointer
    # @param firstLine first line which mach with header
    # @param startIndex it defines where we are in self.lines array
    # @param maxlen defines how many iterations we can do from startIndex
    def __transformDatasetSpecificContent(self, firstLine, startIndex, maxlen):
        tmp_str = ""
        count = 0
        for l in xrange(1, maxlen):
            nextline = self.lines[startIndex+l]
            if self.__isEmpty(nextline):
                count += 1
                continue
            DatasetName = re.match(r".*(\[.*\])\s*:\s*(.*)", nextline)
            if DatasetName:
                count += 1
                tmp_str += "<tr>"+self.__processGroupContentLine(DatasetName.group(2))
                tmp_str += "<td class='tmva_output_hidden_td'></td></tr>"
            else:
                break
        DatasetName = re.match(r".*(\[.*\])\s*:\s*(.*)", firstLine)
        self.__lastDataSetName = DatasetName.group(1).replace("[", "").replace("]", "")
        tbodyclass = ""
        if count > 0:
            tbodyclass = " class='tmva_output_tbody_multiple_row'"
        rstr = "<table class='tmva_output_dataset'><tbody"+tbodyclass+">"
        rstr += "<tr><td rowspan='"+str(count+1)+"'>Dataset: "+self.__lastDataSetName+"</td>"
        rstr += self.__processGroupContentLine(DatasetName.group(2)) + "<td class='tmva_output_hidden_td'></td></tr>"
        rstr += tmp_str
        rstr += "</tbody></table>"
        return (count, rstr)

    ## It transforms number of events related contents.
    # @param self object pointer
    # @param firstLine first line which mach with header
    # @param startIndex it defines where we are in self.lines array
    # @param maxlen defines how many iterations we can do from startIndex
    def __transformNumEvtSpecificContent(self, firstLine, startIndex, maxlen):
        tmp_str = ""
        count = 0
        tmpmap = {}
        for l in xrange(1, maxlen):
            nextline = self.lines[startIndex+l]
            if self.__isEmpty(nextline):
                count += 1
                continue
            NumberOfEvents = re.match(r"\s+:\s*\w+\s*-\s*-\s*((training\sevents)|(testing\sevents)|(training\sand\stesting\sevents))\s*:\s*\d+", nextline)
            if NumberOfEvents:
                lc = re.findall(r"\w+", nextline)
                t = ""
                for i in range(1, len(lc) - 1):
                    t += lc[i] + " "
                t = t[:-1]
                if lc[0] not in tmpmap:
                    tmpmap[lc[0]] = []
                count += 1
                tmpmap[lc[0]].append({"name": t, "value": lc[len(lc) - 1]})
            else:
                break
        rstr = "<table class='tmva_output_traintestevents'>"
        rstr += "<tr><td colspan='3'><center><b>"+firstLine+"</b></center></td></tr>"
        for key in tmpmap:
            rstr += "<tr>"
            rstr += "<td rowspan='"+str(len(tmpmap[key]))+"'>"+key+"</td>"
            rstr += "<td>"+tmpmap[key][0]["name"]+"</td><td>"+tmpmap[key][0]["value"]+"</td>"
            rstr += "<td class='tmva_output_hidden_td'></td>"
            rstr += "</tr>"
            for i in xrange(1, len(tmpmap[key])):
                rstr += "<tr><td>"+tmpmap[key][i]["name"]+"</td><td>"+tmpmap[key][i]["value"]+"</td>"
                rstr += "<td class='tmva_output_hidden_td'></td></tr>"
        rstr += tmp_str
        rstr += "</table>"
        return (count, rstr)

    ## Transform Variable related informations to table.
    # @param self object pointer
    # @param headerMatch re.match object for the first line
    # @param startIndex it defines where we are in self.lines array
    # @param maxlen defines how many iterations we can do from startIndex
    def __transformVariableMeanSpecificContent(self, headerMatch, startIndex, maxlen):
        count = 0
        table = [[]]
        for j in range(1, 6):
            table[0].append(headerMatch.group(j))
        for l in xrange(1, maxlen):
            nextline = self.lines[startIndex + l]
            if self.__isEmpty(nextline):
                count += 1
                continue
            VariableMean = re.match(r"\s*:\s*([\w\d]+)\s*:\s*(-?\d*\.?\d*)\s*(-?\d*\.?\d*)\s*\[\s*(-?\d*\.?\d*)\s*(-?\d*\.?\d*)\s*\]", nextline, re.I)
            if VariableMean:
                count += 1
                tmp = []
                for j in range(1, 6):
                    tmp.append(VariableMean.group(j))
                table.append(tmp)
            else:
                break
        rstr = "<table class='tmva_output_varmeanrms'>"
        for i in xrange(len(table)):
            rstr += "<tr>"
            for j in xrange(len(table[i])):
                rstr += "<td>" + str(table[i][j])
            rstr += "<td class='tmva_output_hidden_td'></td>"
            rstr += "</tr>"
        rstr += "</table>"
        return (count, rstr)

    ## This function creates a correlation matrix from python matrix.
    # @param self object pointer
    # @param title first part of the title of histogram
    # @param className Signal/Background: it will be on title
    # @param varNames array with variable names
    # @param matrix the correlation matrix
    def __correlationMatrix(self, title, className, varNames, matrix):
        id = "jsmva_outputtansformer_events_"+str(self.__eventsUID)+"_onclick"
        self.__eventsUID += 1
        json = DataLoader.GetCorrelationMatrixInJSON(className, varNames, matrix)
        jsCall = "require(['JsMVA'],function(jsmva){jsmva.outputShowCorrelationMatrix('"+id+"');});"
        rstr = "<div id='"+id+"' style='display: none; width: 600px; height:350'>"+json+"</div>"
        rstr += self.__processGroupContentLine("<a onclick=\""+jsCall+"\" class='tmva_output_corrmat_link'>" + title + " (" + className + ")</a>")
        return rstr

    ## This function add different flag based on the line formation. It add's a specific class for each output type.
    # @param self object pointer
    # @param line current line
    def addClassForOutputFlag(self, line):
        if self.__currentType==None:
            self.__outputFlagClass = ""
            if line.find("\033[1;31m") != -1:
                self.__outputFlagClass = " class='tmva_output_warning'"
            elif line.find("\033[31m") !=-1:
                self.__outputFlagClass = " class='tmva_output_error'"
            elif line.find("\033[37;41;1m") !=-1:
                self.__outputFlagClass = " class='tmva_output_fatal'"
            elif line.find("\033[37;41;1m") !=-1:
                self.__outputFlagClass = " class='tmva_output_silent'"
            elif line.find("\033[34m") !=-1:
                self.__outputFlagClass = " class='tmva_output_debug'"
            return
        if self.__currentType.find("WARNING") != -1 or line.find("\033[1;31m")!=-1:
            self.__outputFlagClass = " class='tmva_output_warning'"
        elif self.__currentType.find("ERROR") !=-1:
            self.__outputFlagClass = " class='tmva_output_error'"
        elif self.__currentType.find("FATAL") !=-1:
            self.__outputFlagClass = " class='tmva_output_fatal'"
        elif self.__currentType.find("SILENT") !=-1:
            self.__outputFlagClass = " class='tmva_output_silent'"
        elif self.__currentType.find("DEBUG") !=-1:
            self.__outputFlagClass = " class='tmva_output_debug'"
        elif self.__currentType=="impHead":
            self.__outputFlagClass = " class='tmva_output_imphead'"
        else:
            self.__outputFlagClass = ""

    ## This function can transform one group of output. The group has only one header.
    # @param self object pointer
    # @param firstLine the header line, it defines the start
    def __transformOneGroup(self, firstLine):
        tmp_str = ""
        processed_lines = 0
        lineIter = iter(xrange(len(self.lines) - self.lineIndex))
        for j in lineIter:
            if j == 0:
                nextline = firstLine
            else:
                nextline = self.lines[self.lineIndex + j]
            Header = re.match(r"^\s*-*\s*(<\w+>\s*)*\s*(\w+.*\s+)(\s+)(:)\s*(.*)", nextline, re.I)
            EmptyHeader = re.match(r"\s*-*\s*(<\w+>\s*)+\s*:\s*(.*)", nextline, re.I)
            DatasetName = re.match(r".*(\[.*\])\s*:\s*(.*)", nextline)
            NumEvents = re.match(r"(.*)(number\sof\straining\sand\stesting\sevents)", nextline, re.I)
            CorrelationMatrixHeader = re.match(r".*\s*:?\s*(correlation\s*matrix)\s*\((\w+)\)\s*:\s*", nextline, re.I)
            VariableMeanHeader = re.match(r".*\s*:?\s*(variable)\s*(mean)\s*(rms)\s*\[\s*(min)\s*(max)\s*\].*", nextline, re.I)
            WelcomeHeader = re.match(r"^\s*:?\s*(.*ROOT\s*version:.*)", nextline, re.I)
            if (Header == None and EmptyHeader ==None) or j == 0:
                if j != 0:
                    processed_lines += 1
                    self.iterLines.next()
                    tmp_str += "<tr>"
                if DatasetName or NumEvents or VariableMeanHeader:
                    if DatasetName:
                        func = self.__transformDatasetSpecificContent
                        fLine = nextline
                    elif NumEvents:
                        func = self.__transformNumEvtSpecificContent
                        fLine = NumEvents.group(2)
                    else:
                        func = self.__transformVariableMeanSpecificContent
                        fLine = VariableMeanHeader
                    count, tmp = func(fLine, self.lineIndex + j, len(self.lines) - self.lineIndex - j)
                    for x in xrange(count):
                        lineIter.next()
                        self.iterLines.next()
                    tmp_str += self.__processGroupContentLine(tmp)
                elif CorrelationMatrixHeader:
                    self.iterLines.next()
                    lineIter.next()
                    ik = 1
                    matrixLines = []
                    while True:
                        self.iterLines.next()
                        lineIter.next()
                        ik += 1
                        if self.__isEmpty(self.lines[self.lineIndex + j + ik]):
                            break
                        ltmp = re.match(r"\s*:\s*(.*)", self.lines[self.lineIndex + j + ik]).group(1)
                        if ltmp.find(":") != -1:
                            matrixLines.append(ltmp.split(":"))
                    rmatch = "^\s*"
                    for ii in xrange(len(matrixLines)):
                        rmatch += r"(\+\d+\.?\d*|-\d+\.?\d*)\s*"
                    rmatch += "$"
                    matrix = []
                    varNames = []
                    for ii in xrange(len(matrixLines)):
                        varNames.append(matrixLines[ii][0])
                        ll = re.match(rmatch, matrixLines[ii][1])
                        mline = []
                        for jj in xrange(len(matrixLines)):
                            mline.append(float(ll.group(jj + 1)))
                        matrix.append(mline)
                    tmp_str += self.__correlationMatrix(CorrelationMatrixHeader.group(1),
                                                        CorrelationMatrixHeader.group(2), varNames, matrix)
                elif WelcomeHeader:
                    kw = 0
                    while True:
                        nextline = self.lines[self.lineIndex + j + kw]
                        EndWelcome = re.match(r"\s*:?\s*_*\s*(TMVA\s*Version\s*.*)", nextline, re.I)
                        if EndWelcome or re.match(r"[\s_/|]*", nextline) == None:
                            break
                        kw += 1
                        self.iterLines.next()
                        lineIter.next()
                    tmp_str += "<td><b>" + WelcomeHeader.group(1) + "</b></td></tr>"
                    tmp_str += "<tr><td><img src='https://rawgit.com/root-project/root/master/tutorials/tmva/tmva_logo.svg' width='100%' /><br />"
                    tmp_str += "<center><b>" + EndWelcome.group(1) + "</b></center></td></tr>"
                else:
                    lmatch = re.match(r"\s*(<\w+>\s*)*\s*:\s*(.*)", nextline)
                    if lmatch:
                        if lmatch.group(1):
                            self.__currentType = lmatch.group(1)
                        tmp_str += self.__processGroupContentLine(lmatch.group(2))
                    else:
                        tmp_str += self.__processGroupContentLine(nextline)
                tmp_str += "</tr>"
            else:
                break
        tbodyclass = ""
        if processed_lines > 0:
            tbodyclass = " class='tmva_output_tbody_multiple_row'"
        extraHeaderClass=""
        if self.__currentType=="impHead":
            extraHeaderClass = " tmva_output_imphead"
        self.out += "<tbody"+tbodyclass+"><tr><td rowspan='" + str(processed_lines + 1) + "' class='tmva_output_header"+extraHeaderClass+"'>" + self.__currentHeaderName + "</td>"
        self.out += tmp_str + "</tbody>"

    ## This is the main function, it will be registered as transformer function to JupyROOT, it will run every time
    # when some ROOT function produces output. It get's the C++ style output and it returns it in HTML version.
    # @param self object pointer
    # @param output the content of C++ output stream
    # @param error the content of C++ error stream
    def transform(self, output, error):
        self.err = ""
        if str(error).find(", time left:")==-1:
            self.err = error
        self.out = ""
        self.lines = output.splitlines()
        self.iterLines = iter(xrange(len(self.lines)))
        for self.lineIndex in self.iterLines:
            line = self.lines[self.lineIndex]
            Header = re.match(r"^\s*-*\s*(<\w+>\s*)*\s*(\w+.*\s+)(\s+)(:)\s*(.*)", line, re.I)
            EmptyHeader = re.match(r"\s*-*\s*(<\w+>\s*)+\s*:\s*(.*)", line, re.I)
            NewGroup = re.match(r"^\s*:\s*$", line)

            if Header:
                if re.match(r"^\s*-+\s*(<\w+>\s*)*\s*(\w+.*\s+)(\s+)(:)\s*(.*)", line, re.I):
                    self.__currentType = "impHead"
                else:
                    self.__currentType = Header.group(1)
                self.addClassForOutputFlag(Header.group(5))
                self.__currentHeaderName = Header.group(2)
                self.__transformOneGroup(Header.group(5))
            elif EmptyHeader:
                self.__currentType = EmptyHeader.group(1)
                self.addClassForOutputFlag(EmptyHeader.group(2))
                self.__currentHeaderName = ""
                self.__transformOneGroup(EmptyHeader.group(2))
            elif NewGroup:
                kw = 1
                lines = []
                while True:
                    if (self.lineIndex + kw) >=len(self.lines):
                        break
                    nextline = self.lines[self.lineIndex + kw]
                    kw += 1
                    if re.match(r"^\s*:\s*$", nextline):
                        Header = re.match(r"^\s*-*\s*(<\w+>\s*)*\s*(\w+.*\s+)(\s+)(:)\s*(.*)", self.lines[self.lineIndex+kw], re.I)
                        if Header:
                            self.iterLines.next()
                        break
                    self.iterLines.next()
                    if self.__isEmpty(nextline):
                        continue
                    lre = re.match(r"\s*:\s*(.*)", nextline)
                    if lre:
                        lines.append(lre.group(1))
                    else:
                        lines.append(nextline)
                if len(lines)==0:
                    continue
                if len(lines) > 1:
                    tbodyclass = " class='tmva_output_tbody_multiple_row'"
                self.out += "<tbody"+tbodyclass+"><tr><td rowspan='"+str(len(lines))+"'><td>"+lines[0]+"</td></tr>"
                for ii in xrange(1, len(lines)):
                    self.out += "<tr><td>"+lines[ii]+"</td></tr>"
            else:
                lre = re.match(r"\s*:\s*(.*)", line)
                if lre:
                    self.out += "<tr><td>" + lre.group(1) + "</td></tr>"
                else:
                    self.out += "<tr><td>" + line + "</td></tr>"
        if len(self.out) < 1 and len(self.err) < 1:
            return ("", "", "")
        self.out = "<table class='tmva_output_table'>" + self.out + "</table>"
        return (self.out, self.err, "html")
