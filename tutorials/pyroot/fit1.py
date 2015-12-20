## \file
## \ingroup tutorial_pyroot
## Fit example.
##
## \macro_image
## \macro_output
## \macro_code
##
## \author Wim Lavrijsen

from os import path
from ROOT import TCanvas, TFile, TPaveText
from ROOT import gROOT, gBenchmark

c1 = TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
c1.SetGridx()
c1.SetGridy()
c1.GetFrame().SetFillColor( 21 )
c1.GetFrame().SetBorderMode(-1 )
c1.GetFrame().SetBorderSize( 5 )

gBenchmark.Start( 'fit1' )
#
# We connect the ROOT file generated in a previous tutorial
#
fill = TFile( 'py-fillrandom.root' )

#
# The function "ls()" lists the directory contents of this file
#
fill.ls()

#
# Get object "sqroot" from the file.
#

sqroot = gROOT.FindObject( 'sqroot' )
sqroot.Print()

#
# Now fit histogram h1f with the function sqroot
#
h1f = gROOT.FindObject( 'h1f' )
h1f.SetFillColor( 45 )
h1f.Fit( 'sqroot' )

# We now annotate the picture by creating a PaveText object
# and displaying the list of commands in this macro
#
fitlabel = TPaveText( 0.6, 0.3, 0.9, 0.80, 'NDC' )
fitlabel.SetTextAlign( 12 )
fitlabel.SetFillColor( 42 )
fitlabel.ReadFile(path.join(path.dirname(__file__),'fit1_py.py'))
fitlabel.Draw()
c1.Update()
gBenchmark.Show( 'fit1' )
