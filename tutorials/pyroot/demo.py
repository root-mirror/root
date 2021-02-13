## \file
## \ingroup tutorial_pyroot
## To run, do "python <path-to>/demo.py"
##
## \macro_code
##
## \author Wim Lavrijsen, Enric Tejedor

import os, sys
import ROOT

# To run, do "python <path-to>/demo.py"

# enable running from another directory than the one where demo.py resides
workdir = os.path.dirname( sys.argv[0] )
if workdir:
   os.chdir( workdir )

# This macro generates a Controlbar menu.
# To execute an item, click with the left mouse button.
# To see the HELP of a button, click on the right mouse button.

ROOT.gStyle.SetScreenFactor(1)   # if you have a large screen, select 1.2 or 1.4

bar = ROOT.TControlBar( 'vertical', 'Demos', 10, 10 )

# The callbacks to python work by having CLING call the python interpreter through
# the "TPython" class. Note the use of "raw strings."
if sys.version_info[0] > 2:
    to_run = 'exec(open(\'{}\').read())'
else:
    to_run = 'execfile(\'{}\')'

bar.AddButton( 'Help on Demos', r'TPython::Exec( "' + to_run.format('demoshelp.py') + '" );', 'Click Here For Help on Running the Demos' )
bar.AddButton( 'browser',       r'TPython::Exec( "b = ROOT.TBrowser()" );',          'Start the ROOT browser' )
bar.AddButton( 'framework',     r'TPython::Exec( "' + to_run.format('../legacy/pyroot/framework.py') + '" );', 'An Example of Object Oriented User Interface' )
bar.AddButton( 'first',         r'TPython::Exec( "' + to_run.format('../legacy/pyroot/first.py') + '" );',     'An Example of Slide with Root' )
bar.AddButton( 'hsimple',       r'TPython::Exec( "' + to_run.format('hsimple.py') + '" );',   'Creating histograms/Ntuples on file', "button" )
bar.AddButton( 'hsum',          r'TPython::Exec( "' + to_run.format('hsum.py') + '" );',      'Filling Histograms and Some Graphics Options' )
bar.AddButton( 'formula1',      r'TPython::Exec( "' + to_run.format('formula1.py') + '" );',  'Simple Formula and Functions' )
bar.AddButton( 'surfaces',      r'TPython::Exec( "' + to_run.format('surfaces.py') + '" );',  'Surface Drawing Options' )
bar.AddButton( 'fillrandom',    r'TPython::Exec( "' + to_run.format('fillrandom.py') + '" );','Histograms with Random Numbers from a Function' )
bar.AddButton( 'fit1',          r'TPython::Exec( "' + to_run.format('fit1.py') + '" );',      'A Simple Fitting Example' )
bar.AddButton( 'multifit',      r'TPython::Exec( "' + to_run.format('multifit.py') + '" );',  'Fitting in Subranges of Histograms' )
bar.AddButton( 'h1draw',        r'TPython::Exec( "' + to_run.format('h1draw.py') + '" );',    'Drawing Options for 1D Histograms' )
bar.AddButton( 'graph',         r'TPython::Exec( "' + to_run.format('graph.py') + '" );',     'Example of a Simple Graph' )
bar.AddButton( 'gerrors',       r'TPython::Exec( "' + to_run.format('gerrors.py') + '" );',   'Example of a Graph with Error Bars' )
bar.AddButton( 'tornado',       r'TPython::Exec( "' + to_run.format('tornado.py') + '" );',   'Examples of 3-D PolyMarkers' )
bar.AddButton( 'shapes',        r'TPython::Exec( "' + to_run.format('shapes.py') + '" );',    'The Geometry Shapes' )
bar.AddButton( 'geometry',      r'TPython::Exec( "' + to_run.format('geometry.py') + '" );',  'Creation of the NA49 Geometry File' )
bar.AddButton( 'na49view',      r'TPython::Exec( "' + to_run.format('na49view.py') + '" );',  'Two Views of the NA49 Detector Geometry' )
bar.AddButton( 'file',          r'TPython::Exec( "' + to_run.format('../legacy/pyroot/file.py') + '" );',      'The ROOT File Format' )
bar.AddButton( 'fildir',        r'TPython::Exec( "' + to_run.format('../legacy/pyroot/fildir.py') + '" );',    'The ROOT File, Directories and Keys' )
bar.AddButton( 'tree',          r'TPython::Exec( "' + to_run.format('../legacy/pyroot/tree.py') + '" );',      'The Tree Data Structure' )
bar.AddButton( 'ntuple1',       r'TPython::Exec( "' + to_run.format('ntuple1.py') + '" );',   'Ntuples and Selections' )
bar.AddButton( 'rootmarks',     r'TPython::Exec( "' + to_run.format('../legacy/pyroot/rootmarks.py') +'" );', 'Prints an Estimated ROOTMARKS for Your Machine' )
bar.AddSeparator()       # not implemented
bar.AddButton( 'make ntuple',   r'TPython::Exec( "' + to_run.format('mrt.py') + '" );',       'Convert a text file to an ntuple' )

bar.Show()

ROOT.gROOT.SaveContext()


## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
      # Check if we are in Python 2 or 3
      if sys.version_info[0] > 2:
         rep = input( 'enter "q" to quit: ' )
      else:
         rep = raw_input( 'enter "q" to quit: ' )
      if 1 < len(rep):
         rep = rep[0]
