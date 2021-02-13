## \file
## \ingroup tutorial_pyroot
## \notebook -nodraw
## example of macro to read data from an ascii file and
## create a root file with a Tree.
##
## NOTE: comparing the results of this macro with those of staff.C, you'll
## notice that the resultant file is a couple of bytes smaller, because the
## code below strips all white-spaces, whereas the .C version does not.
##
## \macro_code
##
## \author Wim Lavrijsen

import re, array, os
import ROOT
from ROOT import TFile, TTree, gROOT, addressof

## A C/C++ structure is required, to allow memory based access
gROOT.ProcessLine(
"struct staff_t {\
   Int_t           Category;\
   UInt_t          Flag;\
   Int_t           Age;\
   Int_t           Service;\
   Int_t           Children;\
   Int_t           Grade;\
   Int_t           Step;\
   Int_t           Hrweek;\
   Int_t           Cost;\
   Char_t          Division[4];\
   Char_t          Nation[3];\
};" );

## Function to read in data from ASCII file and fill the ROOT tree
def staff():

    staff = ROOT.staff_t()

    # The input file cern.dat is a copy of the CERN staff data base
    # from 1988

    f = TFile( 'staff.root', 'RECREATE' )
    tree = TTree( 'T', 'staff data from ascii file' )
    tree.Branch( 'staff', staff, 'Category/I:Flag:Age:Service:Children:Grade:Step:Hrweek:Cost' )
    tree.Branch( 'Divisions', addressof( staff, 'Division' ), 'Division/C' )
    tree.Branch( 'Nation', addressof( staff, 'Nation' ), 'Nation/C' )

    # note that the branches Division and Nation cannot be on the first branch
    fname = os.path.join(str(ROOT.gROOT.GetTutorialDir()), 'tree', 'cernstaff.dat')
    for line in open(fname).readlines():
        t = list(filter( lambda x: x, re.split( '\s+', line ) ) )
        staff.Category = int(t[0])             # assign as integers
        staff.Flag     = int(t[1])
        staff.Age      = int(t[2])
        staff.Service  = int(t[3])
        staff.Children = int(t[4])
        staff.Grade    = int(t[5])
        staff.Step     = int(t[6])
        staff.Hrweek   = int(t[7])
        staff.Cost     = int(t[8])
        staff.Division = t[9]                  # assign as strings
        staff.Nation   = t[10]

        tree.Fill()

    tree.Print()
    tree.Write()

#### run fill function if invoked on CLI
if __name__ == '__main__':
   staff()
