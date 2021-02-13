import unittest
from array import array

import ROOT
from libcppyy import SetOwnership


class TTreeIterable(unittest.TestCase):
    """
    Test for the pythonization that makes TTree instances iterable in Python. 
    For example, this allows to do:
    `for event in mytree:`
         `...`

    Since this pythonization is common to TTree and its subclasses, TChain, TNtuple
    and TNtupleD are also tested here.
    """

    filename  = 'treeiterable.root'
    treename  = 'mytree'
    tuplename = 'mytuple'
    nentries  = 10
    arraysize = 10
    more      = 10

    # Setup
    @classmethod
    def setUpClass(cls):
        ROOT.gInterpreter.Declare('#include "TreeHelper.h"')                                          
        ROOT.CreateTTree(cls.filename,
                         cls.treename,
                         cls.nentries,
                         cls.arraysize,
                         cls.more,
                         'RECREATE')
        ROOT.CreateTNtuple(cls.filename,
                           cls.tuplename,
                           cls.nentries,
                           cls.more,
                           'UPDATE')

    # Helpers
    def get_tree_and_chain(self):
        f = ROOT.TFile(self.filename)
        t = f.Get(self.treename)
        # Prevent double deletion of the tree (Python and C++ TFile)
        SetOwnership(t, False)

        c = ROOT.TChain(self.treename)
        c.Add(self.filename)
        c.Add(self.filename)

        return f,t,c

    def get_ntuples(self):
        f = ROOT.TFile(self.filename)

        nt = f.Get(self.tuplename)
        SetOwnership(nt, False)

        ntd = f.Get(self.tuplename + 'D')
        SetOwnership(ntd, False)

        return f,nt,ntd

    # Tests
    def test_basic_type_branch(self):
        f,t,c = self.get_tree_and_chain()

        for ds in t,c:
            n = array('f', [ 0. ])
            ds.SetBranchAddress('floatb', n)

            i = 0
            for entry in ds:
                self.assertEqual(n[0], i+self.more)
                i = (i + 1) % self.nentries

    def test_array_branch(self):
        f,t,c = self.get_tree_and_chain()

        for ds in t,c:
            a = array('d', self.arraysize*[ 0. ])
            ds.SetBranchAddress('arrayb', a)

            i = 0
            for entry in ds:
                for j in range(self.arraysize):
                    self.assertEqual(a[j], i+j)
                i = (i + 1) % self.nentries

    def test_struct_branch(self):
        f,t,c = self.get_tree_and_chain()

        for ds in t,c:
            ms = ROOT.MyStruct()
            ds.SetBranchAddress('structleaflistb', ms)

            i = 0
            for entry in ds:
                self.assertEqual(ms.myint1, i+self.more)
                self.assertEqual(ms.myint2, i*self.more)
                i = (i + 1) % self.nentries

    def test_ntuples(self):
        f,nt,ntd = self.get_ntuples()

        colnames = ['x','y','z']
        cols  = [ array('f', [ 0. ]) for _ in colnames ]
        colsd = [ array('d', [ 0. ]) for _ in colnames ]
        ncols = len(cols)

        for ds,cs in [(nt,cols),(ntd,colsd)]:
            for i in range(ncols):
                ds.SetBranchAddress(colnames[i], cs[i])

            numentry = 0
        
            for entry in ds:
                for i in range(ncols):
                    self.assertEqual(cs[i][0], numentry + i*self.more)
                numentry += 1


if __name__ == '__main__':
    unittest.main()
