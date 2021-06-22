import os
import unittest
from array import array

import ROOT
from DistRDF.HeadNode import Factory


class FriendInfoTest(unittest.TestCase):
    """Unit test for the FriendInfo class"""

    def create_parent_tree(self, treename, filename):
        """Creates a .root file with the parent TTree"""
        f = ROOT.TFile(filename, "recreate")
        t = ROOT.TTree(treename, "parent tree")

        x = array("i", [0])
        t.Branch("x", x, "x/I")

        for i in range(9):
            x[0] = i
            t.Fill()

        f.Write()
        f.Close()

    def create_friend_tree(self, treename, filename):
        """Creates a .root file with the friend TTree"""
        ff = ROOT.TFile(filename, "recreate")
        tf = ROOT.TTree(treename, "friend tree")

        y = array("i", [0])
        tf.Branch("y", y, "y/I")

        # The friend will have a gaussian distribution with mean 20 and
        # standard deviation 1
        for i in range(3):
            y[0] = i
            tf.Fill()

        ff.Write()
        ff.Close()

    def test_friend_info_with_ttree(self):
        """
        Check that FriendInfo correctly stores information about the friend
        trees
        """
        # Parent Tree
        base_tree_name = "T"
        base_tree_filename = "treeparent.root"
        self.create_parent_tree(base_tree_name, base_tree_filename)
        basetree = ROOT.TChain(base_tree_name)
        basetree.Add(base_tree_filename)

        # Friend Tree
        friend_tree_name = "TF"
        friend_tree_alias = "TF"
        friend_tree_filename = "treefriend.root"
        self.create_friend_tree(friend_tree_name, friend_tree_filename)
        friendtree = ROOT.TChain(friend_tree_name)
        friendtree.Add(friend_tree_filename)

        # Add friendTree to the parent
        basetree.AddFriend(friendtree, friend_tree_alias)

        # Instantiate head node of the graph with the base TTree
        headnode = Factory.get_headnode(basetree)

        # Retrieve information about friends
        friendnamesalias, friendfilenames, friendchainsubnames = headnode.get_friendinfo()

        # Check that FriendInfo has non-empty lists
        self.assertIsNotNone(friendnamesalias)
        self.assertIsNotNone(friendfilenames)
        self.assertIsNotNone(friendchainsubnames)

        # Check that the three lists with treenames, filenames and subnames are populated
        # as expected.
        self.assertTupleEqual(friendnamesalias, ((friend_tree_name, friend_tree_alias),))
        self.assertTupleEqual(friendfilenames,  ((friend_tree_filename,),))
        self.assertTupleEqual(friendchainsubnames,  ((friend_tree_name,),))

        # Remove unnecessary .root files
        os.remove(base_tree_filename)
        os.remove(friend_tree_filename)


    def test_friend_info_chain_with_subnames(self):
        """
        Check that FriendInfo correctly stores information about the friend
        trees
        """
        # Parent Tree
        parent_name = "treeparent"
        parent_filename = "treeparent.root"
        self.create_parent_tree(parent_name, parent_filename)
        parenttree = ROOT.TChain(parent_name)
        parenttree.Add(parent_filename)

        # Friend chain
        friendchainname = "treefriend"
        friendchain = ROOT.TChain(friendchainname)
        actualfriendchainfilenames = []
        actualfriendchainsubnames = []
        for i in range(1,4):
            friend_name = "treefriend" + str(i)
            friend_filename = friend_name + ".root"
            actualfriendchainfilenames.append(friend_filename)
            actualfriendchainsubnames.append(friend_name)
            self.create_friend_tree(friend_name, friend_filename)
            friendchain.Add(friend_filename + "/" + friend_name)

        # Add friendTree to the parent
        parenttree.AddFriend(friendchain)

        # Instantiate head node of the graph with the base TTree
        headnode = Factory.get_headnode(parenttree)

        # Retrieve information about friends
        friendnamesalias, friendfilenames, friendchainsubnames = headnode.get_friendinfo()

        print(friendnamesalias)
        print(friendfilenames)
        print(friendchainsubnames)
        # Check that FriendInfo has non-empty lists
        self.assertIsNotNone(friendnamesalias)
        self.assertIsNotNone(friendfilenames)
        self.assertIsNotNone(friendchainsubnames)

        # Check that the three lists with treenames, filenames and subnames are populated
        # as expected.
        self.assertTupleEqual(friendnamesalias, ((friendchainname, friendchainname),))
        self.assertTupleEqual(friendfilenames,  (tuple(actualfriendchainfilenames),) )
        self.assertTupleEqual(friendchainsubnames,  (tuple(actualfriendchainsubnames),)   )

        # Remove unnecessary .root files
        os.remove(parent_filename)
        for filename in actualfriendchainfilenames:
            os.remove(filename)
