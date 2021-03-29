# @author Vincenzo Eduardo Padulano
#  @author Enric Tejedor
#  @date 2021-02

################################################################################
# Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from __future__ import print_function

import collections
import glob
import logging
import warnings

import ROOT

from DistRDF.Operation import Operation

logger = logging.getLogger(__name__)

TreeInfo = collections.namedtuple("TreeInfo", ["treenames", "treefilenames", "friendnamesalias", "friendfilenames"])
Range = collections.namedtuple("Range", ["start", "end", "treename",
                                         "treefilenames", "friendnamesalias", "friendfilenames"])


def _n_even_chunks(iterable, n_chunks):
    """
    Yield `n_chunks` as even chunks as possible from `iterable`. Though generic,
    this function is used in _get_clustered_ranges to split a list of clusters
    into multiple sublists. Each sublist will hold the clusters that should fit
    in a single partition of the distributed dataset::

        [
            # Partition 1 will process the following clusters
            [
                (start_0_0, end_0_0, offset_0, (filename_0, 0)),
                (start_0_1, end_0_1, offset_0, (filename_0, 0)),
                ...,
                (start_1_0, end_1_0, offset_1, (filename_1, 1)),
                (start_1_1, end_1_1, offset_1, (filename_1, 1)),
                ...,
                (start_n_0, end_n_0, offset_n, (filename_n, n)),
                (start_n_1, end_n_1, offset_n, (filename_n, n)),
                ...
            ],
            # Partition 2 will process these other clusters
            [
                (start_n+1_0, end_n+1_0, offset_n+1, (filename_n+1, n+1)),
                (start_n+1_1, end_n+1_1, offset_n+1, (filename_n+1, n+1)),
                ...,
                (start_m_0, end_m_0, offset_m, (filename_m, m)),
                (start_m_1, end_m_1, offset_m, (filename_m, m)),
                ...
            ],
            ...
        ]

    """
    last = 0
    itlenght = len(iterable)
    for i in range(1, n_chunks + 1):
        cur = int(round(i * (itlenght / n_chunks)))
        yield iterable[last:cur]
        last = cur


class Node(object):
    """
    A Class that represents a node in RDataFrame operations graph. A Node
    houses an operation and has references to children nodes.
    For details on the types of operations supported, try :

    Example::

        import DistRDF
        DistRDF.use(...) # Choose your backend
        print(DistRDF.current_backend.supported_operations)

    Attributes:
        get_head (function): A lambda function that returns the head node of
            the current graph.

        operation: The operation that this Node represents.
            This could be :obj:`None`.

        children (list): A list of :obj:`DistRDF.Node` objects which represent
            the children nodes connected to the current node.

        _new_op_name (str): The name of the new incoming operation of the next
            child, which is the last child node among the current node's
            children.

        value: The computed value after executing the operation in the current
            node for a particular DistRDF graph. This is permanently :obj:`None`
            for transformation nodes and the action nodes get a
            :obj:`ROOT.RResultPtr` after event-loop execution.

        pyroot_node: Reference to the PyROOT object that implements the
            functionality of this node on the cpp side.

        has_user_references (bool): A flag to check whether the node has
            direct user references, that is if it is assigned to a variable.
            Default value is :obj:`True`, turns to :obj:`False` if the proxy
            that wraps the node gets garbage collected by Python.
    """

    def __init__(self, get_head, operation, *args):
        """
        Creates a new node based on the operation passed as argument.

        Args:
            get_head (function): A lambda function that returns the head node
                of the current graph. This value could be `None`.

            operation (DistRDF.Operation.Operation): The operation that this Node
                represents. This could be :obj:`None`.
        """
        if get_head is None:
            # Function to get 'head' Node
            self.get_head = lambda: self
        else:
            self.get_head = get_head

        self.operation = operation
        self.children = []
        self._new_op_name = ""
        self.value = None
        self.pyroot_node = None
        self.has_user_references = True

    def __getstate__(self):
        """
        Converts the state of the current node
        to a Python dictionary.

        Returns:
            dictionary: A dictionary that stores all instance variables
            that represent the current DistRDF node.

        """
        state_dict = {'children': self.children}
        if self.operation:
            state_dict['operation_name'] = self.operation.name
            state_dict['operation_args'] = self.operation.args
            state_dict['operation_kwargs'] = self.operation.kwargs

        return state_dict

    def __setstate__(self, state):
        """
        Retrieves the state dictionary of the current
        node and sets the instance variables.

        Args:
            state (dict): This is the state dictionary that needs to be
                converted to a `Node` object.

        """
        self.children = state['children']
        if state.get('operation_name'):
            self.operation = Operation(state['operation_name'],
                                       *state['operation_args'],
                                       **state["operation_kwargs"])
        else:
            self.operation = None

    def is_prunable(self):
        """
        Checks whether the current node can be pruned from the computational
        graph.

        Returns:
            bool: True if the node has no children and no user references or
            its value has already been computed, False otherwise.
        """
        if not self.children:
            # Every pruning condition is written on a separate line
            if not self.has_user_references or \
               (self.operation and self.operation.is_action() and self.value):

                # ***** Condition 1 *****
                # If the node is wrapped by a proxy which is not directly
                # assigned to a variable, then it will be flagged for pruning

                # ***** Condition 2 *****
                # If the current node's value was already
                # computed, it should get pruned only if it's
                # an Action node.

                # Logger debug statements
                logger.debug("{} node can be pruned".format(
                    self.operation.name
                ))

                return True

        # Logger debug statements
        if self.operation:  # Node has an operation
            logger.debug("{} node shouldn't be pruned".format(
                self.operation.name
            ))
        else:  # Node is the RDataFrame
            logger.debug("Graph pruning completed")
        return False

    def graph_prune(self):
        """
        Prunes nodes from the current DistRDF graph under certain conditions.
        The current node will be pruned if it has no children and the user
        application does not hold any reference to it. The children of the
        current node will get recursively pruned.

        Returns:
            bool: True if the current node has to be pruned, False otherwise.
        """
        children = []

        # Logger debug statements
        if self.operation:
            logger.debug("Checking {} node for pruning".format(
                self.operation.name
            ))
        else:
            logger.debug("Starting computational graph pruning")

        for n in self.children:
            # Logger debug statement
            # Select children based on pruning condition
            if not n.graph_prune():
                children.append(n)

        self.children = children
        return self.is_prunable()


class RangesBuilder(object):

    def __init__(self, headnode):
        """Initialization"""

        self._headnode = headnode

    @property
    def npartitions(self):
        """Retrieve the npartitions value from the head node of the graph."""
        return self._headnode.npartitions

    @npartitions.setter
    def npartitions(self, value):
        """Change the npartitions value of the head node."""
        self._headnode.npartitions = value

    def get_clusters(self, treename, filelist):
        """
        Extract a list of cluster boundaries for the given tree and files

        Args:
            treename (str): Name of the TTree split into one or more files.

            filelist (list): List of one or more ROOT files.

        Returns:
            list: List of tuples defining the cluster boundaries. Each tuple
            contains four elements: first entry of a cluster, last entry of
            cluster (exclusive), offset of the cluster and file where the
            cluster belongs to::

                [
                    (0, 100, 0, ("filename_1.root", 0)),
                    (100, 200, 0, ("filename_1.root", 0)),
                    ...,
                    (10000, 10100, 10000, ("filename_2.root", 1)),
                    (10100, 10200, 10000, ("filename_2.root", 1)),
                    ...,
                    (n, n+100, n, ("filename_n.root", n)),
                    (n+100, n+200, n, ("filename_n.root", n)),
                    ...
                ]
        """

        clusters = []
        cluster = collections.namedtuple(
            "cluster", ["start", "end", "offset", "filetuple"])
        fileandindex = collections.namedtuple("fileandindex",
                                              ["filename", "index"])
        offset = 0
        fileindex = 0

        for filename in filelist:
            f = ROOT.TFile.Open(filename)
            t = f.Get(treename)

            entries = t.GetEntriesFast()
            it = t.GetClusterIterator(0)
            start = it()
            end = 0

            while start < entries:
                end = it()
                clusters.append(cluster(start + offset, end + offset, offset,
                                        fileandindex(filename, fileindex)))
                start = end

            fileindex += 1
            offset += entries

        logger.debug("Returning files with their clusters:\n%s",
                     "\n\n".join(map(str, clusters)))

        return clusters

    def _get_balanced_ranges(self, nentries):
        """
        Builds range pairs from the given values of the number of entries in
        the dataset and number of partitions required. Each range contains the
        same amount of entries, except for those cases where the number of
        entries is not a multiple of the partitions.

        Args:
            nentries (int): The number of entries in a dataset.

        Returns:
            list: List of :obj:`Range`s objects.
        """
        partition_size = int(nentries / self.npartitions)

        i = 0  # Iterator

        ranges = []

        remainder = nentries % self.npartitions

        while i < nentries:
            # Start value of current range
            start = i
            end = i = start + partition_size

            if remainder:
                # If the modulo value is not
                # exhausted, add '1' to the end
                # of the current range
                end = i = end + 1
                remainder -= 1

            ranges.append(Range(start, end, None, None, None, None))

        return ranges

    def _get_clustered_ranges(self, treeinfo):
        """
        Builds ``Range`` objects taking into account the clusters of the
        dataset. Each range will represent the entries processed within a single
        partition of the distributed dataset.

        Args:
            treeinfo (TreeInfo): namedtuple holding information about the tree

        Returns:
            list[collections.namedtuple]: List containinig the ranges in which
            the dataset has been split for distributed execution. Each ``Range``
            contains a starting entry, an ending entry, the list of files
            that are traversed to get all the entries and information about
            friend trees::

                [
                    Range(start=0,
                        end=42287856,
                        treename="Events",
                        treefilenames=['Run2012B_TauPlusX.root',
                                  'Run2012C_TauPlusX.root'],
                        friendnamesalias=None,
                        friendfilenames=None),
                    Range(start=6640348,
                        end=51303171,
                        treename="Events",
                        treefilenames=['Run2012C_TauPlusX.root'],
                        friendnamesalias=None,
                        friendfilenames=None)
                ]

        """

        # Retrieve a list of clusters for all files of the tree
        clustersinfiles = self.get_clusters(treeinfo.treenames[0], treeinfo.treefilenames)
        numclusters = len(clustersinfiles)

        # Restrict `npartitions` if it's greater than clusters of the dataset
        if self.npartitions > numclusters:
            msg = ("Number of partitions is greater than number of clusters "
                   "in the dataset. Using {} partition(s)".format(numclusters))
            warnings.warn(msg, UserWarning, stacklevel=2)
            self.npartitions = numclusters

        logger.debug("%s clusters will be split along %s partitions.",
                     numclusters, self.npartitions)

        """
        This list comprehension builds ``Range`` tuples with the following
        elements:

        1. ``start``: The minimum entry among all the clusters considered in a
           given partition. The offset of the first cluster of the list is
           subtracted. This is useful to keep the reference of the range with
           respect to the current files (see below).
        2. ``end``: The maximum entry among all the clusters considered in a
           given partition. The offset of the first cluster of the list is
           subtracted. This is useful to keep the reference of the range with
           respect to the current files (see below).
        3. ``treename``: The name of the current tree.
        3. ``treefilenames``: The list of files that are span between entries
           ``start`` and ``end``::

                Filelist: [file_1,file_2,file_3,file_4]

                Clustered range: [0,150]
                file_1 holds entries [0, 100]
                file_2 holds entries [101, 200]
                Then the clustered range should open [file_1, file_2]

                Clustered range: [150,350]
                file_3 holds entris [201, 300]
                file_4 holds entries [301, 400]
                Then the clustered range should open [file_2, file_3, file_4]

           Each ``cluster`` namedtuple has a ``fileandindex`` namedtuple. The
           second element of this tuple corresponds to the index of the file in
           the input `TChain`. This way all files can be uniquely identified,
           even if there is some repetition (e.g. when building a TChain with
           multiple instances of the same file). The algorithm to retrieve the
           correct files for each range takes the unique filenames from the list
           of clusters and sorts them by their index to keep the original order.

           In each file only the clusters needed to process the clustered range
           will be read.
        4. ``friendnamesalias`` and ``friendfilenames``: Information about
           friend trees.

        In each range, the offset of the first file is always subtracted to the
        ``start`` and ``end`` entries. This is needed to maintain a reference of
        the entries of the range with respect to the list of files that hold
        them. For example, given the following files::

            tree10000entries10clusters.root --> 10000 entries, 10 clusters
            tree20000entries10clusters.root --> 20000 entries, 10 clusters
            tree30000entries10clusters.root --> 30000 entries, 10 clusters

        Building 2 ranges will lead to the following tuples::

            Range(start=0,
                  end=20000,
                  treename="mytree",
                  treefilenames=['tree10000entries10clusters.root',
                                 'tree20000entries10clusters.root'],
                  friendnames=None,
                  friendfilenames=None)

            Range(start=10000,
                  end=50000,
                  treename="mytree",
                  treefilenames=['tree20000entries10clusters.root',
                                 'tree30000entries10clusters.root'],
                  friendnames=None,
                  friendfilenames=None)

        The first ``Range`` will read the first 10000 entries from the first
        file, then switch to the second file and read the first 10000 entries.
        The second ``Range`` will start from entry number 10000 of the second
        file up until the end of that file (entry number 20000), then switch to
        the third file and read the whole 30000 entries there.
        """
        clustered_ranges = [
            Range(
                min(clusters)[0] - clusters[0].offset,  # type: int
                max(clusters)[1] - clusters[0].offset,  # type: int
                treeinfo.treenames[0],  # type: str
                [
                    filetuple.filename
                    for filetuple in sorted(set([
                        cluster.filetuple for cluster in clusters
                    ]), key=lambda curtuple: curtuple[1])
                ],  # type: list[str]
                treeinfo.friendnamesalias,  # type: list[str]
                treeinfo.friendfilenames  # type: list[str]
            )  # type: collections.namedtuple
            for clusters in _n_even_chunks(clustersinfiles, self.npartitions)
        ]

        logger.debug("Created following clustered ranges:\n%s",
                     "\n\n".join(map(str, clustered_ranges)))

        return clustered_ranges

    def build_ranges(self):
        """
        Define two type of ranges based on the arguments passed to the
        RDataFrame head node.
        """
        # Create variables here to call properties only once
        nentries = self._headnode.get_num_entries()
        treeinfo = self._headnode.get_treeinfo()

        # Empty trees cannot be processed distributedly
        if not nentries:
            raise RuntimeError(
                ("No entries in the TTree. "
                 "Distributed computation will fail. "
                 "Please make sure your dataset is not empty."))

        if self.npartitions > nentries:
            # Restrict 'npartitions' if it's greater than 'nentries'
            self.npartitions = nentries

        if treeinfo is not None:

            if treeinfo.treenames and treeinfo.treefilenames:
                logger.debug("Building clustered ranges for tree %s with the "
                             "following input files:\n%s",
                             treeinfo.treenames[0],
                             treeinfo.treefilenames
                             )
                return self._get_clustered_ranges(treeinfo)
            elif not treeinfo.treefilenames:
                raise RuntimeError(("In-memory-only TTrees are not supported. "
                                    "Please make sure your tree is saved to a file first."))
        else:
            logger.debug(
                "Building balanced ranges for %d entries.", nentries)
            return self._get_balanced_ranges(nentries)


class HeadNode(Node):
    """
    The Python equivalent of ROOT C++'s
    RDataFrame class.

    Attributes:
        args (list): A list of arguments that were provided to construct
            the RDataFrame object.


    DistRDF's RDataFrame constructor accepts the same arguments as the ROOT's
    RDataFrame constructor (see
    `RDataFrame <https://root.cern/doc/master/classROOT_1_1RDataFrame.html>`_)
    """

    def __init__(self, *args):
        """
        Creates a new RDataFrame instance for the given arguments.

        Args:
            *args (list): Variable length argument list to construct the
                RDataFrame object.
        """
        super(HeadNode, self).__init__(None, None, *args)
        # Early check of the arguments to RDataFrame constructor
        df_check_args = ROOT.RDataFrame(*args)
        self.args = list(args)  # Make args mutable

        # Set at creation of the dataframe, might be optimized by the backend
        # in optimize_partitions
        self.npartitions = 2

    def build_ranges(self):
        """Build the ranges for this dataset."""
        return RangesBuilder(self).build_ranges()

    def get_branches(self):
        """Gets list of default branches if passed by the user."""
        # ROOT Constructor:
        # RDataFrame(TTree& tree, defaultBranches = {})
        if len(self.args) == 2 and isinstance(self.args[0], ROOT.TTree):
            return self.args[1]
        # ROOT Constructors:
        # RDataFrame(treeName, filenameglob, defaultBranches = {})
        # RDataFrame(treename, filenames, defaultBranches = {})
        # RDataFrame(treeName, dirPtr, defaultBranches = {})
        if len(self.args) == 3:
            return self.args[2]

        return None

    def get_num_entries(self):
        """
        Gets the number of entries in the given dataset.

        Returns:
            int: This is the computed number of entries in the input dataset.

        """
        first_arg = self.args[0]
        if isinstance(first_arg, int):
            # If there's only one argument
            # which is an integer, return it.
            return first_arg
        elif isinstance(first_arg, ROOT.TTree):
            # If the argument is a TTree or TChain,
            # get the number of entries from it.
            return first_arg.GetEntries()

        second_arg = self.args[1]

        # Construct a ROOT.TChain object
        chain = ROOT.TChain(first_arg)

        if isinstance(second_arg, str):
            # If the second argument is a string
            chain.Add(second_arg)
        else:
            # If the second argument is a list or vector
            for fname in second_arg:
                chain.Add(str(fname))  # Possible bug in conversion of string

        return chain.GetEntries()

    def get_treeinputfiles(self, tree):
        """
        Get list of input files for a TTree or TChain.

        Returns:
            (list, None): list of input files for the TTree used to construct
            this dataframe (may contain globbing characters), None otherwise.
        """
        if isinstance(tree, ROOT.TChain):
            return [chainElem.GetTitle()
                    for chainElem in tree.GetListOfFiles()]
        elif isinstance(tree, ROOT.TTree):
            # Retrieve the associated file
            treefile = tree.GetCurrentFile()
            if not treefile:
                # The tree has no associated input file
                return None
            else:
                return [treefile.GetName()]

        # RDataFrame may have been created with no input files
        return None

    def get_treeinfo(self):
        """
        Retrieve information about the tree used to construct this RDataFrame:
        the name, the list of files and info about friends.

        Returns:
            (TreeInfo, None): A `collections.namedtuple` holding information
            about the tree. If the dataframe was constructed without a TTree,
            returns None.
        """
        # Dispatch according to RDataFrame constructor
        firstarg = self.args[0]
        if isinstance(firstarg, int):
            # RDataFrame(ULong64_t numEntries)
            # No TTree was used to build the dataframe
            return None
        elif isinstance(firstarg, ROOT.TTree):
            # RDataFrame(TTree &tree, const ColumnNames_t &defaultBranches = {})
            # The first argument to the constructor is a TTree or TChain
            treefullpaths = ROOT.Internal.GetTreeFullPaths(firstarg)
            treeinputfiles = self.get_treeinputfiles(firstarg)
            treefriendinfo = ROOT.Internal.GetFriendInfo(firstarg)
            return TreeInfo(treefullpaths, treeinputfiles, treefriendinfo.fFriendNames, treefriendinfo.fFriendFileNames)
        elif isinstance(firstarg, str):
            # Get file(s) from second argument (may contain globbing characters)
            secondarg = self.args[1]
            if isinstance(secondarg, ROOT.TDirectory):
                # RDataFrame(treeName, dirPtr, defaultBranches = {})
                return TreeInfo([firstarg], [secondarg.GetName()], None, None)
            elif isinstance(secondarg, (ROOT.std.vector("string"), list)):
                # RDataFrame(treename, filenames, defaultBranches = {})
                return TreeInfo([firstarg], list(secondarg), None, None)
            elif isinstance(secondarg, str):
                # RDataFrame(treeName, filenameglob, defaultBranches = {})
                return TreeInfo([firstarg], glob.glob(secondarg), None, None)
        else:
            return None
