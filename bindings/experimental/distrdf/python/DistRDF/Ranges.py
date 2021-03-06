import logging

from functools import total_ordering

import ROOT

logger = logging.getLogger(__name__)


class EmptySourceRange(object):
    """
    Empty source range of entries

    Attributes:

    id (int): Sequential counter to identify this range. It is used to assign
        a filename to a partial Snapshot in case it was requested. The id is
        assigned in `get_balanced_ranges` to ensure that each distributed
        RDataFrame run has a list of ranges with sequential ids starting from
        zero.

    start (int): Starting entry of this range.

    end (int): Ending entry of this range.
    """

    def __init__(self, rangeid, start, end):
        """set attributes"""
        self.id = rangeid
        self.start = start
        self.end = end

class TreeRange(object):
    """
    TTree range of entries. The entries are local with respect to the list of
    files that are processed with this range. These files are a subset of the
    global list of input files of the original dataset.

    Attributes:

    id (int): Sequential counter to identify this range. It is used to assign
        a filename to a partial Snapshot in case it was requested. The id is
        assigned in `get_clustered_ranges` to ensure that each distributed
        RDataFrame run has a list of ranges with sequential ids starting from
        zero.


    start (int): Starting entry of this range.

    end (int): Ending entry of this range.

    filelist (list[str]): List of files to be processed with this range.

    friend_info (ROOT.Internal.TreeUtils.RFriendInfo): Information about friend trees.
    """

    def __init__(self, rangeid, start, end, filelist, friend_info):
        """set attributes"""
        self.id = rangeid
        self.start = start
        self.end = end
        self.filelist = filelist
        self.friend_info = friend_info


@total_ordering
class FileAndIndex(object):
    """
    This is a pair (filename, index) that represents the index of the current
    filename in the list of input files of the dataset

    Attributes:

    filename (str): The name of the file.

    index (int): The index of the file in the list of input files.
    """

    def __init__(self, filename, fileindex):
        """set attributes"""
        self.filename = filename
        self.fileindex = fileindex

    def __lt__(self, other):
        """Total ordering is defined to be able to sort a set[FileAndIndex]."""
        return self.fileindex < other.fileindex

    def __eq__(self, other):
        """Defined in compliance with total_ordering decorator"""
        return self.fileindex == other.fileindex

    def __hash__(self):
        """
        This type needs to be hashable to be part of a set in
        get_clustered_ranges
        """
        return hash(self.filename)


@total_ordering
class ChainCluster(object):
    """
    Descriptor of a cluster of entries in a TChain. Uses global entries rather
    than local.

    Attributes:

    start (int): The starting global entry of this cluster in the chain.

    end (int): The ending global entry of this cluster in the chain.

    offset (int): The offset of this cluster in the chain. That is, the starting
        entry of the file this cluster belongs to in the chain.

    filetuple (FileAndIndex): A pair with the name of the file this cluster
        belongs to and the index of that file in the chain.
    """

    def __init__(self, start, end, offset, filetuple):
        """set attributes"""
        self.start = start
        self.end = end
        self.offset = offset
        self.filetuple = filetuple

    def __lt__(self, other):
        """
        In `get_clustered_ranges` we need to retrieve the minimum and maximum
        entries in a certain list of clusters.
        """
        return self.start < other.start and self.end < other.end

    def __eq__(self, other):
        """Defined in compliance with total_ordering decorator"""
        return (self.start == other.start and
               self.end == other.end and
               self.offset == other.offset and
               self.filetuple.filename == other.filetuple.filename and
               self.filetuple.fileindex == other.filetuple.fileindex)


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


def get_clusters(treename, filelist):
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
            clusters.append(ChainCluster(start + offset, end + offset, offset,
                                    FileAndIndex(filename, fileindex)))
            start = end

        fileindex += 1
        offset += entries

    logger.debug("Returning files with their clusters:\n%s",
                 "\n\n".join(map(str, clusters)))

    return clusters


def get_balanced_ranges(nentries, npartitions):
    """
    Builds range pairs from the given values of the number of entries in
    the dataset and number of partitions required. Each range contains the
    same amount of entries, except for those cases where the number of
    entries is not a multiple of the partitions.

    Args:
        nentries (int): The number of entries in a dataset.

        npartitions (int): The number of partititions the sequence of entries
            should be split in.

    Returns:
        list[DistRDF.Ranges.EmptySourceRange]: Each element of the list contains
            the start and end entry of the corresponding range.
    """

    partition_size = nentries // npartitions

    i = 0  # Iterator

    ranges = []

    remainder = nentries % npartitions

    rangeid = 0 # Keep track of the current range id
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

        ranges.append(EmptySourceRange(rangeid, start, end))
        rangeid += 1

    return ranges


def get_clustered_ranges(clustersinfiles, npartitions, treename, friend_info):
    """
    Builds ``TreeRange`` objects taking into account the clusters of the
    dataset. Each range will represent the entries processed within a single
    partition of the distributed dataset.

    Args:
        clustersinfiles (list): List of namedtuples representing clusters in
            the input files of the current dataset.

        npartitions (int): Number of ranges that will be produced.

        treename (str): Name of the tree.

        friend_info (ROOT.Internal.TreeUtils.RFriendInfo): Information about friend
            trees.

    Returns:
        list[DistRDF.Ranges.TreeRange]: Each element of the list represents one
            range in which the dataset has been split for distributed execution.
            Each `TreeRange` contains a starting entry, an ending entry, the
            list of files that are traversed to get all the entries and
            information about friend trees::

            [
                TreeRange(start=0,
                    end=1500,
                    filelist=['filename_1.root',
                              'filename_2.root'],
                    friend_info=None),
                TreeRange(start=1500,
                    end=3000,
                    filelist=['filename_2..root',
                              'filename_3.root'],
                    friend_info=None)
            ]

    """

    """
    This list comprehension builds ``TreeRange`` tuples with the following
    elements:
    1. ``start``: The minimum entry among all the clusters considered in a
        given partition. The offset of the first cluster of the list is
        subtracted. This is useful to keep the reference of the range with
        respect to the current files (see below).
    2. ``end``: The maximum entry among all the clusters considered in a
        given partition. The offset of the first cluster of the list is
        subtracted. This is useful to keep the reference of the range with
        respect to the current files (see below).
    3. ``filelist``: The list of files that are span between entries
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
    4. ``friend_info``: Information about friend trees.
    In each range, the offset of the first file is always subtracted to the
    ``start`` and ``end`` entries. This is needed to maintain a reference of
    the entries of the range with respect to the list of files that hold
    them. For example, given the following files::
        tree10000entries10clusters.root --> 10000 entries, 10 clusters
        tree20000entries10clusters.root --> 20000 entries, 10 clusters
        tree30000entries10clusters.root --> 30000 entries, 10 clusters
    Building 2 ranges will lead to the following tuples::
        TreeRange(start=0,
                end=20000,
                filelist=['tree10000entries10clusters.root',
                        'tree20000entries10clusters.root'],
                friend_info=None)
        TreeRange(start=10000,
                end=50000,
                filelist=['tree20000entries10clusters.root',
                        'tree30000entries10clusters.root'],
                friend_info=None)
    The first ``TreeRange`` will read the first 10000 entries from the first
    file, then switch to the second file and read the first 10000 entries.
    The second ``TreeRange`` will start from entry number 10000 of the second
    file up until the end of that file (entry number 20000), then switch to
    the third file and read the whole 30000 entries there.
    """

    # TODO: Make this passage more clear. Maybe split the comprehension in more
    # parts.
    clustered_ranges = [
        TreeRange(
            rangeid,  # type: int
            min(clusters).start - clusters[0].offset,  # type: int
            max(clusters).end - clusters[0].offset,  # type: int
            [
                filetuple.filename
                for filetuple in sorted(set([
                    cluster.filetuple for cluster in clusters
                ]), key=lambda curtuple: curtuple.fileindex)
            ],  # type: list[str]
            friend_info  # type: ROOT.Internal.TreeUtils.RFriendInfo
        )  # type: DistRDF.Ranges.TreeRange
        for rangeid, clusters in enumerate(_n_even_chunks(clustersinfiles, npartitions))
    ]

    logger.debug("Created following clustered ranges:\n%s",
                 "\n\n".join(map(str, clustered_ranges)))

    return clustered_ranges
