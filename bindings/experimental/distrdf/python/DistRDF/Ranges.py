import collections
import logging

import ROOT

logger = logging.getLogger(__name__)

EmptySourceRange = collections.namedtuple("EmptySourceRange", ["start", "end"])
TreeRange = collections.namedtuple("TreeRange", ["start", "end", "filelist", "friend_info"])


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
        list: List of `EmptySourceRange` namedtuples. Each tuple contains the
            start and end entry of the corresponding range.
    """
    partition_size = nentries // npartitions

    i = 0  # Iterator

    ranges = []

    remainder = nentries % npartitions

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

        ranges.append(EmptySourceRange(start, end))

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

        friend_info (DistRDF.HeadNode.FriendInfo): Information about friend
            trees.

    Returns:
        list: List of `TreeRange` namedtuples. Each tuple represents one range
            in which the dataset has been split for distributed execution.
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
    clustered_ranges = [
        TreeRange(
            min(clusters)[0] - clusters[0].offset,  # type: int
            max(clusters)[1] - clusters[0].offset,  # type: int
            [
                filetuple.filename
                for filetuple in sorted(set([
                    cluster.filetuple for cluster in clusters
                ]), key=lambda curtuple: curtuple[1])
            ],  # type: list[str]
            friend_info  # type: DistRDF.HeadNode.FriendInfo
        )  # type: collections.namedtuple
        for clusters in _n_even_chunks(clustersinfiles, npartitions)
    ]

    logger.debug("Created following clustered ranges:\n%s",
                 "\n\n".join(map(str, clustered_ranges)))

    return clustered_ranges
