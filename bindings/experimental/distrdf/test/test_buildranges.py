from DistRDF.HeadNode import get_headnode
from DistRDF import Ranges
import warnings
import unittest


def emptysourceranges_to_tuples(ranges):
    """Convert EmptySourceRange objects to tuples with the shape (start, end)"""
    return [(r.start, r.end) for r in ranges]

def treeranges_to_tuples(ranges):
    """Convert TreeRange objects to tuples with the shape (start, end, filelist)"""
    return [(r.globalstart, r.globalend, r.localstarts, r.localends, r.filelist) for r in ranges]

def create_dummy_headnode(*args):
    """Create dummy head node instance needed in the test"""
    # Pass None as `npartitions`. The tests will modify this member
    # according to needs
    return get_headnode(None, *args)


class BuildRangesTest(unittest.TestCase):
    """
    Test the building of ranges with information from the head node of the graph
    adn the RangesBuilder class.
    """

    def test_nentries_multipleOf_npartitions(self):
        """
        `BuildRanges` method when the number of entries is a multiple of the
        number of partitions.

        """

        nentries_small = 10
        npartitions_small = 5
        nentries_large = 100
        npartitions_large = 10

        # First case
        rng = Ranges.get_balanced_ranges(nentries_small, npartitions_small)
        ranges_small = emptysourceranges_to_tuples(rng)

        # Second case
        rng = Ranges.get_balanced_ranges(nentries_large, npartitions_large)
        ranges_large = emptysourceranges_to_tuples(rng)

        ranges_small_reqd = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)]
        ranges_large_reqd = [
            (0, 10),
            (10, 20),
            (20, 30),
            (30, 40),
            (40, 50),
            (50, 60),
            (60, 70),
            (70, 80),
            (80, 90),
            (90, 100)
        ]

        self.assertListEqual(ranges_small, ranges_small_reqd)
        self.assertListEqual(ranges_large, ranges_large_reqd)

    def test_nentries_not_multipleOf_npartitions(self):
        """
        `BuildRanges` method when then number of entries is not a multiple of
        the number of partitions.

        """

        nentries_1 = 10
        nentries_2 = 9
        npartitions = 4

        # Example in which fractional part of
        # (nentries/npartitions) >= 0.5
        rng = Ranges.get_balanced_ranges(nentries_1, npartitions)
        ranges_1 = emptysourceranges_to_tuples(rng)

        # Example in which fractional part of
        # (nentries/npartitions) < 0.5
        rng = Ranges.get_balanced_ranges(nentries_2, npartitions)
        ranges_2 = emptysourceranges_to_tuples(rng)

        # Required output pairs
        ranges_1_reqd = [(0, 3), (3, 6), (6, 8), (8, 10)]
        ranges_2_reqd = [(0, 3), (3, 5), (5, 7), (7, 9)]

        self.assertListEqual(ranges_1, ranges_1_reqd)
        self.assertListEqual(ranges_2, ranges_2_reqd)

    def test_nentries_greater_than_npartitions(self):
        """
        `BuildRanges` method when the number of entries is smaller than the
        number of partitions.
        """

        nentries = 5
        npartitions = 7  # > nentries

        rng = Ranges.get_balanced_ranges(nentries, npartitions)
        ranges = emptysourceranges_to_tuples(rng)

        ranges_reqd = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]

        self.assertListEqual(ranges, ranges_reqd)

    def test_clustered_ranges_with_one_cluster(self):
        """
        Check that _get_clustered_ranges returns one range when the dataset
        contains a single cluster and the number of partitions is 1

        """

        treename = "TotemNtuple"
        filelist = ["backend/Slimmed_ntuple.root"]
        npartitions = 1
        clustersinfiles = Ranges.get_clusters(treename, filelist)
        friendinfo = None

        crs = Ranges.get_clustered_ranges(clustersinfiles, npartitions, friendinfo)
        ranges = treeranges_to_tuples(crs)

        ranges_reqd = [(0, 10, [0], [10], ["backend/Slimmed_ntuple.root"])]

        self.assertListEqual(ranges, ranges_reqd)

    def test_warning_when_npartitions_greater_than_clusters(self):
        """
        Check that _get_clustered_ranges raises a warning when the number of
        partitions is bigger than the number of clusters in the dataset.

        """

        treename = "TotemNtuple"
        filelist = ["backend/Slimmed_ntuple.root"]
        headnode = create_dummy_headnode(treename, filelist)
        headnode.npartitions = 2

        ranges_reqd = [(0, 10, [0], [10], ["backend/Slimmed_ntuple.root"])]

        with warnings.catch_warnings(record=True) as w:
            # Trigger warning
            crs = headnode.build_ranges()
            ranges = treeranges_to_tuples(crs)

            # Verify ranges
            self.assertListEqual(ranges, ranges_reqd)

            # Verify warning
            assert issubclass(w[-1].category, UserWarning)

    def test_clustered_ranges_with_two_clusters_two_partitions(self):
        """
        Check that _get_clustered_ranges creates clustered ranges respecting
        the cluster boundaries even if that implies to have ranges with very
        different numbers of entries.

        """

        treename = "myTree"
        filelist = ["backend/2clusters.root"]
        headnode = create_dummy_headnode(treename, filelist)
        headnode.npartitions = 2

        crs = headnode.build_ranges()
        ranges = treeranges_to_tuples(crs)

        ranges_reqd = [
            (0, 777, [0], [777], ["backend/2clusters.root"]),
            (777, 1000, [777], [1000], ["backend/2clusters.root"])
        ]

        self.assertListEqual(ranges, ranges_reqd)

    def test_clustered_ranges_with_four_clusters_four_partitions(self):
        """
        Check that _get_clustered_ranges creates clustered ranges as equal as
        possible for four partitions

        """

        treename = "myTree"
        filelist = ["backend/4clusters.root"]
        headnode = create_dummy_headnode(treename, filelist)
        headnode.npartitions = 4

        crs = headnode.build_ranges()
        ranges = treeranges_to_tuples(crs)

        ranges_reqd = [
            (0, 250, [0], [250],["backend/4clusters.root"]),
            (250, 500, [250], [500],["backend/4clusters.root"]),
            (500, 750, [500], [750],["backend/4clusters.root"]),
            (750, 1000, [750], [1000],["backend/4clusters.root"])
        ]

        self.assertListEqual(ranges, ranges_reqd)

    def test_clustered_ranges_with_many_clusters_four_partitions(self):
        """
        Check that _get_clustered_ranges creates clustered ranges as equal as
        possible for four partitions

        """

        treename = "myTree"
        filelist = ["backend/1000clusters.root"]
        headnode = create_dummy_headnode(treename, filelist)
        headnode.npartitions = 4

        crs = headnode.build_ranges()
        ranges = treeranges_to_tuples(crs)

        ranges_reqd = [
            (0, 250, [0], [250],["backend/1000clusters.root"]),
            (250, 500, [250], [500],["backend/1000clusters.root"]),
            (500, 750, [500], [750],["backend/1000clusters.root"]),
            (750, 1000, [750], [1000],["backend/1000clusters.root"])
        ]

        self.assertListEqual(ranges, ranges_reqd)

    def test_clustered_ranges_with_many_clusters_many_partitions(self):
        """
        Check that _get_clustered_ranges creates clustered ranges as equal as
        possible for the maximum number of possible partitions (number of
        clusters)

        """

        treename = "myTree"
        filelist = ["backend/1000clusters.root"]
        headnode = create_dummy_headnode(treename, filelist)
        headnode.npartitions = 1000

        crs = headnode.build_ranges()
        ranges = treeranges_to_tuples(crs)

        start = 0
        end = 1000
        step = 1

        ranges_reqd = [(a, b, [a], [b], filelist) for a, b in zip(range(start, end, step),
                                              range(step, end + 1, step))]

        self.assertListEqual(ranges, ranges_reqd)

    def test_buildranges_with_clustered_ranges(self):
        """
        Check that build_ranges produces clustered ranges when the dataset
        contains clusters.

        """
        headnode = create_dummy_headnode("myTree", "backend/1000clusters.root")
        headnode.npartitions = 1000

        crs = headnode.build_ranges()
        ranges = treeranges_to_tuples(crs)

        start = 0
        end = 1000
        step = 1

        ranges_reqd = [(a, b, [a], [b], ["backend/1000clusters.root"])
                        for a, b in zip(range(start, end, step), range(step, end + 1, step))]

        self.assertListEqual(ranges, ranges_reqd)

    def test_buildranges_with_balanced_ranges(self):
        """
        Check that build_ranges produces balanced ranges when there are no
        clusters involved.

        """
        headnode = create_dummy_headnode(50)
        headnode.npartitions = 16

        crs = headnode.build_ranges()
        ranges = emptysourceranges_to_tuples(crs)

        ranges_reqd = [
            (0, 4), (4, 8), (8, 11), (11, 14), (14, 17), (17, 20),
            (20, 23), (23, 26), (26, 29), (29, 32), (32, 35), (35, 38),
            (38, 41), (41, 44), (44, 47), (47, 50)
        ]

        self.assertListEqual(ranges, ranges_reqd)

    def test_clustered_ranges_with_two_files(self):
        """
        Create two ranges from two files with a different number of clusters.
        """

        treename = "myTree"
        filelist = ["backend/2clusters.root", "backend/4clusters.root"]
        headnode = create_dummy_headnode(treename, filelist)
        headnode.npartitions = 2

        crs = headnode.build_ranges()
        ranges = treeranges_to_tuples(crs)

        ranges_reqd = [
            (0, 1250, [0, 0], [1000, 250], ["backend/2clusters.root", "backend/4clusters.root"]),
            (250, 1000, [250], [1000], ["backend/4clusters.root"]),
        ]

        self.assertListEqual(ranges, ranges_reqd)
