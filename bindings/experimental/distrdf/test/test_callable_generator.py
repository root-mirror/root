import unittest

from DistRDF import ComputationGraphGenerator, Node, HeadNode, Proxy
from DistRDF.Backends import Base


class ComputationGraphGeneratorTest(unittest.TestCase):
    """
    Check mechanism to create a callable function that returns a PyROOT object
    per each DistRDF graph node. This callable takes care of the grape pruning.
    """

    class TestBackend(Base.BaseBackend):
        """Dummy backend."""

        def ProcessAndMerge(self, ranges, mapper, reducer):
            """Dummy implementation of ProcessAndMerge."""
            pass

        def distribute_unique_paths(self, includes_list):
            """
            Dummy implementation of distribute_unique_paths. Does nothing.
            """
            pass

        def make_dataframe(self, *args, **kwargs):
            """Dummy make_dataframe"""
            pass

    class Temp(object):
        """A Class for mocking RDF CPP object."""

        def __init__(self):
            """
            Creates a mock instance. Each mock method adds an unique number to
            the `ord_list` so we can check the order in which they were called.
            """
            self.ord_list = []

        def Define(self):
            """Mock Define method"""
            self.ord_list.append(1)
            return self

        def Filter(self):
            """Mock Filter method"""
            self.ord_list.append(2)
            return self

        def Count(self):
            """Mock Count method"""
            self.ord_list.append(3)
            return self

    def test_mapper_from_graph(self):
        """A simple test case to check the working of mapper."""
        # A mock RDF object
        t = ComputationGraphGeneratorTest.Temp()

        # Head node
        hn = HeadNode.Factory.get_headnode(1)
        hn.backend = ComputationGraphGeneratorTest.TestBackend()
        node = Proxy.TransformationProxy(hn)
        # Set of operations to build the graph
        n1 = node.Define()
        n2 = node.Filter().Filter()
        n4 = n2.Count()
        n5 = n1.Count()
        n6 = node.Filter()  # noqa: avoid PEP8 F841

        # Generate and execute the mapper
        generator = ComputationGraphGenerator.ComputationGraphGenerator(
            node.proxied_node)
        mapper_func = generator.get_callable()
        values = mapper_func(t)
        nodes = generator.get_action_nodes()

        reqd_order = [1, 3, 2, 2, 3, 2]

        self.assertEqual(t.ord_list, reqd_order)
        self.assertListEqual(nodes, [n5.proxied_node, n4.proxied_node])
        self.assertListEqual(values, [t, t])

    def test_mapper_with_pruning(self):
        """
        A test case to check that the mapper works even in the case of
        pruning.

        """
        # A mock RDF object
        t = ComputationGraphGeneratorTest.Temp()

        # Head node
        hn = HeadNode.Factory.get_headnode(1)
        hn.backend = ComputationGraphGeneratorTest.TestBackend()
        node = Proxy.TransformationProxy(hn)

        # Set of operations to build the graph
        n1 = node.Define()
        n2 = node.Filter().Filter()
        n4 = n2.Count()
        n5 = n1.Count()
        n6 = node.Filter()  # noqa: avoid PEP8 F841

        # Reason for pruning (change of reference)
        n5 = n1.Filter()  # noqa: avoid PEP8 F841

        # Generate and execute the mapper
        generator = ComputationGraphGenerator.ComputationGraphGenerator(
            node.proxied_node)
        mapper_func = generator.get_callable()
        values = mapper_func(t)
        nodes = generator.get_action_nodes()

        reqd_order = [1, 2, 2, 2, 3, 2]

        self.assertEqual(t.ord_list, reqd_order)
        self.assertListEqual(nodes, [n4.proxied_node])
        self.assertListEqual(values, [t])
