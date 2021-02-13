import unittest
import array

import ROOT


class TGraphGetters(unittest.TestCase):
    """
    Test for the pythonization of TGraph, TGraph2D and their error
    subclasses, in particular of their X,Y,Z coordinates and errors
    getters, which sets the size of the returned buffers.
    """

    # Tests
    def test_graph(self):
        N = 5
        xval, yval = 1, 2

        ax = array.array('d', map(lambda x: x*xval, range(N)))
        ay = array.array('d', map(lambda x: x*yval, range(N)))

        g = ROOT.TGraph(N, ax, ay)

        # x and y are buffers of doubles
        x = g.GetX()
        y = g.GetY()

        # We can get the size of the buffers
        self.assertEqual(len(x), N)
        self.assertEqual(len(y), N)

        # The buffers are iterable
        self.assertEqual(list(x), list(ax))
        self.assertEqual(list(y), list(ay))

    def test_graph2derrors(self):
        N = 5
        xval, yval, zval = 1, 2, 3
        xerrval, yerrval, zerrval = 0.1, 0.2, 0.3

        ax = array.array('d', map(lambda x: x*xval, range(N)))
        ay = array.array('d', map(lambda x: x*yval, range(N)))
        az = array.array('d', map(lambda x: x*zval, range(N)))
        aex = array.array('d', map(lambda x: x*xerrval, range(N)))
        aey = array.array('d', map(lambda x: x*yerrval, range(N)))
        aez = array.array('d', map(lambda x: x*zerrval, range(N)))

        g = ROOT.TGraph2DErrors(N, ax, ay, az, aex, aey, aez)

        # x, y, z, ex, ey and ez are buffers of doubles
        x = g.GetX()
        y = g.GetY()
        z = g.GetZ()
        ex = g.GetEX()
        ey = g.GetEY()
        ez = g.GetEZ()

        # We can get the size of the buffers
        self.assertEqual(len(x), N)
        self.assertEqual(len(y), N)
        self.assertEqual(len(z), N)
        self.assertEqual(len(ex), N)
        self.assertEqual(len(ey), N)
        self.assertEqual(len(ez), N)

        # The buffers are iterable
        self.assertEqual(list(x), list(ax))
        self.assertEqual(list(y), list(ay))
        self.assertEqual(list(z), list(az))
        self.assertEqual(list(ex), list(aex))
        self.assertEqual(list(ey), list(aey))
        self.assertEqual(list(ez), list(aez))


if __name__ == '__main__':
    unittest.main()
