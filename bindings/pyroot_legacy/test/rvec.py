import unittest
import ROOT


class RVec(unittest.TestCase):
    dtype = "float"

    def test_memoryadoption(self):
        x = ROOT.std.vector("float")(1)
        x[0] = 1
        y = ROOT.VecOps.RVec(self.dtype)(x.data(), x.size())
        self.assertEqual(x[0], y[0])

    def test_getset(self):
        x = ROOT.VecOps.RVec(self.dtype)(1)
        x[0] = 1
        self.assertEqual(x[0], 1)
        self.assertEqual(x.at(0), 1)

    def test_iter(self):
        x = ROOT.VecOps.RVec(self.dtype)(3)
        x[0], x[1], x[2] = 1, 2, 3
        for i, y in enumerate(x):
            self.assertEqual(y, x[i])

    def test_push_back(self):
        x = ROOT.VecOps.RVec(self.dtype)()
        x.push_back(1)
        self.assertEqual(x[0], 1)

    def test_len(self):
        x = ROOT.VecOps.RVec(self.dtype)(3)
        self.assertEqual(len(x), 3)


if __name__ == '__main__':
    unittest.main()
