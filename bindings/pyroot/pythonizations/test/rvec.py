import unittest
import ROOT


class RVec(unittest.TestCase):
    fundamental_types = ['int', 'unsigned int', 'long', 'long long', 'Long64_t', 'unsigned long',
                         'unsigned long long', 'ULong64_t', 'float', 'double']
    extra_types = ['bool']


    def test_iteriter(self):
        '''
        Test the iteration over the iterator of an iterator

        This breaks if __iter__ or tp_iter is not defined for the iterator and causes
        issues, e.g., in comparison to numpy arrays.
        '''

        for dtype in self.fundamental_types + self.extra_types:
            rvec = ROOT.VecOps.RVec[dtype](3)
            it = iter(rvec)
            it2 = iter(it)
            self.assertEqual(it, it2)
            c = 0 # Prevent potential optimization of the loop
            for x in it2:
                c += 1
            self.assertEqual(c, 3)


if __name__ == '__main__':
    unittest.main()
