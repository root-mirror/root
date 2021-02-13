# Author: Enric Tejedor CERN  01/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
import cppyy


# Python-list-like methods

def _remove_pyz(self, o):
    # Parameters:
    # - self: collection
    # - o: object to remove from the collection
    res = self.Remove(o)

    if not res:
        raise ValueError('list.remove(x): x not in list')

def _extend_pyz(self, c):
    # Parameters:
    # - self: collection
    # - c: collection to extend self with
    it = iter(c)
    for _ in range(len(c)):
        self.Add(next(it))

def _count_pyz(self, o):
    # Parameters:
    # - self: collection
    # - o: object to be counted in the collection
    # Returns:
    # - Number of occurrences of the object in the collection
    n = 0

    for elem in self:
        if elem == o:
            n += 1

    return n

# Python operators

def _add_pyz(self, c):
    # Parameters:
    # - self: first collection to be added
    # - c: second collection to be added
    # Returns:
    # - self + c
    res = self.__class__()
    _extend_pyz(res, self)
    _extend_pyz(res, c)
    return res

def _mul_pyz(self, n):
    # Parameters:
    # - self: collection to be multiplied
    # - n: factor to multiply the collection by
    # Returns:
    # - self * n
    res = self.__class__()
    for _ in range(n):
        _extend_pyz(res, self)
    return res

def _imul_pyz(self, n):
    # Parameters:
    # - self: collection to be multiplied (in place)
    # - n: factor to multiply the collection by
    # Returns:
    # - self *= n
    c = self.__class__()
    c.AddAll(self)
    for _ in range(n - 1):
        _extend_pyz(self, c)
    return self

# Python iteration

def _iter_pyz(self):
    # Generator function to iterate on TCollections
    # Parameters:
    # - self: collection to be iterated
    it = cppyy.gbl.TIter(self)
    # TIter instances are iterable
    for o in it:
        yield o


@pythonization()
def pythonize_tcollection(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name == 'TCollection':
        # Support `len(c)` as `c.GetEntries()`
        klass.__len__ = klass.GetEntries

        # Add Python lists methods
        klass.append = klass.Add
        klass.remove = _remove_pyz
        klass.extend = _extend_pyz
        klass.count = _count_pyz

        # Define Python operators
        klass.__add__ = _add_pyz
        klass.__mul__ = _mul_pyz
        klass.__rmul__ = _mul_pyz
        klass.__imul__ = _imul_pyz

        # Make TCollections iterable
        klass.__iter__ = _iter_pyz

    return True
