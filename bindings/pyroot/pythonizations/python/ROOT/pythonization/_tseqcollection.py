# Author: Enric Tejedor CERN  02/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
from libcppyy import SetOwnership

import sys


# Item access

def _check_type(idx, msg):
    # Parameters:
    # - idx: index whose type needs to be checked
    # - msg: message to show in case of type issue

    # Python2 also allows long indices
    if sys.version_info >= (3,0):
        allowed_types = (int,)
    else:
        allowed_types = (int, long)

    t = type(idx)
    if not t in allowed_types:
        raise TypeError(msg.format(t.__name__))

def _check_index(self, idx):
    # Parameters:
    # - self: collection
    # - idx: index to be checked
    # Returns:
    # - An index >= 0, equivalent to the input idx, which is verified
    # to be an integer and within the boundaries of the collection

    _check_type(idx, 'list indices must be integers, not {}')

    lcol = len(self)
    if idx < 0:
        idx += lcol

    if idx < 0 or idx >= lcol:
        raise IndexError('list assignment index out of range')

    return idx

def _getitem_pyz(self, idx):
    # Parameters:
    # - self: collection to get the item/s from
    # - idx: index/slice of the item/s
    # Returns:
    # - self[idx]

    # Slice
    if isinstance(idx, slice):
        res = self.__class__()
        indices = idx.indices(len(self))
        for i in range(*indices):
            res.Add(self.At(i))
    # Number
    else:
        idx = _check_index(self, idx)
        res = self.At(idx)

    return res

def _remove_at(self, idx):
    # Parameters:
    # - self: collection to remove the item from
    # - idx: index of the item, always positive
    lnk = self.FirstLink()
    for i in range(idx):
        lnk = lnk.Next()
    return self.Remove(lnk)

def _setitem_pyz(self, idx, val):
    # Parameters:
    # - self: collection where to set item/s
    # - idx: index/slice of the item/s
    # - val: value to assign

    # Slice
    if isinstance(idx, slice):
        # The value we assign has to be iterable
        try:
            _ = iter(val)
        except TypeError:
            raise TypeError('can only assign an iterable')

        indices = idx.indices(len(self))
        it = iter(range(*indices))
        for elem in val:
            # Prevent this new Python proxy from owning the C++ object
            # Otherwise we get an 'already deleted' error in
            # TList::Clear when the application ends
            SetOwnership(elem, False)
            try:
                i = next(it)
                self[i] = elem
            except StopIteration:
                # No more indices in range, just append
                self.append(elem)

        # If range is longer than the number of elements in val,
        # we need to remove the remaining elements of the range
        try:
            for i in it:
                del self[i]
        except StopIteration:
            # No more indices in range, we are done
            pass

    # Number
    else:
        idx = _check_index(self, idx)
        _remove_at(self, idx)
        self.AddAt(val, idx)

def _delitem_pyz(self, idx):
    # Parameters:
    # - self: collection to delete item from
    # - idx: index of the item

    # Slice
    if isinstance(idx, slice):
        indices = idx.indices(len(self))
        rg = range(*indices)

        step = indices[2]
        if step > 0:
            # Need to remove starting from the end
            rg = reversed(rg)

        for i in rg:
            _remove_at(self, i)
    # Number
    else:
        idx = _check_index(self, idx)
        _remove_at(self, idx)

# Python-list-like methods

def _insert_pyz(self, idx, val):
    # Parameters:
    # - self: collection where to insert a new item
    # - idx: index where to insert
    # - val: value to insert

    _check_type(idx, 'integer argument expected, got {}')

    # Check index
    lcol = len(self)
    if (idx < 0):
        idx += lcol

    if idx < 0:
        idx = 0
    elif idx > lcol:
        idx = lcol

    self.AddAt(val, idx)

def _pop_pyz(self, *args):
    # Parameters:
    # - self: collection where to pop an item from
    # - args: either empty or index to pop
    # Returns:
    # - If args is empty, it returns the last element of
    # the collection, else it returns the element for
    # which the index was specified.

    # Check arguments
    nargs = len(args)
    if nargs == 0:
        idx = len(self) - 1
    elif nargs > 1:
        raise TypeError('pop() takes at most 1 argument ({} given)'.format(nargs))
    else:
        idx = args[0]
        _check_type(idx, 'integer argument expected, got {}')

    if len(self) == 0:
        raise IndexError('pop from empty list')

    # Check index
    lcol = len(self)
    if idx < 0:
        idx += lcol

    if idx < 0 or idx >= lcol:
        raise IndexError('pop index out of range')

    return _remove_at(self, idx)

def _reverse_pyz(self):
    # Parameters:
    # - self: collection to be reversed

    if len(self) == 0:
        return

    t = tuple(self)
    self.Clear()
    for elem in t:
        self.AddAt(elem, 0)

def _sort_pyz(self, *args, **kwargs):
    # Parameters:
    # - self: collection to be reversed
    # - args: positional arguments
    # - kwargs: keyword arguments
    # For both args and kwargs, the Python sort()
    # arguments are accepted: key, reverse

    if len(self) == 0:
        return

    if not args and not kwargs:
        # No arguments -> rely on ROOT's Sort
        self.Sort()
    else:
        # Sort in a Python list copy
        l = list(self)
        l.sort(*args, **kwargs)
        self.Clear()
        self.extend(l)

def _index_pyz(self, val):
    # Parameters:
    # - self: collection
    # - val: element to find the index of
    # Returns:
    # - Index of the element in the collection

    idx = self.IndexOf(val)

    if idx < 0:
        raise ValueError('{} is not in list'.format(val))

    return idx


@pythonization()
def pythonize_tseqcollection(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name == 'TSeqCollection':
        # Item access methods
        klass.__getitem__ = _getitem_pyz
        klass.__setitem__ = _setitem_pyz
        klass.__delitem__ = _delitem_pyz

        # Python lists methods
        klass.insert  = _insert_pyz
        klass.pop     = _pop_pyz
        klass.reverse = _reverse_pyz
        klass.sort    = _sort_pyz
        klass.index   = _index_pyz

    return True
