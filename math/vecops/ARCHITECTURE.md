# Architecture of the VecOps library

> This document is meant for ROOT developers, to quickly get their bearings around the VecOps library.

The main type in the library is `RVec`. Besides `RVec`, the library only contains helper types and functions.

`RVec` is a vector type that tries to be as `std::vector`-like as possible while adding a
few important features, namely:
- the ability to act as a view over an existing memory buffer (see "Memory adoption" below)
- a small-buffer optimization
- vectorized operator overloads
- a vectorized `operator[](mask)` to allow quick element selection together with vectorized operators
  (e.g. `etas[etas > k]` returns a new `RVec` with all elements greater than `k`)
- helper functions such as `InvariantMass`, `DeltaR`, `Argsort` are also provided

The current implementation of `RVec` is based on LLVM's SmallVector, extracted
from the head of LLVM's repo around December 2020.
We are not tracking the upstream implementation.

Compared to LLVM's SmallVectors:

- memory adoption capabilities have been added
- patches have been applied to make RVec work with (ROOT's version of) cppyy (notably `using` declarations had to be
  lowered in the inheritance hierarchy for cppyy to pick them up)
- `operator[](mask)` has been added, as well as several other "numpy-like" helper
  functions (these latter ones are free functions)
- logical operators `==`, `<`, `>` etc. return vectors rather than booleans
- the type of fSize and fCapacity is signed rather than unsigned, and fixed to 32 bits
- a number of minor patches have been applied for backward compatibility with the previous
  implementation of RVec (which did not have a small buffer optimization and was implemented
  in terms of `std::vector` with a custom allocator) and to make the code more consistent
  with ROOT's coding conventions

## RVec design

`SmallVectorBase`
   - `fBeginX`
   - `fSize`
   - `fCapacity`

   Basically the same as the corresponding LLVM class, with the template parameter removed: LLVM's SmallVectorBase
   is templated over the type of fSize and fCapacity. It contains the parts of `RVec` that do not depend on the value
   type.
   No other classes in the hierarchy can contain data members! We expect the memory after `SmallVectorBase` to be
   occupied by the small inline buffer.

`SmallVectorTemplateCommon<T>`
   - `getFirstEl()`: returns the address of the beginning of the small buffer
   - `begin()`, `end()`, `front()`, `back()`, etc.

   Basically the same as the corresponding LLVM class.
   It contains the parts that are independent of whether T is a POD or not.

`SmallVectorTemplateBase<T, bool TriviallyCopiable>` and the specialization `SmallVectorTemplateBase<T, true>`
   - `grow()`, `uninitialized_copy`, `uninitialized_move`, `push_back()`, `pop_back()`

   This class contains the parts of `RVec` that can be optimized for trivially copiable types.
   In particular, destruction can be skipped and memcpy can be used in place of copy/move construction.
   These optimizations are inherited from LLVM's SmallVector.

`RVecImpl<T>`
   The analogous of LLVM's `SmallVectorImpl`, it factors out of `RVec` the parts that are independent of
   the small buffer size, to limit the amount of code generated and provide a way to slice the small buffer
   size when passing around `RVec` objects.

`RVec<T, N = SensibleDefaultBufferSize<T>>`
   It aggregates `RVecImpl` and `SmallVectorStorage` (see below) through public inheritance.
   `N` is the small buffer size and defaults to a sensible value that depends on `sizeof(T)`.
   We expect most users to use the default and only rarely tweak the small buffer size.

### Helper types

- `SmallVectorAlignmentAndSize`: used to figure out the offset of the first small-buffer element in
  `SmallVectorTemplateCommon::getFirstEl`
- `SmallVectorStorage`: properly aligned "small buffer" storage. It's a separate type so that it can be specialized to
  be properly aligned also for the case of small buffer size = 0

## Memory adoption

We need RVec to be able to act as a view over an existing buffer rather than use its own
to save copies and allocations when reading ROOT data into `RVec`s, e.g. in `RDataFrame`.

The feature is exposed via a dedicated constructor: `RVec(pointer, size)`.
`RVec` then switches to its own storage as soon as a resize is requested.
`fCapacity == -1` indicates that we are in "memory adoption mode".

## Exception safety guarantees

As per [its docs](https://llvm.org/doxygen/classllvm_1_1SmallVector.html), LLVM's
`SmallVector` implementation "does not attempt to be exception-safe".
In its current implementation, `RVec` does not attempt to fix that.
This should not be a problem for `RVec`'s usecases (an exception thrown during
construction of an `RVec` typically means there is a bug to fix in the analysis code),
and we expect to be able to revisit the implementation and fix broken behavior if it
ever turns out to be problematic.

Relevant discussion on GitHub:

- https://github.com/root-project/root/pull/7502#issuecomment-818864506
- https://github.com/root-project/root/pull/7502#issuecomment-818905333
- https://github.com/root-project/root/pull/7502#issuecomment-821054757


