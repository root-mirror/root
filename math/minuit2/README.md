This is the Minuit2 fitter standalone edition, from the [ROOT] toolkit. It uses [CMake] 3.1+ to build.
For information about the Minuit2 fitter, please see the [documentation in ROOT][minuitdoc].

## Source

There are two ways to get Minuit2; you can checkout the [ROOT] source, then just build or use `add_subdirectory` with `<ROOT_SOURCE>/math/minuit2`, or you can get a Minuit2 source distribution which contains all the needed files to build with [CMake]. See [DEVELOP.md] for more information about extracting the source files from [ROOT].


## Building

To build, use the standard [CMake] procedure; on most systems, this looks like:

```bash
mkdir PATH_TO_MINIUT2_BUILD
cd PATH_TO_MINUIT2_BUILD
cmake PATH_TO_MINUIT2_SOURCE
cmake --build .
```

Of course, GUIs, IDEs, etc. that work with [CMake] will work with this package. The standard method of CMake building, with a build directory inside the Minuit2 source directory and using the makefile generator, would look like:

```bash
cd PATH_TO_MINUIT2_SOURCE
mkdir build
cd build
cmake ..
make
```


The standard [CMake] variables, such as `CMAKE_BUILD_TYPE` and `CMAKE_INSTALL_PREFIX`, work with Minuit2.  There are two other options:

* `minuit2_mpi` activates the (outdated C++) MPI bindings.
* `minuit2_omp` activates OpenMP (make sure all FCNs are threadsafe).

## Testing

You can run `ctest` or `make test` to run the Minuit2 test suite.

## Installing or using in another package

You can install the package using `cmake --build --target install .` (or `make install` if directly using the make system), or you can use it from the build directory. You can also include it in another CMake project using `add_subdirectory()` and linking to the `Minuit2` target. Since this package also exports targets, `find_package(Minuit2)` will also work once this package is built or installed. (For the curious, CMake adds a config script to `~/.cmake/packages` when building or
`$CMAKE_INSTALL_PREFIX/share/cmake/Modules` when installing a package that has export commands.)

To repeat; using this in your own CMake project usually amounts to:

```cmake
find_package(Minuit2)
# OR
add_subdirectory(Minuit2)

target_link_libraries(MyExeOrLib PUBLIC Minuit2::Minuit2)
```

You do not need to add include directories or anything else for Minuit2; the CMake target system handles all of this for you.

## Packaging

To build a binary package (add other generators with `-G`):
```bash
make package
```


[DEVELOP.md]: ./DEVELOP.md
[ROOT]: https://root.cern
[minuitdoc]: https://root.cern/root/htmldoc/guides/users-guide/ROOTUsersGuide.html#minuit2-package
[CMake]: https://cmake.org
