This is the Minuit2 fitter standalone extractor, from the [ROOT] toolkit. It uses [CMake] 3.1+ to build.
See `README.md` for information about building Minuit2.

## Extracting from the ROOT source

To extract, run the following commands from the `math/minuit2/build` directory:

```bash
cmake .. -Dminuit2_standalone=ON
```

This will fill in the `math/minuit2` directory with all the files needed for Minuit2, copied from the corresponding ROOT files, as part of the configure step.
At this point, you could continue to build (using `make`). Note that the CMake option `minuit2_inroot` will automatically be set to `ON` if you are inside the ROOT source tree. Setting `minuit2_standalone` requires that this be inside the ROOT source tree. As always, any manual setting of a cached variable in CMake will be remembered as long as the `CMakeCache.txt` file is not removed.

Remember that after building a tarball or a binary package you should remove the copied files using:

```bash
make purge
```

Otherwise git shows the file as untracked, unless you explicitly remove their tracking yourself with a .gitignore file


## Building a tarball

Minuit2 standalone also has support for CPack to make installers for different platforms. To build a source package:

```bash
make package_source
```


This will create a source file in several formats that you can distribute. Reminder: You **must** have used `-Dminuit2_standalone=ON` when you configured CMake, or many of the files will be missing.

## Building a binary

To build a binary package (add other generators with `-G`):
```bash
make
make package
```

## Maintenance

If new files are needed by Minuit2 due to additions to [ROOT], they should be added to the source files lists in `src/Math/CMakeLists.txt` and `src/Minuit2/CMakeLists.txt` (depending on if it's a new Math or Minuit2 requirement).

For testing, the main `test/CMakeLists.txt` is used by ROOT, and the `test/*/CMakeLists.txt` files are used by the standalone build.

## How it works

Extracting from the ROOT sources is made possible through a few careful design features:

* A CMake variable `minuit2_inroot` lets the build system know we are inside ROOT (it looks for `../../build/version_info`)
* All files that are not part of the minuit2 directory are passed into `copy_standalone`, and that handles selecting the correct location
* `copy_standalone` copies the files into the minuit2 source directory if `minuit2_standalone` is `ON`

After this happens, all the standard CMake machinery can produce the source distribution. And, CMake correctly builds and installs in either mode, since all source and header files are explicitly listed.


[ROOT]: https://root.cern.ch
[minuitdoc]: https://root.cern.ch/root/htmldoc/guides/users-guide/ROOTUsersGuide.html#minuit2-package
[CMake]: https://cmake.org
