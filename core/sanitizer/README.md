## Configuration for the address sanitizer.
When built with `-Dasan=ON`, build flags for address sanitizer are added to ROOT's build setup. In this directory, an additional library is created
that holds default sanitizer configs for ROOT executables. It e.g. contains suppressions for leak sanitizer, which automatically runs with address
sanitizer.

When asan starts up, it checks if somebody defined the symbols in `SanitizerConfig.cxx`. In the standard asan runtime, these
functions are weak symbols, i.e. one can just override them with the desired configuration. That's what's happening here.

This can be achieved in two ways:
1. `LD_PRELOAD`: A micro library `libROOTSanitizerConfig.<dylib|so>` is created with the setup in this folder, and can be found in `<builddir>/lib/`.
   Loading it with `LD_PRELOAD` will bring ROOT's default sanitiser config into any non-sanitised exectuable, e.g. python.
2. ROOT exectuables will get the config automatically, using a static version of the config library, `libROOTStaticSanitizerConfig.a`.
   All ROOT executables statically link against it, so they start up without reporting lots of unfixable memory leaks (e.g. llvm).


#### Small linker magic to get the config symbols into ROOT's executables
When linking a ROOT exectuable, the setup functions from the sanitiser config library might get ignored, because they are not used in any of our exectuables.
In `cmake/modules/SetUp{Linux|MacOS}.cmake`, the functions are therefore marked as "undefined" for the linker, so it starts copying
them into all ROOT executables.
This way, root.exe, cling, ... can start up with a sane default config.


### Use your own address/leak sanitizer configuration
The default configurations can be overridden using the environment variables`ASAN_OPTIONS` and `LSAN_OPTIONS`. Refer to the
[address sanitizer documentation](https://github.com/google/sanitizers/wiki/AddressSanitizer) or use `ASAN_OPTIONS=help=1` when starting
up a sanitised exectuable (e.g. `root.exe`). A template for a leak suppression file can be found in `$ROOTSYS/etc/lsan-root.supp`.


## Create your own sanitised executable
ROOT exports a library target called `ROOT::ROOTStaticSanitizerConfig` that can be used to create sanitised executables with ROOT's default
address sanitizer config. Linking against this target will add the above setup functions and also add the address sanitizer flags that
ROOT was built with. It should be sufficient to
- Link all exectuables against this target
- And have `-fsanitize=address` in the `CXXFLAGS`.


## Use sanitised ROOT libraries from a non-sanitised executable (e.g. `python`)
When ROOT libraries are built with sanitizers, the address sanitizer runtime needs to be loaded at startup. However, when calling into ROOT
functions from python, that won't happen, since python is not sanitised. Therefore, the address sanitizer runtime has to be preloaded with
    LD_PRELOAD=<pathToRuntime>:libROOTSanitizerConfig.<so|dylib> pythonX ROOTScript.py

Preloading the shared sanitizer config as above is optional, but recommended, because it adds leak sanitizer suppressions.

On Mac, preloading is theoretically possible, but code signing and many other barriers might make it difficult.
