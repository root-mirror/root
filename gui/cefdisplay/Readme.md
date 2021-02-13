\page cefpage ROOT installation with CEF

## Compilation with CEF support

See details about [Chromimum Embedded Framework](https://bitbucket.org/chromiumembedded/cef)

1. Current code tested with CEF3 branch 4147, Chromium 84 (August 2020)
   Older CEF versions are no longer supported.

2. Download binary code from [http://opensource.spotify.com/cefbuilds/index.html](http://opensource.spotify.com/cefbuilds/index.html) and unpack it in directory without spaces and special symbols:

~~~
     $ mkdir /d/cef
     $ cd /d/cef/
     $ wget http://opensource.spotify.com/cefbuilds/cef_binary_84.4.1+gfdc7504+chromium-84.0.4147.105_linux32.tar.bz2
     $ tar xjf cef_binary_84.4.1+gfdc7504+chromium-84.0.4147.105_linux32.tar.bz2
~~~


3. Set `CEF_ROOT` shell variable to unpacked directory:

~~~
     $ export CEF_ROOT=/d/cef/cef_binary_84.4.1+gfdc7504+chromium-84.0.4147.105_linux64
~~~


4. As it is on 21.08.2020, CEF has problem to compile with latest gcc, issuing error because of warnings in some standard libraries.
   Therefore one has to modify `$CEF_ROOT/cmake/cef_variables.cmake`, removing code at line ~90 for Linux platforms:

~~~
   -Werror                         # Treat warnings as errors
~~~

5. Install prerequisites - see comments in `$CEF_ROOT/CMakeLists.txt`.
   For the linux these are: `build-essential`, `libgtk2.0-dev`, `libgtkglext1-dev`

6. Compile CEF to produce `libcef_dll_wrapper`:

~~~
     $ cd $CEF_ROOT
     $ mkdir build
     $ cd build
     $ cmake $CEF_ROOT
     $ make -j libcef_dll_wrapper cefsimple
~~~

7. When configure ROOT compilation with `cmake -Dwebgui=ON -Dcefweb=ON ...`, CEF_ROOT shell variable should be set appropriately.
   During compilation library `$ROOTSYS/lib/libROOTCefDisplay.so` and executable `$ROOTSYS/bin/cef_main`
   should be created. Also check that several files like `icudtl.dat`, `v8_context_snapshot_blob.bin`, `snapshot_blob.bin`
   copied into ROOT library directory

8. Run ROOT with `--web=cef` argument to use CEF web display like

~~~
   $ root --web=cef $ROOTSYS/tutorials/v7/draw_rh2.cxx
~~~


## Compile libcef_dll_wrapper on Windows

1. Download binary win32 build like http://opensource.spotify.com/cefbuilds/cef_binary_85.2.11+g0202816+chromium-85.0.4183.83_windows32.tar.bz2

2. Extract in directory without spaces like `C:\Soft\cef`

3. Modify `cmake/cef_variables.cmake` to set dynamic linking, replace "/MT" by "/MD" in approx line 389

4. Start "x86 Native tools Command Prompt for VS 2019". Do:
~~~
   $ cd C:\Soft\cef
   $ mkdir build
   $ cd build
   $ cmake -G"Visual Studio 16 2019" -A Win32 -Thost=x64 ..
   $ cmake --build . --config Release --target libcef_dll_wrapper
~~~

5. Before compiling ROOT, `set CEF_ROOT=C:\Soft\cef` variable


## Using plain CEF in ROOT batch mode on Linux

Default CEF builds, provided by [http://opensource.spotify.com/cefbuilds/index.html](http://opensource.spotify.com/cefbuilds/index.html), do
not include support of Ozone framework, which the only support headless mode in CEF. To run ROOT in headless (or batch) made with such CEF distribution,
one can use `Xvfb` server. Most simple way is to use `xvfb-run` utility like:

~~~
      $ xvfb-run --server-args='-screen 0, 1024x768x16'  root.exe -l --web=cef $ROOTSYS/tutorials/v7/line.cxx -q
~~~

Or run `Xvfb` before starting ROOT:

~~~
     $ Xvfb :99 &
     $ export DISPLAY=:99
     $ root.exe -l --web=cef $ROOTSYS/tutorials/v7/line.cxx -q
~~~


## Compile CEF with ozone support

Since March 2019 one can compile [CEF without X11](https://bitbucket.org/chromiumembedded/cef/issues/2296/), but such builds not provided.
Therefore to be able to use real headless mode in CEF, one should compile it from sources.
On [CEF build tutorial](https://bitbucket.org/chromiumembedded/cef/wiki/AutomatedBuildSetup.md) one can find complete compilation documentation.
Several Ubuntu distributions are supported by CEF, all others may require extra work. Once all depndencies are installed,
CEF with ozone support can be compiled with following commands:

~~~
   $ export GN_DEFINES="is_official_build=true use_sysroot=true use_allocator=none symbol_level=1 is_cfi=false use_thin_lto=false use_ozone=true"
   $ python automate-git.py --download-dir=/home/user/cef --branch=4147 --minimal-distrib --client-distrib --force-clean --x64-build --build-target=cefsimple
~~~

With little luck one get prepared tarballs in `/home/user/cef/chromium/src/cef/binary_distrib`.
Just install it in the same way as described before in this document.
ROOT will automatically detect that CEF build with `ozone` support and will use it for both interactive and headless modes.


