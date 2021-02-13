dicts   := functioncallsDict.so runvectorDict.so
libs    := libfunctioncalls.so librunvector.so
modules :=
execs   := functioncalls_main runvector_main
CHECK_PYBIND11 := $(shell python -m pybind11 1>/dev/null 2>/dev/null; echo $$?)
ifeq ($(CHECK_PYBIND11),0)
  modules += py11_functioncalls.so
endif
CHECK_SWIG := $(shell swig -help 1>/dev/null 2>/dev/null; echo $$?)
ifeq ($(CHECK_SWIG),0)
  modules += _swig_functioncalls.so swig_functioncalls.py
  modules += _swig_runvector.so swig_runvector.py
endif
all : $(dicts) $(libs) $(execs) $(modules)

ifeq ($(ROOTSYS),)
  clingconfig := cling-config
else
  clingconfig := root-config
endif

cppflags:=-O3 -fPIC
cppyy_cppflags=$(shell $(clingconfig) --cflags) $(cppflags)

PLATFORM := $(shell uname -s)
ifeq ($(PLATFORM),Darwin)
  cppflags+=-dynamiclib -single_module -arch x86_64 -undefined dynamic_lookup
endif

py11_%.so: py11_%.cxx lib%.so
	$(CXX) $(cppflags) -shared -std=c++11 -fPIC `python -m pybind11 --includes` py11_$*.cxx -o $@ -L. -l$*

swig_%.py: %.h %.i
	swig -python -c++ -builtin $*.i

_swig_%.so: swig_%.py lib%.so
	$(CXX) $(cppflags) -shared `python-config --includes` -std=c++11 -fPIC $*_wrap.cxx -o $@ -L. -l$*

%Dict.so: %_rflx.cpp lib%.so
	$(CXX) $(cppyy_cppflags) -shared -o $@ $*_rflx.cpp -L. -l$*

%_rflx.cpp: %.h %.xml
	genreflex $< --selection=$*.xml --rootmap=$*Dict.rootmap --rootmap-lib=$*Dict.so

lib%.so: %.cxx
	$(CXX) $(cppflags) -shared -o lib$*.so $*.cxx

# note -O2 for cxx to make sure code actually runs instead of being optimized out of existence
%_main: %_main.cxx lib%.so
	$(CXX) -O2 -fPIC -o $@ $*_main.cxx -L. -l$*

.PHONY: bench clean

bench: all
	pytest -s bench_runvector.py --benchmark-sort=mean

clean:
	-rm -f $(dicts) $(libs) $(execs) $(modules) $(wildcard *.rootmap) $(wildcard *_rdict.pcm)
