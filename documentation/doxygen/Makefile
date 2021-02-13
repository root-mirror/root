
.PHONY: filter folders mathjax images doxygen

PYTHON_EXECUTABLE ?= /usr/bin/python
export PYTHON_EXECUTABLE

DOXYGEN_OUTPUT_DIRECTORY ?= $(HOME)/rootdoc
export DOXYGEN_OUTPUT_DIRECTORY
export DOXYGEN_SOURCE_DIRECTORY := ../..
export DOXYGEN_ROOT_VERSION := $(shell root -l -b -q -e 'printf("%s\n",gROOT->GetGitBranch());')
export DOXYGEN_ROOT_VERSION := $(subst -00-patches,,$(DOXYGEN_ROOT_VERSION))
export DOXYGEN_ROOT_VERSION := $(subst heads/,,$(DOXYGEN_ROOT_VERSION))
export DOXYGEN_EXAMPLE_PATH := $(DOXYGEN_OUTPUT_DIRECTORY)/macros
export DOXYGEN_IMAGE_PATH := $(DOXYGEN_OUTPUT_DIRECTORY)/html
export DOXYGEN_NOTEBOOK_PATH := $(DOXYGEN_OUTPUT_DIRECTORY)/notebooks
export DOXYGEN_STRIP_PATH := $(shell cd $(PWD)/../..; pwd)
export DOXYGEN_INCLUDE_PATH := $(shell find $(DOXYGEN_STRIP_PATH) -type d -name dictpch -prune -o -type d -name inc)

PYZ_DIR = ../../bindings/pyroot/pythonizations/python/ROOT/pythonization

define MkDir
	+@[ -d $1 ] || mkdir -p $1
endef

all: filter folders mathjax images doxygen

filter:
	`root-config --cxx` -o filter filter.cxx -std=c++11

mathjax:
	$(call MkDir,$(DOXYGEN_IMAGE_PATH))
	tar xfz mathjax.tar.gz -C $(DOXYGEN_IMAGE_PATH)

images:
	$(call MkDir,$(DOXYGEN_IMAGE_PATH))
	cp images/* $(DOXYGEN_IMAGE_PATH)

pyzdoc:
	$(PYTHON_EXECUTABLE) extract_docstrings.py $(PYZ_DIR)

doxygen: filter pyzdoc
	$(call MkDir,$(DOXYGEN_OUTPUT_DIRECTORY))
	$(call MkDir,$(DOXYGEN_EXAMPLE_PATH))
	$(call MkDir,$(DOXYGEN_NOTEBOOK_PATH))
	./makehtmlfooter.sh > htmlfooter.html
	doxygen
	./listlibs.sh
	rm -rf files c1* *.ps *.png *.svg *.pdf *.root *.xpm *.out *.dat *.dtd *.dot *.txt listofclass.sh
	rm -f $(DOXYGEN_NOTEBOOK_PATH)/*.root

clean:
	rm -r $(DOXYGEN_OUTPUT_DIRECTORY) filter htmlfooter.html dataset* $(PYZ_DIR)/*.pyzdoc
