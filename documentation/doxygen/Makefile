
.PHONY: filter folders mathjax images doxygen

NJOB ?= $(shell nproc)
PYTHON_EXECUTABLE ?= python3
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

all: filter mathjax images doxygen replaceCollaborationDiagrams notebooks rootWork

filter:
	`root-config --cxx` -o filter filter.cxx -std=c++11 -O2

mathjax:
	$(call MkDir,$(DOXYGEN_IMAGE_PATH))
	tar xfz mathjax.tar.gz -C $(DOXYGEN_IMAGE_PATH)

images:
	$(call MkDir,$(DOXYGEN_IMAGE_PATH))
	cp images/* $(DOXYGEN_IMAGE_PATH)

pyzdoc:
	$(PYTHON_EXECUTABLE) extract_docstrings.py $(PYZ_DIR)

doxygen: filter pyzdoc
	rm -f tutorialWorklist_py tutorialWorklist_root
	$(call MkDir,$(DOXYGEN_OUTPUT_DIRECTORY))
	$(call MkDir,$(DOXYGEN_EXAMPLE_PATH))
	$(call MkDir,$(DOXYGEN_NOTEBOOK_PATH))
	./makehtmlfooter.sh > htmlfooter.html
	doxygen
	rm -rf files c1* *.ps *.png *.svg *.pdf *.root *.xpm *.out *.dat *.dtd *.dot *.txt listofclass.sh

collaborationDiagrams:
	bash ./makeCollaborationDiagrams.sh

replaceCollaborationDiagrams: doxygen collaborationDiagrams
	bash ./modifyClassWebpages.sh -j$(NJOB)

tutorialWorklist_py: doxygen
tutorialWorklist_root: doxygen

notebooks: tutorialWorklist_py makeNotebooks.sh
	cp ../../tutorials/unfold/*xml ../../tutorials/unfold/*dtd $(DOXYGEN_NOTEBOOK_PATH)/
	bash ./makeNotebooks.sh $< -j$(NJOB)
	rm -f $(DOXYGEN_NOTEBOOK_PATH)/*.root

rootWork: tutorialWorklist_root
	xargs -L 1 -P 12 root < $<

clean:
	rm -r $(DOXYGEN_OUTPUT_DIRECTORY) filter htmlfooter.html dataset* $(PYZ_DIR)/*.pyzdoc tutorialWorklist_py tutorialWorklist_root

