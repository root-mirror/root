# Module.mk for utils module
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

# see also ModuleVars.mk

MODNAME := utils

ifneq ($(HOST),)

.PHONY: all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME):

else # ifneq ($(HOST),)

.PHONY: all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

.SECONDARY: $(ROOTCLINGTMPS)

CLINGMETAUTILSO    = $(METAUTILSTO) $(METAUTILSOLLVM)
ROOTCLINGEXEEXTRAO = $(COREO) $(COREDO) $(IOO) $(IODO) $(THREADO) $(THREADDO) $(METAOLLVM)
# The dependency on $(CLINGLIB) was added to prevent $(CLINGLIB) and
# $(ROOTCLINGEXE) from being linked in parallel.
$(ROOTCLINGEXE): $(ROOTCLINGO) $(ROOTCLINGUTILO) $(ROOTCLINGTCLINGO) \
	   $(CLINGMETAUTILSO) $(SNPRINTFO) $(CLINGO) $(ROOTCLINGEXEEXTRAO) \
           $(PCREDEP) $(CORELIBDEP) $(CLINGLIB)
	$(LD) $(LDFLAGS) $(OSTHREADLIBDIR) $(OSTHREADLIB) -o $@ $(ROOTCLINGO) $(ROOTCLINGUTILO) \
	   $(ROOTCLINGTCLINGO) $(CLINGMETAUTILSO) \
	   $(SNPRINTFO)  $(CLINGO) $(ROOTCLINGEXEEXTRAO) $(CLINGLIBEXTRA) \
	   $(RPATH) $(CILIBS) $(CORELIBEXTRA) $(PCRELDFLAGS) $(PCRELIB) \
	   $(CRYPTLIBS)

$(ROOTCLINGTMPEXE): $(CINTTMPO) $(ROOTCLINGTMPO) $(ROOTCLINGUTILO) \
	   $(METAUTILSO) $(CLINGMETAUTILSO) $(SNPRINTFO) $(STRLCPYO) $(CLINGO)
	$(LD) $(LDFLAGS) $(OSTHREADLIBDIR) $(OSTHREADLIB) -o $@ $(ROOTCLINGTMPO) $(ROOTCLINGUTILO) \
	   $(METAUTILSO) $(CLINGMETAUTILSO) $(SNPRINTFO) $(STRLCPYO) \
	   $(CINTTMPLIBS) $(CLINGO) $(CLINGLIBEXTRA) $(CILIBS)

$(ROOTCINTEXE): $(ROOTCLINGEXE)
	ln -f $(ROOTCLINGEXE) $(ROOTCINTEXE)

$(GENREFLEXEXE): $(ROOTCLINGEXE)
	ln -f $(ROOTCLINGEXE) $(GENREFLEXEXE)

all-$(MODNAME): $(ROOTCLINGTMPEXE) $(ROOTCLINGEXE) $(ROOTCINTEXE) \
                $(GENREFLEXEXE)

clean-$(MODNAME):
	@rm -f $(ROOTCLINGTMPO) $(ROOTCLINGO) $(ROOTCLINGUTILO)

clean:: clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)
	@rm -f $(ROOTCLINGDEP) $(ROOTCLINGTMPEXE) $(ROOTCLINGEXE) \
	   $(ROOTCINTEXE) $(GENREFLEXEXE) \
	   $(call stripsrc,$(UTILSDIRS)/*.exp $(UTILSDIRS)/*.lib \
	      $(UTILSDIRS)/*_tmp.cxx)

distclean:: distclean-$(MODNAME)

##### extra rules ######
$(call stripsrc,$(UTILSDIRS)/%_tmp.cxx): $(UTILSDIRS)/%.cxx
	$(MAKEDIR)
	cp $< $@

$(call stripsrc,$(UTILSDIRS)/rootcling_tmp.o): $(call stripsrc,\
	   $(UTILSDIRS)/rootcling_tmp.cxx)

$(call stripsrc,$(UTILSDIRS)/RStl_tmp.o): $(call stripsrc,\
	   $(UTILSDIRS)/RStl_tmp.cxx)

$(ROOTCLINGTMPO): $(LLVMDEP)
$(ROOTCLINGTMPO): CXXFLAGS += -UR__HAVE_CONFIG -DROOT_STAGE1_BUILD -I$(UTILSDIRS) -I$(METAUTILSDIRS) \
	   $(ROOTCLINGCXXFLAGS)
$(ROOTCLINGO): $(LLVMDEP)
$(ROOTCLINGO): CXXFLAGS += -UR__HAVE_CONFIG -I$(UTILSDIRS) -I$(METAUTILSDIRS) $(ROOTCLINGCXXFLAGS)
$(ROOTCLINGUTILO): $(LLVMDEP)
$(ROOTCLINGUTILO): CXXFLAGS += -UR__HAVE_CONFIG -I$(UTILSDIRS) -I$(METAUTILSDIRS) \
	   $(ROOTCLINGCXXFLAGS)
$(ROOTCLINGTCLINGO): CXXFLAGS += -I$(METADIRS)

# the -rdynamic flag is needed on cygwin to make symbols visible to dlsym
ifneq (,$(filter $(ARCH),win32gcc win64gcc))
$(ROOTCLINGEXE): LDFLAGS += -rdynamic
endif

endif # ifneq ($(HOST),)
