# Module.mk for sql module
# Copyright (c) 2005 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 7/12/2005

MODNAME      := sql
MODDIR       := $(ROOT_SRCDIR)/io/$(MODNAME)
MODDIRS      := $(MODDIR)/src
MODDIRI      := $(MODDIR)/inc

SQLDIR       := $(MODDIR)
SQLDIRS      := $(SQLDIR)/src
SQLDIRI      := $(SQLDIR)/inc

##### libSQL #####
SQLL         := $(MODDIRI)/LinkDef.h
SQLDS        := $(call stripsrc,$(MODDIRS)/G__SQLIO.cxx)
SQLDO        := $(SQLDS:.cxx=.o)
SQLDH        := $(SQLDS:.cxx=.h)

SQLH         := $(filter-out $(MODDIRI)/LinkDef%,$(wildcard $(MODDIRI)/*.h))
SQLS         := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.cxx))
SQLO         := $(call stripsrc,$(SQLS:.cxx=.o))

SQLDEP       := $(SQLO:.o=.d) $(SQLDO:.o=.d)

SQLLIB       := $(LPATH)/libSQLIO.$(SOEXT)
SQLMAP       := $(SQLLIB:.$(SOEXT)=.rootmap)

# used in the main Makefile
SQLH_REL     := $(patsubst $(MODDIRI)/%.h,include/%.h,$(SQLH))
ALLHDRS      += $(SQLH_REL)
ALLLIBS      += $(SQLLIB)
ALLMAPS      += $(SQLMAP)
ifeq ($(CXXMODULES),yes)
  CXXMODULES_HEADERS := $(patsubst include/%,header \"%\"\\n,$(SQLH_REL))
  CXXMODULES_MODULEMAP_CONTENTS += module Io_$(MODNAME) { \\n
  CXXMODULES_MODULEMAP_CONTENTS += $(CXXMODULES_HEADERS)
  CXXMODULES_MODULEMAP_CONTENTS += "export * \\n"
  CXXMODULES_MODULEMAP_CONTENTS += link \"$(SQLLIB)\" \\n
  CXXMODULES_MODULEMAP_CONTENTS += } \\n
endif

# include all dependency files
INCLUDEFILES += $(SQLDEP)

##### local rules #####
.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

include/%.h:    $(SQLDIRI)/%.h
		cp $< $@

$(SQLLIB):      $(SQLO) $(SQLDO) $(ORDER_) $(MAINLIBS) $(SQLLIBDEP)
		@$(MAKELIB) $(PLATFORM) $(LD) "$(LDFLAGS)" \
		   "$(SOFLAGS)" libSQLIO.$(SOEXT) $@ "$(SQLO) $(SQLDO)" \
		   "$(SQLLIBEXTRA)"

$(call pcmrule,SQL)
	$(noop)

$(SQLDS):       $(SQLH) $(SQLL) $(ROOTCLINGEXE) $(call pcmdep,SQL)
		$(MAKEDIR)
		@echo "Generating dictionary $@..."
		$(ROOTCLINGSTAGE2) -f $@ $(call dictModule,SQL) -c $(SQLH) $(SQLL)

$(SQLMAP):      $(SQLH) $(SQLL) $(ROOTCLINGEXE) $(call pcmdep,SQL)
		$(MAKEDIR)
		@echo "Generating rootmap $@..."
		$(ROOTCLINGSTAGE2) -r $(SQLDS) $(call dictModule,SQL) -c $(SQLH) $(SQLL)

all-$(MODNAME): $(SQLLIB)

clean-$(MODNAME):
		@rm -f $(SQLO) $(SQLDO)

clean::         clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)
		@rm -f $(SQLDEP) $(SQLDS) $(SQLDH) $(SQLLIB) $(SQLMAP)

distclean::     distclean-$(MODNAME)
