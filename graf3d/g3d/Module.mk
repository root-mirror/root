# Module.mk for g3d module
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

MODNAME      := g3d
MODDIR       := $(ROOT_SRCDIR)/graf3d/$(MODNAME)
MODDIRS      := $(MODDIR)/src
MODDIRI      := $(MODDIR)/inc

G3DDIR       := $(MODDIR)
G3DDIRS      := $(G3DDIR)/src
G3DDIRI      := $(G3DDIR)/inc

##### libGraf3d #####
G3DL         := $(MODDIRI)/LinkDef.h
G3DDS        := $(call stripsrc,$(MODDIRS)/G__Graf3d.cxx)
G3DDO        := $(G3DDS:.cxx=.o)
G3DDH        := $(G3DDS:.cxx=.h)

G3DH1        := $(filter-out $(MODDIRI)/LinkDef%,$(wildcard $(MODDIRI)/*.h))
G3DH1        := $(filter-out $(MODDIRI)/X3DBuffer%,$(G3DH1))
G3DH1        := $(filter-out $(MODDIRI)/X3DDefs%,$(G3DH1))
G3DH2        := $(MODDIRI)/X3DBuffer.h 
G3DH         := $(G3DH1) $(G3DH2)
G3DS1        := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.cxx))
G3DS2        := $(wildcard $(MODDIRS)/*.c)
G3DO         := $(call stripsrc,$(G3DS1:.cxx=.o) $(G3DS2:.c=.o))

G3DDEP       := $(G3DO:.o=.d) $(G3DDO:.o=.d)

G3DLIB       := $(LPATH)/libGraf3d.$(SOEXT)
G3DMAP       := $(G3DLIB:.$(SOEXT)=.rootmap)

# used in the main Makefile
G3DH1_REL    := $(patsubst $(MODDIRI)/%.h,include/%.h,$(G3DH1))
G3DH2_REL    := $(patsubst $(MODDIRI)/%.h,include/%.h,$(G3DH2))
ALLHDRS     += $(G3DH1_REL) $(G3DH2_REL)
ALLLIBS     += $(G3DLIB)
ALLMAPS     += $(G3DMAP)
ifeq ($(CXXMODULES),yes)
  CXXMODULES_HEADERS := $(patsubst include/%,header \"%\"\\n,$(G3DH1_REL))
  CXXMODULES_MODULEMAP_CONTENTS += module Graf3d_$(MODNAME) { \\n
  CXXMODULES_MODULEMAP_CONTENTS += $(CXXMODULES_HEADERS)
  CXXMODULES_MODULEMAP_CONTENTS += "export * \\n"
  CXXMODULES_MODULEMAP_CONTENTS += link \"$(G3DLIB)\" \\n
  CXXMODULES_MODULEMAP_CONTENTS += } \\n

  CXXMODULES_HEADERS := $(patsubst include/%,header \"%\"\\n,$(G3DH2_REL))
  CXXMODULES_MODULEMAP_CONTENTS += module Graf3d_$(MODNAME)_C [extern_c] { \\n
  CXXMODULES_MODULEMAP_CONTENTS += $(CXXMODULES_HEADERS)
  CXXMODULES_MODULEMAP_CONTENTS += "export * \\n"
  CXXMODULES_MODULEMAP_CONTENTS += link \"$(G3DLIB)\" \\n
  CXXMODULES_MODULEMAP_CONTENTS += } \\n
endif

# include all dependency files
INCLUDEFILES += $(G3DDEP)

##### local rules #####
.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

include/%.h:    $(G3DDIRI)/%.h
		cp $< $@

$(G3DLIB):      $(G3DO) $(G3DDO) $(ORDER_) $(MAINLIBS) $(G3DLIBDEP)
		@$(MAKELIB) $(PLATFORM) $(LD) "$(LDFLAGS)" \
		   "$(SOFLAGS)" libGraf3d.$(SOEXT) $@ "$(G3DO) $(G3DDO)" \
		   "$(G3DLIBEXTRA)"

$(call pcmrule,G3D)
	$(noop)

$(G3DDS):       $(G3DH1) $(G3DL) $(ROOTCLINGEXE) $(call pcmdep,G3D)
		$(MAKEDIR)
		@echo "Generating dictionary $@..."
		$(ROOTCLINGSTAGE2) -f $@ $(call dictModule,G3D) -c $(G3DH1) $(G3DL)

$(G3DMAP):      $(G3DH1) $(G3DL) $(ROOTCLINGEXE) $(call pcmdep,G3D)
		$(MAKEDIR)
		@echo "Generating rootmap $@..."
		$(ROOTCLINGSTAGE2) -r $(G3DDS) $(call dictModule,G3D) -c $(G3DH1) $(G3DL)

all-$(MODNAME): $(G3DLIB)

clean-$(MODNAME):
		@rm -f $(G3DO) $(G3DDO)

clean::         clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)
		@rm -f $(G3DDEP) $(G3DDS) $(G3DDH) $(G3DLIB) $(G3DMAP)

distclean::     distclean-$(MODNAME)
