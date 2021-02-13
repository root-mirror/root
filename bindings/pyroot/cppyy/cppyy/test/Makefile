dicts = advancedcppDict.so \
        advancedcpp2Dict.so \
        conversionsDict.so \
        cpp11featuresDict.so \
        crossinheritanceDict.so \
        datatypesDict.so \
        doc_helperDict.so \
        example01Dict.so \
        fragileDict.so \
        operatorsDict.so \
        overloadsDict.so \
        pythonizablesDict.so \
        std_streamsDict.so \
        stltypesDict.so \
        templatesDict.so
all : $(dicts)

genreflex_flags := $(shell genreflex --cppflags)
cppflags=$(shell cling-config --cppflags) $(genreflex_flags) -O3 -fPIC -I$(shell python -c 'import distutils.sysconfig as ds; print(ds.get_config_var("INCLUDEPY"))') -Wno-register

PLATFORM := $(shell uname -s)
ifeq ($(PLATFORM),Darwin)
  cppflags+=-dynamiclib -single_module -arch x86_64 -undefined dynamic_lookup -Wno-delete-non-virtual-dtor
endif

%Dict.so: %_rflx.cpp %.cxx
	$(CXX) $(cppflags) -shared -o $@ $^

%_rflx.cpp: %.h %.xml
	genreflex $< --selection=$*.xml --rootmap=$*Dict.rootmap --rootmap-lib=$*Dict.so

.PHONY: test clean

test:
	pytest test_*.py

clean:
	-rm -f $(dicts) $(subst .so,.rootmap,$(dicts)) $(subst Dict.so,_rflx_rdict.pcm,$(dicts)) $(subst Dict.so,_rflx.cpp,$(dicts)) $(wildcard *.pyc)
