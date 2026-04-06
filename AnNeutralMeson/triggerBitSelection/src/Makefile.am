AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libtriggerbitselection.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

pkginclude_HEADERS = \
  DstPhotonSelection.h \
  DstMBDSelection.h \
  DstPhotonMBDSelection.h

pcmdir = $(libdir)
nobase_dist_pcm_DATA = 

libtriggerbitselection_la_SOURCES = \
  $(ROOTDICTS) \
  DstPhotonSelection.cc \
  DstMBDSelection.cc \
  DstPhotonMBDSelection.cc

libtriggerbitselection_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lffarawobjects \
  -lfun4all \
  -luspin \
  -lodbc++

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libtriggerbitselection.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f $(BUILT_SOURCES)
