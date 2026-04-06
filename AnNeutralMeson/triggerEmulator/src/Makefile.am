AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir` \
  -isystem$(OPT_SPHENIX)/include

lib_LTLIBRARIES = \
   libtriggeremulator_io.la \
   libtriggeremulator.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtriggeremulator_io_la_LIBADD = \
  -lphool \
  -lodbc++ \
  -lcalo_io \
  -lmbd_io \
  -lzdcinfo_io

libtriggeremulator_la_LIBADD = \
  libtriggeremulator_io.la \
  -lcalotrigger

pkginclude_HEADERS = \
  CaloTriggerEmulatorAnNeutral.h \
  TriggerTile.h

ROOTDICTS = \
  TriggerTile_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libtriggeremulator_io_la_SOURCES = \
  $(ROOTDICTS) \
  TriggerTile.cc

libtriggeremulator_la_SOURCES = \
  CaloTriggerEmulatorAnNeutral.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_triggeremulator_io \
  testexternals_triggeremulator

testexternals_triggeremulator_io_SOURCES = testexternals.cc
testexternals_triggeremulator_io_LDADD = libtriggeremulator_io.la

testexternals_triggeremulator_SOURCES = testexternals.cc
testexternals_triggeremulator_LDADD = libtriggeremulator.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
