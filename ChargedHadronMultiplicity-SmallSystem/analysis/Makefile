CFLAGS = -g3 -pthread -m64 -Wno-deprecated -I src/ 
CFLAGS += $(shell root-config --cflags --libs)
CFLAGS += -lEG

all : TrackletAna INTTVtxZ TrackletAna_Combinatoric

TrackletAna : src/TrackletAna.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/TrackletAna.cxx src/Hit.h src/Tracklet.h -o TrackletAna

INTTVtxZ : src/INTTVtxZ.cxx
	g++ $(CFLAGS) src/INTTVtxZ.cxx -o INTTVtxZ

TrackletAna_Combinatoric : src/TrackletAna_Combinatoric.cxx src/Hit.h src/Tracklet.h src/Tracklet_combinatoric.h
	g++ ${CFLAGS} src/TrackletAna_Combinatoric.cxx src/Hit.h src/Tracklet.h src/Tracklet_combinatoric.h -o TrackletAna_Combinatoric

clean :
	rm TrackletAna INTTVtxZ TrackletAna_Combinatoric