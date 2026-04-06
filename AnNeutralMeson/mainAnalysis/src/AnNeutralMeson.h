#ifndef AN_NEUTRAL_MESON_H__
#define AN_NEUTRAL_MESON_H__

#include "ClusterSmallInfoContainer.h"

#include <fun4all/SubsysReco.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawClusterContainer.h>

// GL1p Trigger mapping
#include <ffarawobjects/Gl1Packet.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>

#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <vector>
#include <utility>
#include <map>
#include <fstream>
// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;

class TFile;
class TH2F;
class TH1F;
class TH1;
class TF1;


const int nbunches = 120;

class AnNeutralMeson : public SubsysReco
{
 public:

  //! constructor
  AnNeutralMeson(const std::string &name = "AnNeutralMeson", const int& runnumber = 48748, const std::string &fname = "file.root", const std::string& folder = "");

  //! destructor
  virtual ~AnNeutralMeson();

  //! full initialization
  int Init(PHCompositeNode *);

  //! Run initialization
  int InitRun(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);
  int process_towers(PHCompositeNode *);

  //! Create output DST nodes
  void CreateNodes(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

 protected:
  int _runnumber;
  std::string outfilename;
  std::string outfolder;
  TFile *outfile_histograms = nullptr;

  // Read nodes:
  Gl1Packet *gl1packet = nullptr;
  GlobalVertexMap *vertexmap = nullptr;
  RawClusterContainer *clusterContainer = nullptr;
  TowerInfoContainer *towers = nullptr;

  // event level QA histograms
  TH1F* h_event_vtx_z;

  // Cluster histograms (before any cut)
  TH1F* h_clusE;
  TH1F* h_clus_eta;
  TH1F* h_clus_phi;
  TH2F *h_clus_eta_phi;
  TH2F *h_clus_eta_E;
  TH2F *h_clus_eta_vtxz;
  TH1F* h_clus_pt;
  TH1F* h_clus_chisq;
  TH2F* h_etaphi_clus;

  // Event level cuts
  float EMCal_E_cut = 0;
  
  // cluster level cuts
  float clus_chisq_cut = 1000;
  float clus_E_cut = 0.5;

  // Counters
  int _eventcounter; // Number of events in DST
  int _eventcounter_selection1 = 0; // Nb events after MBD live trigger cut
  int _eventcounter_selection2 = 0; // Nb events after globalVtxMap cut
  int _clustercounter = 0;
  int _clustercounter_selection = 0;
  unsigned long _countentry = 0; // Number of events with at least one good cluster

  bool getVtx = false;
  
  // Event-level entries for the micro DST
  int bunchnumber = 0;
  unsigned long live_trigger = 0;
  unsigned long scaled_trigger = 0;
  float vtx_z = 0;

  // Recorded node
  std::string ClusterSmallInfoNodeName;
  ClusterSmallInfoContainer *_smallclusters;

  SyncObject *syncobject{nullptr};
};

#endif
