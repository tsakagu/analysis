#ifndef __DST_PHOTON_MBD_SELECTION_H__
#define __DST_PHOTON_MBD_SELECTION_H__

#include <fun4all/SubsysReco.h>
#include <ffarawobjects/Gl1Packet.h>

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <map>

// Forward declarations
class PHCompositeNode;
class TFile;
class TTree;
class TH1;
class TH1F;
class TH2F;
class LorentzVector;

class DstPhotonMBDSelection : public SubsysReco
{
 public:
  //! constructor
  DstPhotonMBDSelection(const std::string &name = "DstPhotonMBDSelection",
                        const int runnb = 48746);

  //! destructor
  virtual ~DstPhotonMBDSelection();

  //! full initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);
 
 private:
  
  // run number
  int runnumber;

  // Global information in clusters' node
  uint64_t live_trigger;
  uint64_t scaled_trigger;
  
  // trigger selection
  bool trigger_mbd_any_vtx = false;
  bool trigger_mbd_vtx_10 = false;
  bool triggerSelection = false;
  bool trigger_mbd_photon_3 = false;
  bool trigger_mbd_photon_4 = false;
  bool trigger_mbd_photon_5 = false;

  // Global Level 1 info
  Gl1Packet *gl1packet = nullptr;
  
  // Scaledown list:
  int scaledown[64] = {0};
  
  // For the SQL access
  std::string db_name = "daq";
  std::string user_name = "phnxrc";
  std::string table_name = "gl1_scaledown";
};

#endif
