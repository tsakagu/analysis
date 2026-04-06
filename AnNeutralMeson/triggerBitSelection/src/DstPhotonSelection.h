#ifndef __DST_PHOTON_SELECTION_H__
#define __DST_PHOTON_SELECTION_H__

#include <fun4all/SubsysReco.h>
#include <ffarawobjects/Gl1Packet.h>

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <map>

// Forward declarations
class PHCompositeNode;

class DstPhotonSelection : public SubsysReco
{
 public:
  //! constructor
  DstPhotonSelection(const std::string &name = "DstPhotonSelection",
                                  const int runnb = 48746);

  //! destructor
  virtual ~DstPhotonSelection();

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
  bool triggerSelection = false;
  bool trigger_mbd_photon_3 = false;
  bool trigger_mbd_photon_4 = false;
  bool trigger_mbd_photon_5 = false;

  // Global Level 1 Trigger Info
  Gl1Packet *gl1packet = nullptr;
  
  // Scaledown list:
  int scaledown[64] = {0};
  
  // For the SQL access
  std::string db_name = "daq";
  std::string user_name = "phnxrc";
  std::string table_name = "gl1_scaledown";
};

#endif
