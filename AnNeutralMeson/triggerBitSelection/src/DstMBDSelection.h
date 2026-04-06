#ifndef __DST_MBD_SELECTION_H__
#define __DST_MBD_SELECTION_H__

#include <fun4all/SubsysReco.h>
#include <ffarawobjects/Gl1Packet.h>

#include <cstdint> // for uint64_t

// Forward declarations
class PHCompositeNode;

class DstMBDSelection : public SubsysReco
{
 public:
  //! constructor
  DstMBDSelection(const std::string &name = "DstMBDSelection",
                  const int runnb = 48746);

  //! destructor
  virtual ~DstMBDSelection();

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
  uint64_t scaled_trigger;
  
  // trigger selection
  bool trigger_mbd_any_vtx = false;
  bool trigger_mbd_vtx_10 = false;

  // Global Level 1 Trigger Packet
  Gl1Packet *gl1packet = nullptr; 
};

#endif
