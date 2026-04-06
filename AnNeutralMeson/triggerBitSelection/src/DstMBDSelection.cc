#include "DstMBDSelection.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

DstMBDSelection::DstMBDSelection(const std::string &name, const int runnb)
  : SubsysReco(name)
  , runnumber(runnb)
{
}

DstMBDSelection::~DstMBDSelection()
{
}

int DstMBDSelection::Init(PHCompositeNode *)
{  
  return 0;
}

int DstMBDSelection::InitRun(PHCompositeNode *topNode)
{
  gl1packet = findNode::getClass<Gl1Packet>(topNode, "14001");
  if (!gl1packet)
  {
    std::cerr << "AnNeutralMeson Gl1Packet node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int DstMBDSelection::process_event(PHCompositeNode *)
{
  // Store event-level entries
  scaled_trigger = gl1packet->getScaledVector();

  trigger_mbd_any_vtx = ((scaled_trigger >> 10 & 0x1U) == 0x1U);
  trigger_mbd_vtx_10 = ((scaled_trigger >> 12 & 0x1U) == 0x1U);
  
  if (!(trigger_mbd_any_vtx || trigger_mbd_vtx_10))
    return Fun4AllReturnCodes::ABORTEVENT;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int DstMBDSelection::End(PHCompositeNode *)
{
  return 0;
}
