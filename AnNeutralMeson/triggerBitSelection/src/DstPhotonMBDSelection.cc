#include "DstPhotonMBDSelection.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

// Spin DB
#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <uspin/SpinDBContent.h>
#include <uspin/SpinDBContentv1.h>
#include <uspin/SpinDBOutput.h>

#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>


DstPhotonMBDSelection::DstPhotonMBDSelection(const std::string &name, const int runnb)
  : SubsysReco(name)
  , runnumber(runnb)
{
}

DstPhotonMBDSelection::~DstPhotonMBDSelection()
{
}

int DstPhotonMBDSelection::Init(PHCompositeNode *)
{  
  return 0;
}

int DstPhotonMBDSelection::InitRun(PHCompositeNode *topNode)
{
  gl1packet = findNode::getClass<Gl1Packet>(topNode, "14001");
  if (!gl1packet)
  {
    std::cerr << "AnNeutralMeson Gl1Packet node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // See how the trigger logic is defined for this run in the CDB.
  odbc::Connection *con = nullptr;
  try
  {
    con = odbc::DriverManager::getConnection(db_name, user_name, "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::stringstream cmd;
  cmd << "SELECT * FROM " << table_name << " WHERE runnumber=" << runnumber << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error: Can't find GL1 scaledown data for run {}", runnumber) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int ncols = rs->getMetaData()->getColumnCount();
  for (int icol = 0; icol < ncols; icol++)
  {
    std::string col_name = rs->getMetaData()->getColumnName(icol + 1);
    // i.e. scaledownXX where XX is the trigger number (00 to 63
    if (col_name.starts_with("scaledown")) 
    {
      int scaledown_index = std::stoi(col_name.substr(9,2));
      std::string cont = rs->getString(col_name);
      scaledown[scaledown_index] = std::stoi(cont);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int DstPhotonMBDSelection::process_event(PHCompositeNode *)
{
  // Store event-level entries
  live_trigger = gl1packet->getLiveVector();
  scaled_trigger = gl1packet->getScaledVector();

  trigger_mbd_photon_3 = false;
  trigger_mbd_photon_4 = false;
  trigger_mbd_photon_5 = false;

  // If (photon 3 GeV + MBD NS >= 1) is enabled for recording
  if (scaledown[25] != -1) 
  {
    trigger_mbd_photon_3 = (trigger_mbd_photon_3 ||
                            ((scaled_trigger >> 25 & 0x1U) == 0x1U));
  }
  // Else if (photon 3 GeV + MBD NS >=1, vtx < 10) is enabled for recording
  if (!trigger_mbd_photon_3 && (scaledown[36] != -1))
  {
    trigger_mbd_photon_3 = (trigger_mbd_photon_3 ||
                            ((scaled_trigger >> 36 & 0x1U) == 0x1U));
  }
  // Else if (photon 3 GeV) is enabled for recording with no prescale
  if (!trigger_mbd_photon_3 && (scaledown[29] == 0))
  {
    trigger_mbd_photon_3 = (trigger_mbd_photon_3 ||
                            ((live_trigger >> 25 & 0x1U) == 0x1U) ||
                            ((live_trigger >> 36 & 0x1U) == 0x1U));
  }

  // If (photon 4 GeV + MBD NS >= 1) is enabled for recording
  if (scaledown[26] != -1) 
  {
    trigger_mbd_photon_4 = (trigger_mbd_photon_4 ||
                            ((scaled_trigger >> 26 & 0x1U) == 0x1U));
  }
  // Else if (photon 4 GeV + MBD NS >=1, vtx < 10) is enabled for recording
  if (!trigger_mbd_photon_4 && (scaledown[37] != -1))
  {
    trigger_mbd_photon_4 = (trigger_mbd_photon_4 ||
                            ((scaled_trigger >> 37 & 0x1U) == 0x1U));
  }
  // Else if (photon 4 GeV) is enabled for recording with no prescale
  if (!trigger_mbd_photon_4 && (scaledown[30] == 0))
  {
    trigger_mbd_photon_4 = (trigger_mbd_photon_4 ||
                            ((live_trigger >> 26 & 0x1U) == 0x1U) ||
                            ((live_trigger >> 37 & 0x1U) == 0x1U));
  }

  // If (photon 5 GeV + MBD NS >= 1) is enabled for recording
  if (scaledown[27] != -1) 
  {
    trigger_mbd_photon_5 = (trigger_mbd_photon_5 ||
                            ((scaled_trigger >> 27 & 0x1U) == 0x1U));
  }
  // Else if (photon 5 GeV + MBD NS >=1, vtx < 10) is enabled for recording
  if (!trigger_mbd_photon_5 && (scaledown[38] != -1))
  {
    trigger_mbd_photon_5 = (trigger_mbd_photon_5 ||
                            ((scaled_trigger >> 38 & 0x1U) == 0x1U));
  }
  // Else if (photon 5 GeV) is enabled for recording with no prescale
  if (!trigger_mbd_photon_5 && (scaledown[31] == 0))
  {
    trigger_mbd_photon_5 = (trigger_mbd_photon_5 ||
                            ((live_trigger >> 27 & 0x1U) == 0x1U) || 
                            ((live_trigger >> 38 & 0x1U) == 0x1U));
  }

  // MBD only selection
  trigger_mbd_any_vtx = ((scaled_trigger >> 10 & 0x1U) == 0x1U);
  trigger_mbd_vtx_10 = ((scaled_trigger >> 12 & 0x1U) == 0x1U);
  
  if (!(trigger_mbd_any_vtx || trigger_mbd_vtx_10 || trigger_mbd_photon_3 || trigger_mbd_photon_4 || trigger_mbd_photon_5))
    return Fun4AllReturnCodes::ABORTEVENT;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int DstPhotonMBDSelection::End(PHCompositeNode *)
{
  return 0;
}
