#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>

#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>

#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloTowerDefs.h>
#include <caloreco/RawClusterBuilderTemplate.h>

#include <phool/recoConsts.h>

#include <ffamodules/CDBInterface.h>
#include <GlobalVariables.C>
#include <Calo_Calib.C>

#include <ctime>
#include <iostream>

R__LOAD_LIBRARY(libfun4all.so)

#include <anneutral/AnNeutralMeson.h>
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libanneutral.so)

void Fun4All_AnNeutralMeson(int nevents = 1,
                            const std::string &inputfile_events = "DST_PHOTON_MBD_TRIGGER_LIST_run2pp_ana502_nocdbtag_v001-00053686-00000.root",
                            const std::string &inputlist1 = "dst_calo.list",
                            const std::string &inputlist2 = "dst_calofitting.list",
                            const std::string &histname = "hist_suffix", const std::string &dstname = "DST_MICRO_ANNEUTRALMESON_1")
{

  // Print the date (important to know what were the calibrations used)
  std::time_t t = std::time(0); // get time now
  std::tm* now = std::localtime(&t);
  std::cout << "date = " << (now->tm_year + 1900) << "-" << (now->tm_mon + 1) << "-" << now->tm_mday << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  cout << "registered F4A" << std::endl;

  // se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  //===============
  // conditions DB flags
  //===============
  ifstream file(inputlist2);
  string first_file;
  getline(file, first_file);
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  cout << "run number = " << runnumber << endl;

  // global tag
  // cf gmattson, apply the tower mask from the CDB
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  rc->set_IntFlag("RUNNUMBER", runnumber);


  std::cout << "Registering selection DST" << std::endl;
  // Only run the code on triggered events (defined in DST_PHOTON_MBD_TRIGGER_LIST)
  // The event DST only contains DST#Sync node
  Fun4AllInputManager *in_events = new Fun4AllDstInputManager("DST_SELECTED_EVENTS");
  in_events->AddFile(inputfile_events);
  se->registerInputManager(in_events);

  if (inputlist1 != "" && inputlist2 == "")
  {
    std::cout << "Registering DST_CALO" << std::endl;
    Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_CALO_CLUSTERS");
    in->AddListFile(inputlist1);
    se->registerInputManager(in);
  }

  // If the calofitting list is given, use Process_Calo_Calib()

  if (inputlist2 != "")
  {
    std::cout << "Registering DST_CALOFITTING" << std::endl;
    Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DST_CALOFITTINGS");
    in2->AddListFile(inputlist2);
    se->registerInputManager(in2);
  }

  // Tower calibrations & masking
  if (inputlist2 != "") {
    // MBD/BBC Reconstruction
    MbdReco *mbdreco = new MbdReco();
    se->registerSubsystem(mbdreco);

    // Official vertex storage
    GlobalVertexReco *gvertex = new GlobalVertexReco();
    se->registerSubsystem(gvertex);
    
    std::cout << "Using Process_Calo_Calib()" << std::endl;
    Process_Calo_Calib();  // this line handles the calibrations, dead/hot tower masking and reruns the clusterizer
  }

  ///////////////////
  // analysis modules
  // foldername = ".";
  std::cout << "Registering AnNeutralMeson" << std::endl;
  AnNeutralMeson *AN = new AnNeutralMeson("calomodulename", runnumber, histname, ".");
  se->registerSubsystem(AN);

  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", dstname);
  out->SplitLevel(99);
  // Keep minimal information in the output node: Event Number + Global Vertex + Small Cluster Info
  out->StripNode("EventHeader");
  out->StripNode("14001");
  out->StripNode("1002");
  out->StripNode("1001");
  out->StripNode("TOWERS_ZDC");
  out->StripNode("TOWERS_CEMC");
  out->StripNode("TOWERINFO_CALIB_CEMC");
  out->StripNode("CLUSTERINFO_CEMC");
  // out->StripNode("CLUSTER_SMALLINFO_CEMC");
  out->StripNode("TOWERS_HCALIN");
  out->StripNode("TOWERINFO_CALIB_HCALIN");
  out->StripNode("TOWERS_HCALOUT");
  out->StripNode("TOWERINFO_CALIB_HCALOUT");
  out->StripNode("TOWERS_SEPD");
  out->StripNode("MbdOut");
  out->StripNode("MbdRawContainer");
  out->StripNode("MbdPmtContainer");
  out->StripNode("MbdVertexMap");
  // out->StripNode("GlobalVertexMap");
  out->StripRunNode("TriggerRunInfo");
  out->StripRunNode("RunHeader");
  out->StripRunNode("Flags");
  out->StripRunNode("CdbUrl");
  out->StripRunNode("TOWERGEOM_CEMC_DETAILED");
  out->StripRunNode("TOWERGEOM_CEMC");
  out->StripRunNode("TOWERGEOM_HCALIN");
  out->StripRunNode("TOWERGEOM_HCALOUT");
  out->StripRunNode("CYLINDERGEOM_CEMC");
  out->StripRunNode("CYLINDERCELLGEOM_CEMC");
  out->StripRunNode("MbdGeom");
  se->registerOutputManager(out);

  std::cout << "Running F4A" << std::endl;
  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;

  TFile* f_done_signal = new TFile("DONE.root","recreate");
  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);
}

