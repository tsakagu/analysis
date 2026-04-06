#include <anneutral/AnNeutralMeson_micro_dst.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libanneutral.so)

void micro_DST_analysis_complete(const int runnumber = 47289, const int nevents = 0)
{
  std::string inputfilename = "DST_MICRO_ANNEUTRALMESON_run2pp_ana509_2024p022_v001-000" + std::to_string(runnumber) + ".root";
  std::string outputfilename = "analysis_" + std::to_string(runnumber) + ".root";
  std::string outputfiletreename = "diphoton_minimal_" + std::to_string(runnumber) + ".root";
  
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  rc->set_IntFlag("RUNNUMBER", runnumber);

  Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  ingeo->AddFile(geoLocation);
  se->registerInputManager(ingeo);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_MICRO");
  in->fileopen(inputfilename);
  se->registerInputManager(in);

  //Fun4AllInputManager *in_emulator = new Fun4AllDstInputManager("DST_EMULATOR");
  //in_emulator->fileopen(inputemulatorname);
  //se->registerInputManager(in_emulator);
  
  AnNeutralMeson_micro_dst *AN = new AnNeutralMeson_micro_dst("AnNeutralMeson_micro_dst", runnumber, outputfilename, outputfiletreename);
  AN->set_store_tree(true);
  AN->set_store_qa(true);
  std::vector<float> chi2_cuts = {1000}; // {1000, 6, 4, 2};
  //std::vector<float> chi2_cuts = {1000, 6, 4, 2};
  std::vector<float> ecore_cuts = {1.0}; // {0.5, 1.0, 1.5};
  //std::vector<float> ecore_cuts = {0.5, 1.0, 1.5};
  float alpha_cut = 0.7;
  //float pT_cut = 1.0;
  AN->set_chi2cut(chi2_cuts);
  AN->set_ecorecut(ecore_cuts);
  AN->set_alphacut(alpha_cut);
  AN->set_ptcut(1.0, 1000.0);
  //AN->set_pt_threshold(3.0);
  AN->set_mbd_trigger_bit_requirement(true);
  AN->set_photon_trigger_bit_requirement(false);
  AN->set_photon_trigger_efficiency_matching_requirement(false);
  //AN->set_photon_trigger_efficiency_threshold(30);
  AN->set_photon_trigger_emulator_matching_requirement(false);
  AN->set_event_mixing(false);
  se->registerSubsystem(AN);

  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;
  
  gSystem->Exit(0);
}
