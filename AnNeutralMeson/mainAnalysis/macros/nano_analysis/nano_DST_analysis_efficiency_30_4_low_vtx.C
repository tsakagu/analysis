#include <anneutral/AnNeutralMeson_nano.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libanneutral.so)

void nano_DST_analysis_efficiency_30_4_low_vtx(const int seednumber = 0)
{
  std::string inputlistname = "inputruns.txt";
  std::string inputfiletemplate = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/trees_complete_ana509_01312026_efficiency_30/diphoton_minimal_";
  std::string outputfolder = "analysis_complete_ana509_01312026_low_vtx_efficiency_30_pt_g4_0mrad/";
  gSystem->Exec(("mkdir -p " + outputfolder).c_str());
  std::string outputfiletemplate = outputfolder + "/analysis_";
  std::string outputbunchtemplate = outputfolder + "/fast_analysis_";
  
  Fun4AllServer *se = Fun4AllServer::instance();

  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024"); // Required for the spin DB
  rc->set_uint64Flag("TIMESTAMP", 48746); // I am not using the TIMESTAMP anyway.

  AnNeutralMeson_nano *AN = new AnNeutralMeson_nano("AnNeutralMeson_nano", inputlistname, inputfiletemplate, outputfiletemplate);
  AN->set_store_bunch_yields(false, "");
  AN->set_trigger_photon(true);
  AN->set_ptcut(4.0, 1000.0);
  AN->set_low_vtx_cut(true);
  se->registerSubsystem(AN);
  se->run(1);
  se->End();
  delete se;
  
  gSystem->Exit(0);
}
