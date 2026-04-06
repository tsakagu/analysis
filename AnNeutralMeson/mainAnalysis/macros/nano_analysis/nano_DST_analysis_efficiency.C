#include <anneutral/AnNeutralMeson_nano.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libanneutral.so)

void nano_DST_analysis_efficiency(const int seednumber = 0)
{
  std::string inputlistname = "inputruns.txt";
  std::string inputfiletemplate = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/trees_complete_ana509_01312026_efficiency_30/diphoton_minimal_";
  std::string outputfiletemplate = "analysis_complete_ana509_01312026_backward_efficiency_30_pt_g3_0mrad/analysis_";
  
  Fun4AllServer *se = Fun4AllServer::instance();
  //se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024"); // Required for the spin DB
  rc->set_uint64Flag("TIMESTAMP", 48746); // I am not using the TIMESTAMP anyway.

  AnNeutralMeson_nano *AN = new AnNeutralMeson_nano("AnNeutralMeson_nano", inputlistname, inputfiletemplate, outputfiletemplate);
  AN->set_store_tree(true);
  //AN->set_forward_cut(true);
  AN->set_backward_cut(true);
  AN->set_ptcut(3.0, 1000.0);
  se->registerSubsystem(AN);
  se->run(1);
  se->End();
  delete se;
  
  gSystem->Exit(0);
}
