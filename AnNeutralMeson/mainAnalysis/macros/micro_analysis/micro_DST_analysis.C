#include <anneutral/AnNeutralMeson_micro.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libanneutral.so)

void micro_DST_analysis(const int runnumber = 47289)
{
  std::string inputfilename = "/sphenix/tg/tg01/coldqcd/vmahaut/AnNeutralMeson/combine_out_run_ana450_2024p009/OUTDIPHOTON_" + std::to_string(runnumber) + ".bin";
  std::string outputfilename = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/analysis_per_run_allruns/analysis_" + std::to_string(runnumber) + ".root";
  std::string outputfiletreename = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/analysis_per_run_allruns/diphoton_minimal_" + std::to_string(runnumber) + ".root";
  //std::string inputfilename = "../combine_out_run_48746_05_1000/OUTDIPHOTON_" + std::to_string(runnumber) + ".root";
  //std::string outputfilename = "analysis_per_run_48746_asymmetry/analysis_" + std::to_string(runnumber) + ".root";
  
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1000);

  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  AnNeutralMeson_micro *AN = new AnNeutralMeson_micro("AnNeutralMeson_micro", runnumber, inputfilename, outputfilename,outputfiletreename);
  std::vector<float> chi2_cuts = {1000}; // {1000, 6, 4, 2};
  //std::vector<float> chi2_cuts = {1000, 6, 4, 2};
  std::vector<float> ecore_cuts = {1.0}; // {0.5, 1.0, 1.5};
  //std::vector<float> ecore_cuts = {0.5, 1.0, 1.5};
  float alpha_cut = 0.7;
  float pT_cut = 1.0;
  AN->set_chi2cut(chi2_cuts);
  AN->set_ecorecut(ecore_cuts);
  AN->set_alphacut(alpha_cut);
  AN->set_ptcut(pT_cut);
  AN->set_triggerselection(true);
  se->registerSubsystem(AN);

  se->run(1);
  se->End();
  delete se;
  
  gSystem->Exit(0);
}
