#include <string>
#include <sstream>
#include <iostream>
#include <format>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <ffamodules/CDBInterface.h>
#include "Calo_CalibLocal.C"

#include <caloreco/RawClusterBuilderTopo.h>

#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetCalib.h>
#include <jetbackground/FastJetAlgoSub.h>

#include <g4jets/TruthJetInput.h>
//#include <globalvertex/GlobalVertexReco.h>
//#include <GlobalVertex.h>                                                            
//#include <MbdDigitization.h>
//#include <MbdReco.h>

#include <jetbackground/TimingCut.h>

#include <dummy/dummy.h>
#include <jetbackground/RetowerCEMC.h>
#include <fstream>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4dst.so)

void Fun4All_VandySkimmerTruthCounter(const std::string caloDSTlist, const std::string jetDSTlist, const std::string g4HitsDSTlist, const std::string globalDSTlist, const std::string outDir = "/sphenix/tg/tg01/jets/bkimelman/wEEC/", const std::string nfiles="25")
{

  bool doSim = true;
	
  Fun4AllServer* se=Fun4AllServer::instance();
  se->Verbosity(0);


  int runnumber = 0;
  int seg = 0;
  int n=0;
  n = std::stoi(nfiles);

  std::ifstream ifs(caloDSTlist);
  std::string filepath;
  std::getline(ifs,filepath);
  std::pair<int,int> runseg = Fun4AllUtils::GetRunSegment(filepath);
  runnumber = runseg.first;
  seg = runseg.second;

  std::string sample_name {"Jet20"};
  std::stringstream fullfilename (caloDSTlist);
  std::string temp1, temp2;
  while(std::getline(fullfilename, temp1, '/'))
  {
	  if(temp1.find("Jet") == std::string::npos) continue;
	  std::stringstream filetag (temp1);
	  while(std::getline(filetag, temp2, '_'))
	  {
		  if(temp2.find("data") == std::string::npos){
			  sample_name = temp2;
			  break;
		  }
	  }
	  break;
  }
  Fun4AllInputManager *inCalo = new Fun4AllDstInputManager("InputManagerCalo");
  inCalo->AddListFile(caloDSTlist);
  se->registerInputManager(inCalo);
  
  Fun4AllInputManager *inJet = new Fun4AllDstInputManager("InputManagerJet");
  inJet->AddListFile(jetDSTlist);
  se->registerInputManager(inJet);
  
  Fun4AllInputManager *inTruth = new Fun4AllDstInputManager("InputManagerG4Hits");
  inTruth->AddListFile(g4HitsDSTlist);
  se->registerInputManager(inTruth);
  
  Fun4AllInputManager *inGlobal = new Fun4AllDstInputManager("InputManagerGlobal");
  inGlobal->AddListFile(globalDSTlist);
  se->registerInputManager(inGlobal);
  
  dummy* dm = new dummy(); 
  se->registerSubsystem(dm);
  se->run(12);
  se->End();
  int n_events=dm->n_evts;
  std::string end_text=std::format("~/{}_starting_seg_{}.txt", sample_name, seg);
  std::fstream outfile (end_text);
  outfile<<n_events<<std::endl;
  std::cout<<dm->n_evts<<std::endl;
  outfile.close();
 // se->PrintTimer();
//  delete se;
//  gSystem->Exit(0);
}
