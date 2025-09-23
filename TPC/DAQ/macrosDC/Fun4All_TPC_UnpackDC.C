
#ifndef MACRO_FUN4ALL_READDST_C
#define MACRO_FUN4ALL_READDST_C

#include <tpcrawtottree/TPCRawDataDCChecker.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libTPCRawDataDCChecker.so)

//int Fun4All_TPC_UnpackPRDF(const int nEvents = 100,
int Fun4All_TPC_UnpackDC(const int nEvents = 10,
//                           const string &inputFile = "/sphenix/lustre01/sphnxpro/rawdata/commissioning/tpc/beam/TPC_ebdc*_beam-00011012-0000.prdf"  //
                           const string &inputFile = "/bbox/bbox0/tpc/calib/TPC_ebdc00_calib-00045887-0000.evt"

)
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
//  RecoConst rc;
//  rc->SetIntFlag("RUNNUMBER",2);

  string outDir = "/sphenix/u/jamesj3j3/tpc/sPHENIXProjects/beam-run-11012/";

  string fileName = inputFile;
  size_t pos = fileName.find("TPC_ebdc");
  fileName.erase(fileName.begin(),fileName.begin()+pos);
  
  //TPCRawDataTree *r2tree = new TPCRawDataTree(outDir + fileName + "_TPCRawDataTree_skip100.root");/////////////////////////////
  TPCRawDataDCChecker *r2tree = new TPCRawDataDCChecker( fileName + "_TPCRawDataTree_skip10.root");/////////////////////////////
//  TPCRawDataTree *r2tree = new TPCRawDataTree( fileName + "_TPCRawDataTree_skip100.root");/////////////////////////////

  // add all possible TPC packet that we need to analyze
  for (int packet = 4000; packet<=4230; packet+=10)
  {
    r2tree->AddPacket(packet);
    r2tree->AddPacket(packet+1);
  }
//  r2tree->AddPacket(6000);
//  r2tree->AddPacket(6001);
  r2tree->AddPacket(5000);
  r2tree->AddPacket(5001);

  se->registerSubsystem(r2tree);

  Fun4AllPrdfInputManager *in1 = new Fun4AllPrdfInputManager("PRDF1");
  in1->AddFile(inputFile);
  se->registerInputManager(in1);

  se->skip(10);
  se->run(nEvents);

  se->End();

  delete se;
  std::cout << "All done processing" << std::endl;
  gSystem->Exit(0);
  return 0;
}
#endif
