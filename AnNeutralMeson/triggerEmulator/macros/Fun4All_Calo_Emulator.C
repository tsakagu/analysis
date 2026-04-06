#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
//#include <calotrigger/CaloTriggerEmulator.h>
#include <triggeremulator/CaloTriggerEmulatorAnNeutral.h>
#include <phool/recoConsts.h>

#include <regex>

R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libtriggeremulator.so)
#endif

void Fun4All_Calo_Emulator(int nEvents = 120,
                           const std::string& dst_event_filename = "DST_PHOTON_TRIGGER_LIST_run2pp_ana502_nocdbtag_v001-00053686-00000.root",
                           const std::string& dst_triggered_event_filename = "DST_TRIGGERED_EVENT_seb00_run2pp_ana502_nocdbtag_v001-00053686-00000.root",
                           const std::string& dst_outname = "DST_TRIGGER_EMULATOR_run2pp_ana502_nocdbtag_v001-00053686-00000.root")
{ 
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(dst_event_filename);
  int runnumber = runseg.first;

  std::cout << "runnumber = " << runnumber << std::endl;
  
  Fun4AllServer *se = Fun4AllServer::instance();

  recoConsts *rc = recoConsts::instance();

  //===============
  // conditions DB flags
  //===============

  //CDB to grab Lookup tables
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  rc->set_IntFlag("RUNNUMBER", runnumber);
  
  int verbosity = 0;
  se->Verbosity(verbosity);

  // Read micro-DSTs
  // To treat fewer events
  Fun4AllDstInputManager *in = nullptr;
  in = new Fun4AllDstInputManager("in_micro");
  in->AddFile(dst_event_filename);
  se->registerInputManager(in);

  // Read Calo Packets
  std::string dst_file;
  std::string inputname = "";
  std::string seb_string = "seb00";
  std::stringstream seb_stringstream;
  for (int iseb = 0; iseb < 16; iseb++) // i.e. all EMCal sebs
  {
    seb_stringstream.str("");
    seb_stringstream << "seb" << std::setfill('0') << std::setw(2) << iseb;
    seb_string = seb_stringstream.str();
    dst_file =std::regex_replace(dst_triggered_event_filename, std::regex("seb00"), seb_string);
    inputname = "in_" + std::to_string(iseb);
    in = new Fun4AllDstInputManager(inputname);
    in->AddFile(dst_file);
    se->registerInputManager(in);
  }

  CaloTriggerEmulatorAnNeutral *te = new CaloTriggerEmulatorAnNeutral("CALOTRIGGEREMULATOR");
  te->Verbosity(10);
  te->setTriggerType("PHOTON"); // photon trigger
  te->setNSamples(12);
  te->setTriggerSample(6);
  //te->setJetThreshold(8, 13, 17, 22);
  te->setPhotonThreshold(8, 13, 17, 22); // Actually 8, 13, 17, 22
  // subtraction delay of the post and pre sample
  te->setTriggerDelay(5);
  te->SetIsData(true);
  te->SetUseOffline(true);
  te->Verbosity(0);
  std::stringstream OptMaskFilename;
  OptMaskFilename << "/sphenix/user/dlis/Projects/macros/CDBTest/ll1config/optmask_" << runnumber << ".root";
  std::stringstream EmcalLUTFilename;
  EmcalLUTFilename << "/sphenix/user/dlis/Projects/macros/CDBTest/adc_lut_root/emcal_luts_run" << runnumber << ".root";
  te->setOptMaskFile(OptMaskFilename.str().c_str());
  te->setEmcalLUTFile(EmcalLUTFilename.str().c_str());
  //te->setEmcalLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/emcal_ll1_lut_0.50tr_new.root");
  //te->setHcalinLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/hcalin_ll1_lut_0.50tr_new.root");
  //te->setHcaloutLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/hcalout_ll1_lut_0.50tr_new.root");
  se->registerSubsystem(te);

  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("TriggerOut", dst_outname);
  out->SplitLevel(99);
  // Only keep TRIGGERTILE in output DST
  out->StripNode("14001"); // GL1 Packet
  for (int ipacket = 6001; ipacket <= 6128; ipacket++)
  {
    out->StripNode(std::to_string(ipacket));
  }
  out->StripNode("EventHeader");
  out->StripNode("LL1OUT_PHOTON");
  out->StripNode("LL1OUT_JET");
  out->StripNode("LL1OUT_PAIR");
  out->StripNode("TRIGGERPRIMITIVES_PHOTON");
  out->StripNode("TRIGGERPRIMITIVES_JET");
  out->StripNode("TRIGGERPRIMITIVES_PAIR");
  out->StripNode("TRIGGERPRIMITIVES_EMCAL");
  out->StripNode("TRIGGERPRIMITIVES_EMCAL_LL1");
  out->StripNode("TRIGGERPRIMITIVES_EMCAL_2x2_LL1");
  out->StripNode("TRIGGERPRIMITIVES_HCALOUT");
  out->StripNode("TRIGGERPRIMITIVES_HCAL_LL1");
  out->StripNode("TRIGGERPRIMITIVES_HCALIN");
  out->StripRunNode("RunHeader");
  out->StripRunNode("Flags");
  out->StripRunNode("CdbUrl");
  out->StripRunNode("TriggerRunInfo");
  
  se->registerOutputManager(out);
  

  // Fun4All
  se->Print("NODETREE");
  se->run(nEvents);
  se->End();
  se->PrintTimer();
  gSystem->Exit(0);
}
