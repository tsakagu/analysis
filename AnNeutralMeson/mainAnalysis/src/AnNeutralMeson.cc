#include "AnNeutralMeson.h"

// Spin DB
#include <uspin/SpinDBContent.h>
#include <uspin/SpinDBOutput.h>

// Tower includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/MbdVertex.h>

// MBD
#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>

#include <Math/Vector4D.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <Event/Event.h>
#include <Event/packet.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

// To check virtual and resident memory usage
//#include "monitorMemory.h"

#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <utility>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TRandom3.h"

#include <cdbobjects/CDBHistos.h>  // for CDBHistos
#include <cdbobjects/CDBTTree.h>   // for CDBTTree

#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <g4main/PHG4TruthInfoContainer.h>

using namespace std;

AnNeutralMeson::AnNeutralMeson(const std::string &name, const int &runnumber,
                               const std::string &filename, const std::string& foldername)
  : SubsysReco(name)
  , _runnumber(runnumber)
  , outfilename(filename)
  , outfolder(foldername)
{
  _eventcounter = 0;
}

AnNeutralMeson::~AnNeutralMeson()
{
  // delete hm;
  delete outfile_histograms;
}

int AnNeutralMeson::Init(PHCompositeNode *)
{
  // create and register your histos (all types) here
  outfile_histograms =
      new TFile((outfolder + "/OUTHIST_" + outfilename).c_str(),
                "RECREATE");

  // Event QA
  h_event_vtx_z =
      new TH1F(
          "h_event_vtx_z",
          ";counts; vtx_z [cm]", 200, -200, 200);

  // cluster QA
  h_clusE = new TH1F(
      "h_clusE",
      "; E [GeV]; counts", 300, 0, 15);
  h_clus_eta = new TH1F(
      "h_clus_eta",
      ";eta;counts", 200, -1.2, 1.2);
  h_clus_phi =
      new TH1F(
          "h_clus_phi",
          ";#phi; counts", 200, -TMath::Pi(), TMath::Pi());
  h_clus_eta_phi = new TH2F(
      "h_clus_eta_phi",
      ";#eta;#phi;counts", 200, -2, 2,
      160, -TMath::Pi(), TMath::Pi());
  h_clus_eta_E =
      new TH2F(
          "h_clus_eta_E",
          ";#eta;E;counts", 200, -2, 2, 150, 0, 15);
  h_clus_eta_vtxz = new TH2F(
      "h_clus_eta_vtxz",
      ";#eta;vertex z [cm];counts",
      200, -2, 2, 200, -200, 200);
  h_clus_pt = new TH1F(
      "h_clus_pt",
      ";p_{T};counts", 100, 0, 10);
  h_clus_chisq = new TH1F(
      "h_clus_chisq",
      ";#chi^{2};counts", 100, 0, 10);

  return 0;
}

int AnNeutralMeson::InitRun(PHCompositeNode *topNode)
{
  syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  if (!syncobject)
  {
    // Not used here but sync object should be present to keep the event number
    // in the output µDST
    std::cout << "syncobject not found" << std::endl;
  }
  
  gl1packet = findNode::getClass<Gl1Packet>(topNode, "14001");
  if (!gl1packet)
  {
    std::cerr << "AnNeutralMeson Gl1Packet node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  vertexmap =
    findNode::getClass<GlobalVertexMap>(topNode,
                                        "GlobalVertexMap");
  if (!vertexmap)
  {
      std::cout << "AnNeutralMeson GlobalVertexMap node is missing"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  clusterContainer =
    findNode::getClass<RawClusterContainer>(topNode,
                                            "CLUSTERINFO_CEMC");
  if (!clusterContainer)
    {
      std::cout << PHWHERE << "AnNeutralMeson - Fatal Error - "
        "CLUSTER_CEMC node is missing. "
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }
  
  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }
  
  return 0;
}

int AnNeutralMeson::process_event(PHCompositeNode *topNode)
{
  _eventcounter++;

  return process_towers(topNode);

  return 0;
}

int AnNeutralMeson::process_towers(PHCompositeNode */*topNode*/)
{
  // Show progression
  if ((_eventcounter % 1000) == 0)
    std::cout << _eventcounter << std::endl;

  //-----------------------check MBD
  // trigger----------------------------------------//

  // see Joe Mead's GL1/GTM manual (Firmware version 51)
  // Live counter increments only when the "busy" signal is inactive
  live_trigger = gl1packet->getLiveVector();
  // Scaledown counter only increments when the number of triggers reaches the
  // scaledown vector (scaledown is the ratio between live and scaled values)
  scaled_trigger = gl1packet->getScaledVector();

  // 0x1U is 1 in hexadecimal notation. U means it is unsigned
  if ((live_trigger >> 10 & 0x1U) != 0x1U)
  {
    // Ensure that the 10th bit is equal to 1 (MBD N&S >= 1)
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _eventcounter_selection1++;

  //-----------------------spin
  // information----------------------------------------//

  // Only keep the bunchnumber at this point.
  // The spin direction and average polarization can be
  // deduced later at the run level (or at the fill level).
  bunchnumber = gl1packet->getBunchNumber();

  //-----------------------get vertex----------------------------------------//

  vtx_z = 0;
    
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      vtx_z = vtx->get_z();
      h_event_vtx_z->Fill(vtx_z);
    }
    else
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  _eventcounter_selection2++;

  // Make sure ClusterSmallInfo container is empty before filling it with new clusters
  _smallclusters->Reset();
  _smallclusters->set_live_trigger(live_trigger);
  _smallclusters->set_scaled_trigger(scaled_trigger);
  _smallclusters->set_bunchnumber(bunchnumber);

  // Select good clusters
  int num_good_clusters = 0;
  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second;
       clusterIter++)
  {
    RawCluster *recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster =
        RawClusterUtility::GetECoreVec(*recoCluster, vertex);
    CLHEP::Hep3Vector pos_vec_cluster = recoCluster->get_position() - vertex;

    float clusE = E_vec_cluster.mag();
    float clus_eta = E_vec_cluster.pseudoRapidity();
    float clus_phi = E_vec_cluster.phi();
    float clus_pt = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();

    h_clusE->Fill(clusE);
    h_clus_eta->Fill(clus_eta);
    h_clus_phi->Fill(clus_phi);
    h_clus_eta_phi->Fill(clus_eta, clus_phi);
    h_clus_eta_E->Fill(clus_eta, clusE);
    h_clus_eta_vtxz->Fill(clus_eta, vtx_z);
    h_clus_pt->Fill(clus_pt);
    h_clus_chisq->Fill(clus_chisq);

    _clustercounter++; // Cluster counter before cuts

    if (clusE < clus_E_cut)
      continue;
    if (clus_chisq > clus_chisq_cut)
      continue;

    _clustercounter_selection++; // Cluster counter after cuts
    num_good_clusters++;

    // Don't store events with more than 50 clusters
    // Too noisy
    if (!(_smallclusters->add_cluster(clus_eta,
                                      clus_phi,
                                      clusE,
                                      recoCluster->get_energy(),
                                      clus_chisq)))
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }  // cluster loop

  // Compress the ClusterSmallInfoContainer object
  _smallclusters->compress();

  // Don't store events with no clusters
  if (num_good_clusters == 0)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Store any other events
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnNeutralMeson::End(PHCompositeNode * /*topNode*/)
{
  // Save QA histograms (before good cluster selection)
  outfile_histograms->cd();
  outfile_histograms->Write();
  outfile_histograms->Close();
  delete outfile_histograms;

  // Count events and clusters
  std::cout
      << "selection (all > MBDTriggerLive > GlobalVertex): "
      << _eventcounter << " > " << _eventcounter_selection1 << " > "
      << _eventcounter_selection2 << std::endl;

  std::cout << "selection clusters (all > (E>0.5, chi2<1000)): "
            << _clustercounter << " > " << _clustercounter_selection
            << std::endl;

  return 0;
}

void AnNeutralMeson::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  
  // Grab the CEMC node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  // Get the _det_name subnode
  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));

  // Check that it is there
  if (!cemcNode)
  {
    cemcNode = new PHCompositeNode("CEMC");
    dstNode->addNode(cemcNode);
  }
  ClusterSmallInfoNodeName = "CLUSTER_SMALLINFO_CEMC";
  
  _smallclusters = findNode::getClass<ClusterSmallInfoContainer>(dstNode, ClusterSmallInfoNodeName);
  
  if (!_smallclusters)
  {
    _smallclusters = new ClusterSmallInfoContainer();
  }

  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_smallclusters, ClusterSmallInfoNodeName, "PHObject");

  cemcNode->addNode(clusterNode);
}
