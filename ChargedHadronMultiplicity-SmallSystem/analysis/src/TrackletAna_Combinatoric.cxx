//

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TObjString.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TTreeIndex.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "GenHadron.h"
#include "LightHelixFit.h"
#include "Tracklet_combinatoric.h"
#include "Vertex.h"
#include "pdgidfunc.h"

bool verbose_debug = true;
bool savedetail = true;
bool dotruthmatching = false;
bool use_sPHENIX_primary_definition = true; // false: use original AuAu dNdEta primary definition

int main(int argc, char *argv[])
{
    if (argc < 8 || argc > 10)
    {
        cout << "Usage: ./TrackletAna_Combinatoric [isdata] [infile] [outfile] [NevtToRun] [pairCut] [clusadccutset] [clusphisizecutset] [pairCutMode(optional: dR|dPhi|absdPhi)] [minCrossingSilSeeds(optional)]" << endl;
        exit(1);
    }

    std::cout << "Primary definition = " << (use_sPHENIX_primary_definition ? "sPHENIXPrimary" : "PrimaryPHG4Ptcl (AuAu dNdEta)") << std::endl;

    // print out the input parameters
    for (int i = 0; i < argc; i++)
    {
        cout << "argv[" << i << "] = " << argv[i] << endl;
    }

    bool IsData = (TString(argv[1]).Atoi() == 1) ? true : false;
    TString infilename = TString(argv[2]);
    TString outfilename = TString(argv[3]);
    int NevtToRun_ = TString(argv[4]).Atoi();
    float pairCut = TString(argv[5]).Atof();
    int clusadccutset = TString(argv[6]).Atoi();
    int clusphisizecutset = TString(argv[7]).Atoi();

    TString idxstr = (IsData) ? "gl1bco" : "counter";

    TrackletData tkldata = {};

    // Setup for combinatoric tracklet method
    TrackletCombiDRConfig cfg;
    cfg.pair_cut = pairCut;
    auto parse_int_arg = [](const char *arg, const char *name) -> int
    {
        char *endptr = nullptr;
        const long value = std::strtol(arg, &endptr, 10);
        if (!endptr || *endptr != '\0')
        {
            throw std::invalid_argument(std::string("Invalid integer for ") + name + ": " + arg);
        }
        return static_cast<int>(value);
    };

    if (argc >= 9)
    {
        const std::string arg8 = argv[8];
        try
        {
            cfg.pair_cut_mode = ParsePairCutMode(arg8);
        }
        catch (const std::exception &e)
        {
            if (argc == 9)
            {
                try
                {
                    cfg.min_vertex_associated_silseeds = parse_int_arg(argv[8], "minCrossingSilSeeds");
                }
                catch (const std::exception &seed_e)
                {
                    std::cerr << e.what() << std::endl;
                    std::cerr << seed_e.what() << std::endl;
                    return 1;
                }
            }
            else
            {
                std::cerr << e.what() << std::endl;
                return 1;
            }
        }
    }
    if (argc == 10)
    {
        try
        {
            cfg.min_vertex_associated_silseeds = parse_int_arg(argv[9], "minCrossingSilSeeds");
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    // Output (eta,phi) grid binning: uniform or edges
    cfg.etaAxis = AxisBinning(-1.1, 1.1, 11);
    cfg.phiAxis = AxisBinning(-3.2, 3.2, 16);
    // Background rotation settings
    cfg.enable_eta_preselect = false;
    cfg.rotate_search_layer = true; // rotate logical search layer ({5,6})
    cfg.bkg_scale = 1.0;
    cfg.abs_dphi_count_cut = 0.15; // default PHOBOS-like scale; tune as needed
    // Pair-variable hist settings per (eta,phi) cell
    cfg.dR_nbins = 120;
    if (cfg.pair_cut_mode == TrackletCombiDRConfig::PairCutMode::kDeltaPhi)
    {
        cfg.dR_min = -cfg.pair_cut;
        cfg.dR_max = cfg.pair_cut;
    }
    else
    {
        cfg.dR_min = 0.0;
        cfg.dR_max = cfg.pair_cut;
    }
    cfg.enable_zvtx_cut = true;
    cfg.zvtx_cut_min = -10.0;
    cfg.zvtx_cut_max = 10.0;
    cfg.use_truth_vertex_for_zvtx_cut = false; // set to true for simulation study with truth vertex only for now
    if (cfg.use_truth_vertex_for_zvtx_cut)
    {
        std::cout << "Using truth vertex for z vertex cut with range [" << cfg.zvtx_cut_min << ", " << cfg.zvtx_cut_max << "] cm" << std::endl;
    }

    if (IsData && cfg.use_truth_vertex_for_zvtx_cut)
    {
        std::cout << "[WARNING] use_truth_vertex_for_zvtx_cut currently set to true but running with data. Setting to false" << std::endl;
        cfg.use_truth_vertex_for_zvtx_cut = false;
    }

    const std::string output_prefix = "tkl_Combinatoric";
    const std::string multiplicity_boundary_prefix = IsData ? "Data_" : "MC_";
    const std::string multiplicity_percentile_boundary_file = "figure/figure-NInttClusterCrossing/" + multiplicity_boundary_prefix + "NInttClusterPercentileBoundaries.root";
    const std::vector<MultiplicityPercentileBin> multiplicity_percentile_bins = LoadMultiplicityPercentileBins(multiplicity_percentile_boundary_file);

    std::cout << "Loaded " << multiplicity_percentile_bins.size() << " multiplicity percentile bins from " << multiplicity_percentile_boundary_file << std::endl;
    for (const auto &bin : multiplicity_percentile_bins)
    {
        std::cout << "  multPercentileBin" << bin.percentile_bin << ": [" << bin.cluster_low << ", ";
        if (bin.cluster_high == std::numeric_limits<int>::max())
            std::cout << "INT_MAX";
        else
            std::cout << bin.cluster_high;
        std::cout << ") INTT clusters" << std::endl;
    }

    TrackletCombiDROutput out;
    InitTrackletCombiDROutput(out, cfg, output_prefix, !IsData);

    std::vector<TrackletCombiDROutput> multiplicity_percentile_outputs(multiplicity_percentile_bins.size());
    for (size_t i = 0; i < multiplicity_percentile_bins.size(); ++i)
    {
        InitTrackletCombiDROutput(multiplicity_percentile_outputs[i], cfg, MultiplicityPercentileBinPrefix(output_prefix, multiplicity_percentile_bins[i].percentile_bin), !IsData);
    }

    std::cout << "Pair cut mode = " << PairCutModeLabel(cfg.pair_cut_mode) << ", threshold = " << cfg.pair_cut << std::endl;
    std::cout << "Direct |delta-phi| count cut for merged (eta, phi) subtracted histogram = " << cfg.abs_dphi_count_cut << std::endl;
    std::cout << "Minimum silicon seeds associated to selected crossing vertex = " << cfg.min_vertex_associated_silseeds;
    if (cfg.min_vertex_associated_silseeds >= 3)
        std::cout << " (active)";
    else
        std::cout << " (inactive; vertex finder already requires >=2)";
    std::cout << std::endl;

    // cluster adc cut
    int clusadccut = ConstADCCut(clusadccutset);
    // cluster phi-size cut
    float clusphisizecut = ConstClusPhiCut(clusphisizecutset);
    // print out cut info
    std::cout << "Cluster ADC cut (greater than): " << clusadccut << "; Cluster Phi-size cut (smaller than): " << clusphisizecut << std::endl;

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("VTX");
    t->BuildIndex(idxstr); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    const Long64_t nIndexed = index ? index->GetN() : 0;
    if (nIndexed <= 0)
    {
        std::cerr << "No entries available in tree index." << std::endl;
        f->Close();
        delete f;
        return 1;
    }
    if (NevtToRun_ > nIndexed)
    {
        std::cout << "Requested NevtToRun (" << NevtToRun_ << ") exceeds indexed entries (" << nIndexed << "), clamping." << std::endl;
        NevtToRun_ = static_cast<int>(nIndexed);
    }
    int event, nTruthVertex;
    std::vector<int> *firedTriggers = 0;
    uint64_t gl1bco, bcotr;
    bool is_min_bias;
    bool InttBco_IsToBeRemoved = false; // For this study we don't remove any BCO
    float mbd_south_charge_sum, mbd_north_charge_sum, mbd_charge_sum, mbd_charge_asymm, mbd_z_vtx;

    // vertex information
    int nSvtxVertices;
    std::vector<int> *trackerVertexId = 0;
    std::vector<float> *trackerVertexX = 0;
    std::vector<float> *trackerVertexY = 0;
    std::vector<float> *trackerVertexZ = 0;
    std::vector<int> *trackerVertexNTracks = 0;
    std::vector<std::vector<int>> *trackerVertexTrackIDs = 0;
    std::vector<short int> *trackerVertexCrossing = 0;
    // truth information for simulation
    std::vector<float> *TruthVertexX = 0;
    std::vector<float> *TruthVertexY = 0;
    std::vector<float> *TruthVertexZ = 0;

    std::vector<float> *silseed_x = 0;
    std::vector<float> *silseed_y = 0;
    std::vector<float> *silseed_z = 0;
    std::vector<float> *silseed_eta_vtx = 0;
    std::vector<float> *silseed_phi_vtx = 0;
    std::vector<int> *silseed_nMvtx = 0;
    std::vector<int> *silseed_nIntt = 0;
    std::vector<int> *silseed_crossing = 0;
    std::vector<std::vector<float>> *silseed_cluster_globalX = 0;
    std::vector<std::vector<float>> *silseed_cluster_globalY = 0;
    std::vector<std::vector<float>> *silseed_cluster_globalZ = 0;

    std::vector<int> *ClusLayer = 0;
    std::vector<float> *ClusX = 0, *ClusY = 0, *ClusZ = 0, *ClusPhiSize = 0, *ClusZSize = 0;
    std::vector<unsigned int> *ClusAdc = 0;
    std::vector<int> *cluster_timeBucketID = 0;
    std::vector<int> *ClusLadderZId = 0;

    std::vector<int> *cluster_matchedG4P_trackID = 0, *cluster_matchedG4P_PID = 0;
    std::vector<float> *cluster_matchedG4P_E = 0, *cluster_matchedG4P_pT = 0, *cluster_matchedG4P_eta = 0, *cluster_matchedG4P_phi = 0;

    // sPHENIX primary definition
    std::vector<int> *sPHENIXPrimary_PID = 0, *sPHENIXPrimary_trackID = 0;
    std::vector<float> *sPHENIXPrimary_pT = 0, *sPHENIXPrimary_eta = 0, *sPHENIXPrimary_phi = 0, *sPHENIXPrimary_E = 0;
    std::vector<double> *sPHENIXPrimary_charge = 0;
    std::vector<bool> *sPHENIXPrimary_isChargedHadron = 0;

    // Original AuAu dNdEta primary defintion
    std::vector<int> *PrimaryPHG4Ptcl_PID = 0, *PrimaryPHG4Ptcl_trackID = 0;
    std::vector<float> *PrimaryPHG4Ptcl_pT = 0, *PrimaryPHG4Ptcl_eta = 0, *PrimaryPHG4Ptcl_phi = 0, *PrimaryPHG4Ptcl_E = 0;
    std::vector<double> *PrimaryPHG4Ptcl_charge = 0;
    std::vector<bool> *PrimaryPHG4Ptcl_isChargedHadron = 0;

    // sPHENIX primary truth matching for PHG4Particles (truth clusters and the reco cluster that are matched to the truth cluster)
    std::vector<std::vector<float>> *sPHENIXPrimary_truthcluster_X = 0, *sPHENIXPrimary_truthcluster_Y = 0, *sPHENIXPrimary_truthcluster_Z = 0;
    std::vector<std::vector<float>> *sPHENIXPrimary_recocluster_globalX = 0, *sPHENIXPrimary_recocluster_globalY = 0, *sPHENIXPrimary_recocluster_globalZ = 0;

    std::vector<int> *G4P_trackID = 0;
    std::vector<float> *G4P_Pt = 0, *G4P_Eta = 0, *G4P_Phi = 0, *G4P_E = 0;

    t->SetBranchAddress("counter", &event);
    // t->SetBranchAddress("INTT_BCO", &INTT_BCO);
    t->SetBranchAddress("gl1bco", &gl1bco);
    t->SetBranchAddress("bcotr", &bcotr);
    t->SetBranchAddress("nSvtxVertices", &nSvtxVertices);
    t->SetBranchAddress("trackerVertexId", &trackerVertexId);
    t->SetBranchAddress("trackerVertexX", &trackerVertexX);
    t->SetBranchAddress("trackerVertexY", &trackerVertexY);
    t->SetBranchAddress("trackerVertexZ", &trackerVertexZ);
    t->SetBranchAddress("trackerVertexNTracks", &trackerVertexNTracks);
    t->SetBranchAddress("trackerVertexTrackIDs", &trackerVertexTrackIDs);
    t->SetBranchAddress("trackerVertexCrossing", &trackerVertexCrossing);
    if (!IsData)
    {
        t->SetBranchAddress("nTruthVertex", &nTruthVertex);
        t->SetBranchAddress("TruthVertexX", &TruthVertexX);
        t->SetBranchAddress("TruthVertexY", &TruthVertexY);
        t->SetBranchAddress("TruthVertexZ", &TruthVertexZ);
        if (use_sPHENIX_primary_definition)
        {
            t->SetBranchAddress("sPHENIXPrimary_PID", &sPHENIXPrimary_PID);
            t->SetBranchAddress("sPHENIXPrimary_isChargedHadron", &sPHENIXPrimary_isChargedHadron);
            t->SetBranchAddress("sPHENIXPrimary_trackID", &sPHENIXPrimary_trackID);
            t->SetBranchAddress("sPHENIXPrimary_pT", &sPHENIXPrimary_pT);
            t->SetBranchAddress("sPHENIXPrimary_eta", &sPHENIXPrimary_eta);
            t->SetBranchAddress("sPHENIXPrimary_phi", &sPHENIXPrimary_phi);
            t->SetBranchAddress("sPHENIXPrimary_E", &sPHENIXPrimary_E);
        }
        else
        {
            t->SetBranchAddress("PrimaryPHG4Ptcl_PID", &PrimaryPHG4Ptcl_PID);
            t->SetBranchAddress("PrimaryPHG4Ptcl_isChargedHadron", &PrimaryPHG4Ptcl_isChargedHadron);
            t->SetBranchAddress("PrimaryPHG4Ptcl_trackID", &PrimaryPHG4Ptcl_trackID);
            t->SetBranchAddress("PrimaryPHG4Ptcl_pT", &PrimaryPHG4Ptcl_pT);
            t->SetBranchAddress("PrimaryPHG4Ptcl_eta", &PrimaryPHG4Ptcl_eta);
            t->SetBranchAddress("PrimaryPHG4Ptcl_phi", &PrimaryPHG4Ptcl_phi);
            t->SetBranchAddress("PrimaryPHG4Ptcl_E", &PrimaryPHG4Ptcl_E);
            t->SetBranchAddress("PrimaryPHG4Ptcl_charge", &PrimaryPHG4Ptcl_charge);
        }

        // if (dotruthmatching)
        {
            t->SetBranchAddress("cluster_matchedG4P_trackID", &cluster_matchedG4P_trackID);
            t->SetBranchAddress("cluster_matchedG4P_PID", &cluster_matchedG4P_PID);
            t->SetBranchAddress("cluster_matchedG4P_E", &cluster_matchedG4P_E);
            t->SetBranchAddress("cluster_matchedG4P_pT", &cluster_matchedG4P_pT);
            t->SetBranchAddress("cluster_matchedG4P_eta", &cluster_matchedG4P_eta);
            t->SetBranchAddress("cluster_matchedG4P_phi", &cluster_matchedG4P_phi);
            //     t->SetBranchAddress("G4P_trackID", &G4P_trackID);
            //     t->SetBranchAddress("G4P_Pt", &G4P_Pt);
            //     t->SetBranchAddress("G4P_Eta", &G4P_Eta);
            //     t->SetBranchAddress("G4P_Phi", &G4P_Phi);
            //     t->SetBranchAddress("G4P_E", &G4P_E);
        }

        t->SetBranchAddress("sPHENIXPrimary_charge", &sPHENIXPrimary_charge);
        t->SetBranchAddress("sPHENIXPrimary_truthcluster_X", &sPHENIXPrimary_truthcluster_X);
        t->SetBranchAddress("sPHENIXPrimary_truthcluster_Y", &sPHENIXPrimary_truthcluster_Y);
        t->SetBranchAddress("sPHENIXPrimary_truthcluster_Z", &sPHENIXPrimary_truthcluster_Z);
        t->SetBranchAddress("sPHENIXPrimary_recocluster_globalX", &sPHENIXPrimary_recocluster_globalX);
        t->SetBranchAddress("sPHENIXPrimary_recocluster_globalY", &sPHENIXPrimary_recocluster_globalY);
        t->SetBranchAddress("sPHENIXPrimary_recocluster_globalZ", &sPHENIXPrimary_recocluster_globalZ);

        InttBco_IsToBeRemoved = false;
    }
    if (IsData)
    {
        t->SetBranchAddress("firedTriggers", &firedTriggers);
        // t->SetBranchAddress("InttBco_IsToBeRemoved", &InttBco_IsToBeRemoved);
    }
    // t->SetBranchAddress("is_min_bias", &is_min_bias);
    // t->SetBranchAddress("MBD_centrality", &centrality_mbd);
    // t->SetBranchAddress("MBD_z_vtx", &mbd_z_vtx);
    // t->SetBranchAddress("MBD_south_charge_sum", &mbd_south_charge_sum);
    // t->SetBranchAddress("MBD_north_charge_sum", &mbd_north_charge_sum);
    t->SetBranchAddress("MBD_charge_sum", &mbd_charge_sum);
    // t->SetBranchAddress("MBD_charge_asymm", &mbd_charge_asymm);
    t->SetBranchAddress("cluster_layer", &ClusLayer);
    t->SetBranchAddress("cluster_globalX", &ClusX);
    t->SetBranchAddress("cluster_globalY", &ClusY);
    t->SetBranchAddress("cluster_globalZ", &ClusZ);
    t->SetBranchAddress("cluster_phiSize", &ClusPhiSize);
    t->SetBranchAddress("cluster_zSize", &ClusZSize);
    t->SetBranchAddress("cluster_adc", &ClusAdc);
    t->SetBranchAddress("cluster_timeBucketID", &cluster_timeBucketID);
    t->SetBranchAddress("cluster_ladderZId", &ClusLadderZId);

    t->SetBranchAddress("silseed_eta_vtx", &silseed_eta_vtx);
    t->SetBranchAddress("silseed_phi_vtx", &silseed_phi_vtx);
    t->SetBranchAddress("silseed_nMvtx", &silseed_nMvtx);
    t->SetBranchAddress("silseed_nIntt", &silseed_nIntt);
    const bool has_silseed_crossing_branch = static_cast<bool>(t->GetBranch("silseed_crossing"));
    if (has_silseed_crossing_branch)
    {
        t->SetBranchAddress("silseed_crossing", &silseed_crossing);
    }
    else
    {
        std::cout << "[INFO] silseed_crossing branch is not available; all-seed crossing histograms will remain empty." << std::endl;
    }

    const bool has_silseed_xyz_branches = (t->GetBranch("silseed_x") && t->GetBranch("silseed_y") && t->GetBranch("silseed_z"));
    if (has_silseed_xyz_branches)
    {
        t->SetBranchAddress("silseed_x", &silseed_x);
        t->SetBranchAddress("silseed_y", &silseed_y);
        t->SetBranchAddress("silseed_z", &silseed_z);
    }
    else
    {
        std::cout << "[INFO] silseed_x/y/z branches are not available; seed-PCA display graphs will be skipped." << std::endl;
    }

    const bool has_silseed_cluster_xyz_branches = (t->GetBranch("silseed_cluster_globalX") && t->GetBranch("silseed_cluster_globalY") && t->GetBranch("silseed_cluster_globalZ"));
    if (has_silseed_cluster_xyz_branches)
    {
        t->SetBranchAddress("silseed_cluster_globalX", &silseed_cluster_globalX);
        t->SetBranchAddress("silseed_cluster_globalY", &silseed_cluster_globalY);
        t->SetBranchAddress("silseed_cluster_globalZ", &silseed_cluster_globalZ);
    }
    else
    {
        std::cout << "[INFO] silseed_cluster_globalX/Y/Z branches are not available; per-seed cluster display graphs will be skipped." << std::endl;
    }

    std::unordered_map<int, std::vector<int>> vertexCrossing_idx_map; // key: vertex crossing, value: vertex index in the tree
    for (int i = 0; i < NevtToRun_; i++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);

        std::vector<float> *selectedPrimaryEta = use_sPHENIX_primary_definition ? sPHENIXPrimary_eta : PrimaryPHG4Ptcl_eta;
        std::vector<float> *selectedPrimaryPhi = use_sPHENIX_primary_definition ? sPHENIXPrimary_phi : PrimaryPHG4Ptcl_phi;
        std::vector<bool> *selectedPrimaryIsChargedHadron = use_sPHENIX_primary_definition ? sPHENIXPrimary_isChargedHadron : PrimaryPHG4Ptcl_isChargedHadron;
        bool filledPrimaryChargedAfterZvtxCut = false;
        std::vector<char> filledPrimaryChargedAfterZvtxCutByMultBin(multiplicity_percentile_outputs.size(), 0);

        auto fill_selected_crossing_histograms = [&](TrackletCombiDROutput &target_out, auto &filled_primary_charged, int crossing_value, int vtx_idx_value, const std::vector<size_t> &used_seed_indices_value, int n_associated_silseeds_value, int n_intt_clusters_in_crossing_value, int n_all_silseeds_in_crossing_value, double truth_vtxz_value, bool has_truth_vtxz_value)
        {
            if (target_out.hSeedCrossing_selectedCrossings)
                target_out.hSeedCrossing_selectedCrossings->Fill(crossing_value, 1.0);

            if (target_out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds && n_intt_clusters_in_crossing_value > 4 * n_associated_silseeds_value)
                target_out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds->Fill(crossing_value, 1.0);

            if (!filled_primary_charged && !IsData && target_out.hPrimaryCharged_Eta_Phi && selectedPrimaryEta && selectedPrimaryPhi && selectedPrimaryIsChargedHadron)
            {
                const size_t nprim = std::min({selectedPrimaryEta->size(), selectedPrimaryPhi->size(), selectedPrimaryIsChargedHadron->size()});
                for (size_t iprim = 0; iprim < nprim; ++iprim)
                {
                    if (!selectedPrimaryIsChargedHadron->at(iprim))
                        continue;
                    const double eta_primary = selectedPrimaryEta->at(iprim);
                    double phi_primary = selectedPrimaryPhi->at(iprim);
                    phi_primary = target_out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_primary) : wrapPhi_mpi_pi(phi_primary);
                    target_out.hPrimaryCharged_Eta_Phi->Fill(eta_primary, phi_primary, 1.0);
                }
                filled_primary_charged = true;
            }

            if (target_out.hSilSeed_Eta_Phi && trackerVertexTrackIDs && vtx_idx_value >= 0 && vtx_idx_value < static_cast<int>(trackerVertexTrackIDs->size()) && silseed_eta_vtx && silseed_phi_vtx)
            {
                for (const size_t iseed : used_seed_indices_value)
                {
                    if (iseed >= silseed_eta_vtx->size() || iseed >= silseed_phi_vtx->size())
                        continue;

                    const double eta_seed_vtx = silseed_eta_vtx->at(iseed);
                    double phi_seed_vtx = silseed_phi_vtx->at(iseed);
                    phi_seed_vtx = target_out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_seed_vtx) : wrapPhi_mpi_pi(phi_seed_vtx);
                    target_out.hSilSeed_Eta_Phi->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);

                    if (target_out.hSilSeedDCA3D && has_silseed_xyz_branches && silseed_x && silseed_y && silseed_z && iseed < silseed_x->size() && iseed < silseed_y->size() && iseed < silseed_z->size())
                    {
                        const double dx = silseed_x->at(iseed) - trackerVertexX->at(vtx_idx_value);
                        const double dy = silseed_y->at(iseed) - trackerVertexY->at(vtx_idx_value);
                        const double dz = silseed_z->at(iseed) - trackerVertexZ->at(vtx_idx_value);
                        const double silseed_dca3d = std::sqrt(dx * dx + dy * dy + dz * dz);
                        target_out.hSilSeedDCA3D->Fill(silseed_dca3d, 1.0);

                        if (has_truth_vtxz_value && !target_out.hSilSeedDCA3D_truthVtxZ.empty())
                        {
                            const int truth_vtxz_bin = FindTruthVtxZBin(cfg.truth_vtxz_bin_edges, truth_vtxz_value);
                            if (truth_vtxz_bin >= 0 && truth_vtxz_bin < static_cast<int>(target_out.hSilSeedDCA3D_truthVtxZ.size()) && target_out.hSilSeedDCA3D_truthVtxZ[truth_vtxz_bin])
                                target_out.hSilSeedDCA3D_truthVtxZ[truth_vtxz_bin]->Fill(silseed_dca3d, 1.0);
                        }
                    }

                    if (silseed_nMvtx && silseed_nIntt && iseed < silseed_nMvtx->size() && iseed < silseed_nIntt->size())
                    {
                        const int nMvtx = silseed_nMvtx->at(iseed);
                        const int nIntt = silseed_nIntt->at(iseed);

                        if (nMvtx == 3 && nIntt == 2 && target_out.hSilSeed_Eta_Phi_nMVTX3nINTT2)
                            target_out.hSilSeed_Eta_Phi_nMVTX3nINTT2->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                        else if (nMvtx == 3 && nIntt == 1 && target_out.hSilSeed_Eta_Phi_nMVTX3nINTT1)
                            target_out.hSilSeed_Eta_Phi_nMVTX3nINTT1->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                        else if (nMvtx == 2 && nIntt == 2 && target_out.hSilSeed_Eta_Phi_nMVTX2nINTT2)
                            target_out.hSilSeed_Eta_Phi_nMVTX2nINTT2->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                        else if (nMvtx == 2 && nIntt == 1 && target_out.hSilSeed_Eta_Phi_nMVTX2nINTT1)
                            target_out.hSilSeed_Eta_Phi_nMVTX2nINTT1->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                        else if (target_out.hSilSeed_Eta_Phi_nMVTXnINTTOther)
                            target_out.hSilSeed_Eta_Phi_nMVTXnINTTOther->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                    }
                }
            }

            if (target_out.hINTTClusterVsAssociatedSilSeed)
                target_out.hINTTClusterVsAssociatedSilSeed->Fill(n_intt_clusters_in_crossing_value, n_associated_silseeds_value, 1.0);

            if (target_out.hINTTClusterVsAllSilSeed && n_all_silseeds_in_crossing_value >= 0)
                target_out.hINTTClusterVsAllSilSeed->Fill(n_intt_clusters_in_crossing_value, n_all_silseeds_in_crossing_value, 1.0);
        };

        cout << "i = " << i << " event = " << event << " GL1 BCO = " << gl1bco << " NevtToRun_ = " << NevtToRun_ << " local = " << local << endl;

        // check how many unique vertex crossing in this trigger frame
        for (int j = 0; j < trackerVertexCrossing->size(); j++)
        {
            vertexCrossing_idx_map[trackerVertexCrossing->at(j)].push_back(j);
        }
        cout << "Number of unique vertex crossings in this trigger frame = " << vertexCrossing_idx_map.size() << endl;

        // only consider the crossing with only 1 silicon vertex to avoid in-time pileup contamination
        for (auto &pair : vertexCrossing_idx_map)
        {
            auto vertex_crossing = pair.first;
            // if there is only 1 vertex with this crossing AND this crossing is not short int max, proceed with the analysis
            if (pair.second.size() == 1 && vertex_crossing != std::numeric_limits<short int>::max())
            {
                out.nSelectedCrossings++;

                int vtx_idx = pair.second[0];
                const bool reco_vertex_vectors_ok =
                    trackerVertexNTracks && trackerVertexX && trackerVertexY && trackerVertexZ &&
                    vtx_idx >= 0 &&
                    vtx_idx < static_cast<int>(trackerVertexNTracks->size()) &&
                    vtx_idx < static_cast<int>(trackerVertexX->size()) &&
                    vtx_idx < static_cast<int>(trackerVertexY->size()) &&
                    vtx_idx < static_cast<int>(trackerVertexZ->size());
                if (!reco_vertex_vectors_ok)
                {
                    std::cerr << "[ERROR] Reconstructed vertex index out of range for crossing " << vertex_crossing
                              << ": vtx_idx=" << vtx_idx
                              << ", trackerVertexNTracks size=" << (trackerVertexNTracks ? trackerVertexNTracks->size() : 0)
                              << ", trackerVertexX size=" << (trackerVertexX ? trackerVertexX->size() : 0)
                              << ", trackerVertexY size=" << (trackerVertexY ? trackerVertexY->size() : 0)
                              << ", trackerVertexZ size=" << (trackerVertexZ ? trackerVertexZ->size() : 0)
                              << std::endl;
                    continue;
                }
                std::cout << "Vertex crossing " << vertex_crossing << " has 1 vertex with index " << vtx_idx << " and Ntracks = " << trackerVertexNTracks->at(vtx_idx) << std::endl;

                // get the silicon vertex position
                const bool truth_vertex_vectors_ok =
                    TruthVertexX && TruthVertexY && TruthVertexZ &&
                    vtx_idx >= 0 &&
                    vtx_idx < static_cast<int>(TruthVertexX->size()) &&
                    vtx_idx < static_cast<int>(TruthVertexY->size()) &&
                    vtx_idx < static_cast<int>(TruthVertexZ->size());
                if (cfg.use_truth_vertex_for_zvtx_cut && !truth_vertex_vectors_ok)
                {
                    std::cerr << "[ERROR] Truth vertex index out of range for crossing " << vertex_crossing
                              << ": vtx_idx=" << vtx_idx
                              << ", TruthVertexX size=" << (TruthVertexX ? TruthVertexX->size() : 0)
                              << ", TruthVertexY size=" << (TruthVertexY ? TruthVertexY->size() : 0)
                              << ", TruthVertexZ size=" << (TruthVertexZ ? TruthVertexZ->size() : 0)
                              << ". Falling back to reconstructed vertex for z-vertex selection." << std::endl;
                }

                const bool use_truth_vertex_here = cfg.use_truth_vertex_for_zvtx_cut && truth_vertex_vectors_ok;
                float silvtx_x = use_truth_vertex_here ? TruthVertexX->at(vtx_idx) : trackerVertexX->at(vtx_idx);
                float silvtx_y = use_truth_vertex_here ? TruthVertexY->at(vtx_idx) : trackerVertexY->at(vtx_idx);
                float silvtx_z = use_truth_vertex_here ? TruthVertexZ->at(vtx_idx) : trackerVertexZ->at(vtx_idx);
                const std::array<double, 3> silicon_vtx = {silvtx_x, silvtx_y, silvtx_z};

                std::vector<size_t> used_seed_indices; // indices of seeds associated to this vertex crossing (via trackerVertexTrackIDs)
                if (trackerVertexTrackIDs && vtx_idx >= 0 && vtx_idx < static_cast<int>(trackerVertexTrackIDs->size()) && silseed_eta_vtx && silseed_phi_vtx)
                {
                    const auto &vtx_track_ids = trackerVertexTrackIDs->at(vtx_idx);
                    std::cout << "There are " << vtx_track_ids.size() << " Track IDs associated to this vertex." << std::endl;
                    for (const int track_id : vtx_track_ids)
                    {
                        used_seed_indices.push_back(track_id);
                    }
                }

                std::cout << "Vertex crossing " << vertex_crossing << " has " << used_seed_indices.size() << " associated silicon seeds." << std::endl;

                // Apply z vertex cut before filling selected-crossing diagnostics
                if (cfg.enable_zvtx_cut && (silvtx_z < cfg.zvtx_cut_min || silvtx_z > cfg.zvtx_cut_max))
                {
                    std::cout << "Skipping vertex crossing " << vertex_crossing << " due to z vertex cut: " << silvtx_z << " not in [" << cfg.zvtx_cut_min << ", " << cfg.zvtx_cut_max << "]" << std::endl;
                    continue;
                }

                std::vector<size_t> crossing_cluster_indices;
                crossing_cluster_indices.reserve(ClusLayer->size());
                for (size_t ihit = 0; ihit < ClusLayer->size(); ++ihit)
                {
                    if (cluster_timeBucketID->at(ihit) != vertex_crossing)
                        continue;
                    if (logicalInttLayerIndex(ClusLayer->at(ihit)) < 0)
                        continue;
                    crossing_cluster_indices.push_back(ihit);
                }

                const int n_associated_silseeds = static_cast<int>(used_seed_indices.size());
                const int n_intt_clusters_in_crossing = static_cast<int>(crossing_cluster_indices.size());
                const int n_all_silseeds_in_crossing = (has_silseed_crossing_branch && silseed_crossing) ? static_cast<int>(std::count(silseed_crossing->begin(), silseed_crossing->end(), static_cast<int>(vertex_crossing))) : -1;
                const bool has_truth_vtxz = (!IsData && TruthVertexZ && vtx_idx >= 0 && vtx_idx < static_cast<int>(TruthVertexZ->size()));
                const double truth_vtxz = has_truth_vtxz ? TruthVertexZ->at(vtx_idx) : 0.0;
                const int multiplicity_percentile_bin_index = FindMultiplicityPercentileBinIndex(multiplicity_percentile_bins, n_intt_clusters_in_crossing);

                if (multiplicity_percentile_bin_index >= 0)
                    multiplicity_percentile_outputs[multiplicity_percentile_bin_index].nSelectedCrossings++;
                else
                    std::cout << "Skipping multiplicity-binned fill for crossing " << vertex_crossing << " because N_{INTT clusters} = " << n_intt_clusters_in_crossing << " did not match any configured percentile interval." << std::endl;

                fill_selected_crossing_histograms(out, filledPrimaryChargedAfterZvtxCut, vertex_crossing, vtx_idx, used_seed_indices, n_associated_silseeds, n_intt_clusters_in_crossing, n_all_silseeds_in_crossing, truth_vtxz, has_truth_vtxz);
                if (multiplicity_percentile_bin_index >= 0)
                {
                    fill_selected_crossing_histograms(multiplicity_percentile_outputs[multiplicity_percentile_bin_index], filledPrimaryChargedAfterZvtxCutByMultBin[multiplicity_percentile_bin_index], vertex_crossing, vtx_idx, used_seed_indices, n_associated_silseeds, n_intt_clusters_in_crossing, n_all_silseeds_in_crossing, truth_vtxz, has_truth_vtxz);
                }

                if (cfg.min_vertex_associated_silseeds >= 3 && n_associated_silseeds < cfg.min_vertex_associated_silseeds)
                {
                    std::cout << "Skipping vertex crossing " << vertex_crossing << " because it has only " << used_seed_indices.size() << " associated silicon seeds, below the requested minimum of " << cfg.min_vertex_associated_silseeds << std::endl;
                    continue;
                }

                // const bool make_crossing_display = (n_associated_silseeds >= 4 && n_intt_clusters_in_crossing > 4 * n_associated_silseeds);
                // const bool make_crossing_display = true; // For debugging purposes we want to make the crossing display for every crossing that passes the z vertex cut, regardless of the number of associated seeds or clusters
                const bool make_crossing_display = false; // 
                TrackletCombiDROutput::CrossingDisplayGraphs display;
                std::vector<size_t> nonassociated_seed_indices;
                std::unordered_set<int> sphenix_primary_track_ids;
                if (make_crossing_display)
                {
                    display.crossing_key = IsData ? (std::to_string(gl1bco) + "_" + std::to_string(vertex_crossing)) : (std::to_string(event) + "_" + std::to_string(vertex_crossing));
                    if (has_silseed_crossing_branch && silseed_crossing)
                    {
                        std::vector<size_t> associated_seed_indices_sorted(used_seed_indices.begin(), used_seed_indices.end());
                        std::sort(associated_seed_indices_sorted.begin(), associated_seed_indices_sorted.end());

                        for (size_t iseed = 0; iseed < silseed_crossing->size(); ++iseed)
                        {
                            if (silseed_crossing->at(iseed) != static_cast<int>(vertex_crossing))
                                continue;
                            if (std::binary_search(associated_seed_indices_sorted.begin(), associated_seed_indices_sorted.end(), iseed))
                                continue;
                            nonassociated_seed_indices.push_back(iseed);
                        }
                    }

                    display.grAllINTTClustersXY = new TGraph();
                    display.grAllINTTClustersXY->SetName("gr_all_intt_clusters_xy");
                    display.grAllINTTClustersXY->SetTitle("All INTT clusters in crossing (x-y)");

                    display.grAllINTTClustersZR = new TGraph();
                    display.grAllINTTClustersZR->SetName("gr_all_intt_clusters_zr");
                    display.grAllINTTClustersZR->SetTitle("All INTT clusters in crossing (z-r)");

                    if (!IsData)
                    {
                        display.grPrimaryINTTClustersXY = new TGraph();
                        display.grPrimaryINTTClustersXY->SetName("gr_primary_intt_clusters_xy");
                        display.grPrimaryINTTClustersXY->SetTitle("INTT clusters from sPHENIX primary particles in crossing (x-y)");

                        display.grPrimaryINTTClustersZR = new TGraph();
                        display.grPrimaryINTTClustersZR->SetName("gr_primary_intt_clusters_zr");
                        display.grPrimaryINTTClustersZR->SetTitle("INTT clusters from sPHENIX primary particles in crossing (z-r)");

                        display.grNonPrimaryINTTClustersXY = new TGraph();
                        display.grNonPrimaryINTTClustersXY->SetName("gr_nonprimary_intt_clusters_xy");
                        display.grNonPrimaryINTTClustersXY->SetTitle("INTT clusters not from sPHENIX primary particles in crossing (x-y)");

                        display.grNonPrimaryINTTClustersZR = new TGraph();
                        display.grNonPrimaryINTTClustersZR->SetName("gr_nonprimary_intt_clusters_zr");
                        display.grNonPrimaryINTTClustersZR->SetTitle("INTT clusters not from sPHENIX primary particles in crossing (z-r)");

                        if (sPHENIXPrimary_trackID)
                            sphenix_primary_track_ids.insert(sPHENIXPrimary_trackID->begin(), sPHENIXPrimary_trackID->end());
                    }

                    display.grRecoVertexXY = new TGraph();
                    display.grRecoVertexXY->SetName("gr_reco_vertex_xy");
                    display.grRecoVertexXY->SetTitle("Selected reconstructed vertex (x-y)");

                    display.grRecoVertexZR = new TGraph();
                    display.grRecoVertexZR->SetName("gr_reco_vertex_zr");
                    display.grRecoVertexZR->SetTitle("Selected reconstructed vertex (z-r)");

                    const double reco_r_signed = (silvtx_y < 0.0) ? -std::hypot(silvtx_x, silvtx_y) : std::hypot(silvtx_x, silvtx_y);
                    display.grRecoVertexXY->SetPoint(0, silvtx_x, silvtx_y);
                    display.grRecoVertexZR->SetPoint(0, silvtx_z, reco_r_signed);

                    if (!IsData && TruthVertexX && TruthVertexY && TruthVertexZ && vtx_idx < static_cast<int>(TruthVertexX->size()) && vtx_idx < static_cast<int>(TruthVertexY->size()) && vtx_idx < static_cast<int>(TruthVertexZ->size()))
                    {
                        display.grTruthVertexXY = new TGraph();
                        display.grTruthVertexXY->SetName("gr_truth_vertex_xy");
                        display.grTruthVertexXY->SetTitle("Truth vertex (x-y)");

                        display.grTruthVertexZR = new TGraph();
                        display.grTruthVertexZR->SetName("gr_truth_vertex_zr");
                        display.grTruthVertexZR->SetTitle("Truth vertex (z-r)");

                        const double truth_x = TruthVertexX->at(vtx_idx);
                        const double truth_y = TruthVertexY->at(vtx_idx);
                        const double truth_z = TruthVertexZ->at(vtx_idx);
                        const double truth_r_signed = (truth_y < 0.0) ? -std::hypot(truth_x, truth_y) : std::hypot(truth_x, truth_y);
                        display.grTruthVertexXY->SetPoint(0, truth_x, truth_y);
                        display.grTruthVertexZR->SetPoint(0, truth_z, truth_r_signed);
                    }
                }

                for (const size_t ihit : crossing_cluster_indices)
                {
                    const int logical_layer = logicalInttLayerIndex(ClusLayer->at(ihit));
                    Hit *hit = new Hit(ClusX->at(ihit), ClusY->at(ihit), ClusZ->at(ihit), silvtx_x, silvtx_y, silvtx_z, logical_layer, ClusPhiSize->at(ihit), ClusAdc->at(ihit));
                    hit->SetEdgeFromZId(ClusLadderZId->at(ihit));
                    tkldata.layers[logical_layer].push_back(hit);

                    if (!make_crossing_display)
                        continue;

                    const double x = ClusX->at(ihit);
                    const double y = ClusY->at(ihit);
                    const double z = ClusZ->at(ihit);
                    const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);

                    display.grAllINTTClustersXY->SetPoint(display.grAllINTTClustersXY->GetN(), x, y);
                    display.grAllINTTClustersZR->SetPoint(display.grAllINTTClustersZR->GetN(), z, r_signed);

                    if (!IsData && cluster_matchedG4P_trackID && ihit < cluster_matchedG4P_trackID->size())
                    {
                        const int matched_track_id = cluster_matchedG4P_trackID->at(ihit);
                        const bool is_sphenix_primary_cluster = (matched_track_id != 0) && (sphenix_primary_track_ids.find(matched_track_id) != sphenix_primary_track_ids.end());

                        if (is_sphenix_primary_cluster)
                        {
                            if (display.grPrimaryINTTClustersXY)
                                display.grPrimaryINTTClustersXY->SetPoint(display.grPrimaryINTTClustersXY->GetN(), x, y);
                            if (display.grPrimaryINTTClustersZR)
                                display.grPrimaryINTTClustersZR->SetPoint(display.grPrimaryINTTClustersZR->GetN(), z, r_signed);
                        }
                        else
                        {
                            if (display.grNonPrimaryINTTClustersXY)
                                display.grNonPrimaryINTTClustersXY->SetPoint(display.grNonPrimaryINTTClustersXY->GetN(), x, y);
                            if (display.grNonPrimaryINTTClustersZR)
                                display.grNonPrimaryINTTClustersZR->SetPoint(display.grNonPrimaryINTTClustersZR->GetN(), z, r_signed);
                        }
                    }
                }

                std::cout << "Number of clusters in this crossing (inner, outer) = (" << tkldata.layers[0].size() << ", " << tkldata.layers[1].size() << ")" << std::endl;

                // Prepare crossing display graphs for crossings with high occupancy relative to associated seeds.
                if (make_crossing_display)
                {

                    if (has_silseed_xyz_branches && silseed_x && silseed_y && silseed_z)
                    {
                        display.grAssociatedSilSeedPCAXY = new TGraph();
                        display.grAssociatedSilSeedPCAXY->SetName("gr_associated_silseed_pca_xy");
                        display.grAssociatedSilSeedPCAXY->SetTitle("Associated silicon seed PCA points (x-y)");

                        display.grAssociatedSilSeedPCAZR = new TGraph();
                        display.grAssociatedSilSeedPCAZR->SetName("gr_associated_silseed_pca_zr");
                        display.grAssociatedSilSeedPCAZR->SetTitle("Associated silicon seed PCA points (z-r)");

                        for (const size_t iseed : used_seed_indices)
                        {
                            if (iseed >= silseed_x->size() || iseed >= silseed_y->size() || iseed >= silseed_z->size())
                                continue;

                            const double sx = silseed_x->at(iseed);
                            const double sy = silseed_y->at(iseed);
                            const double sz = silseed_z->at(iseed);
                            const double sr_signed = (sy < 0.0) ? -std::hypot(sx, sy) : std::hypot(sx, sy);

                            display.grAssociatedSilSeedPCAXY->SetPoint(display.grAssociatedSilSeedPCAXY->GetN(), sx, sy);
                            display.grAssociatedSilSeedPCAZR->SetPoint(display.grAssociatedSilSeedPCAZR->GetN(), sz, sr_signed);
                        }
                    }

                    if (has_silseed_cluster_xyz_branches && silseed_cluster_globalX && silseed_cluster_globalY && silseed_cluster_globalZ)
                    {
                        std::vector<size_t> seed_indices_sorted(used_seed_indices.begin(), used_seed_indices.end());
                        std::sort(seed_indices_sorted.begin(), seed_indices_sorted.end());

                        for (const size_t iseed : seed_indices_sorted)
                        {
                            if (iseed >= silseed_cluster_globalX->size() || iseed >= silseed_cluster_globalY->size() || iseed >= silseed_cluster_globalZ->size())
                                continue;

                            const auto &seed_cl_x = silseed_cluster_globalX->at(iseed);
                            const auto &seed_cl_y = silseed_cluster_globalY->at(iseed);
                            const auto &seed_cl_z = silseed_cluster_globalZ->at(iseed);

                            TGraph *gr_seed_xy = new TGraph();
                            gr_seed_xy->SetName(Form("gr_associated_silseed_clusters_xy_idx%zu", iseed));
                            gr_seed_xy->SetTitle(Form("Associated silicon seed idx=%zu clusters (x-y)", iseed));

                            TGraph *gr_seed_zr = new TGraph();
                            gr_seed_zr->SetName(Form("gr_associated_silseed_clusters_zr_idx%zu", iseed));
                            gr_seed_zr->SetTitle(Form("Associated silicon seed idx=%zu clusters (z-r)", iseed));

                            const size_t nseedclus = std::min({seed_cl_x.size(), seed_cl_y.size(), seed_cl_z.size()});
                            for (size_t ic = 0; ic < nseedclus; ++ic)
                            {
                                const double x = seed_cl_x[ic];
                                const double y = seed_cl_y[ic];
                                const double z = seed_cl_z[ic];
                                const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);

                                gr_seed_xy->SetPoint(gr_seed_xy->GetN(), x, y);
                                gr_seed_zr->SetPoint(gr_seed_zr->GetN(), z, r_signed);
                            }

                            display.grAssociatedSilSeedClustersXY.push_back(gr_seed_xy);
                            display.grAssociatedSilSeedClustersZR.push_back(gr_seed_zr);
                        }
                    }

                    if (has_silseed_cluster_xyz_branches && silseed_cluster_globalX && silseed_cluster_globalY && silseed_cluster_globalZ)
                    {
                        for (const size_t iseed : nonassociated_seed_indices)
                        {
                            if (iseed >= silseed_cluster_globalX->size() || iseed >= silseed_cluster_globalY->size() || iseed >= silseed_cluster_globalZ->size())
                                continue;

                            const auto &seed_cl_x = silseed_cluster_globalX->at(iseed);
                            const auto &seed_cl_y = silseed_cluster_globalY->at(iseed);
                            const auto &seed_cl_z = silseed_cluster_globalZ->at(iseed);

                            TGraph *gr_seed_xy = new TGraph();
                            gr_seed_xy->SetName(Form("gr_nonassociated_silseed_clusters_xy_idx%zu", iseed));
                            gr_seed_xy->SetTitle(Form("Non-associated silicon seed idx=%zu clusters (x-y)", iseed));

                            TGraph *gr_seed_zr = new TGraph();
                            gr_seed_zr->SetName(Form("gr_nonassociated_silseed_clusters_zr_idx%zu", iseed));
                            gr_seed_zr->SetTitle(Form("Non-associated silicon seed idx=%zu clusters (z-r)", iseed));

                            const size_t nseedclus = std::min({seed_cl_x.size(), seed_cl_y.size(), seed_cl_z.size()});
                            for (size_t ic = 0; ic < nseedclus; ++ic)
                            {
                                const double x = seed_cl_x[ic];
                                const double y = seed_cl_y[ic];
                                const double z = seed_cl_z[ic];
                                const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);

                                gr_seed_xy->SetPoint(gr_seed_xy->GetN(), x, y);
                                gr_seed_zr->SetPoint(gr_seed_zr->GetN(), z, r_signed);
                            }

                            display.grNonAssociatedSilSeedClustersXY.push_back(gr_seed_xy);
                            display.grNonAssociatedSilSeedClustersZR.push_back(gr_seed_zr);
                        }
                    }

                    if (!IsData && sPHENIXPrimary_truthcluster_X && sPHENIXPrimary_truthcluster_Y && sPHENIXPrimary_truthcluster_Z)
                    {
                        const size_t nprimary_truth = std::min({sPHENIXPrimary_truthcluster_X->size(), sPHENIXPrimary_truthcluster_Y->size(), sPHENIXPrimary_truthcluster_Z->size()});
                        for (size_t iprim = 0; iprim < nprimary_truth; ++iprim)
                        {
                            if (!sPHENIXPrimary_charge || iprim >= sPHENIXPrimary_charge->size() || sPHENIXPrimary_charge->at(iprim) == 0.0)
                                continue;
                            const auto &truth_cl_x = sPHENIXPrimary_truthcluster_X->at(iprim);
                            const auto &truth_cl_y = sPHENIXPrimary_truthcluster_Y->at(iprim);
                            const auto &truth_cl_z = sPHENIXPrimary_truthcluster_Z->at(iprim);

                            TGraph *gr_truth_xy = new TGraph();
                            gr_truth_xy->SetName(Form("gr_sPHENIXPrimary_truthcluster_xy_idx%zu", iprim));
                            gr_truth_xy->SetTitle(Form("sPHENIX primary idx=%zu truth clusters (x-y)", iprim));

                            TGraph *gr_truth_zr = new TGraph();
                            gr_truth_zr->SetName(Form("gr_sPHENIXPrimary_truthcluster_zr_idx%zu", iprim));
                            gr_truth_zr->SetTitle(Form("sPHENIX primary idx=%zu truth clusters (z-r)", iprim));

                            const size_t ntruthclus = std::min({truth_cl_x.size(), truth_cl_y.size(), truth_cl_z.size()});
                            if (ntruthclus == 0)
                                continue;

                            for (size_t ic = 0; ic < ntruthclus; ++ic)
                            {
                                const double x = truth_cl_x[ic];
                                const double y = truth_cl_y[ic];
                                const double z = truth_cl_z[ic];
                                const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);

                                gr_truth_xy->SetPoint(gr_truth_xy->GetN(), x, y);
                                gr_truth_zr->SetPoint(gr_truth_zr->GetN(), z, r_signed);
                            }

                            display.grPrimaryTruthClustersXY.push_back(gr_truth_xy);
                            display.grPrimaryTruthClustersZR.push_back(gr_truth_zr);
                        }
                    }

                    if (!IsData && sPHENIXPrimary_recocluster_globalX && sPHENIXPrimary_recocluster_globalY && sPHENIXPrimary_recocluster_globalZ)
                    {
                        const size_t nprimary_reco = std::min({sPHENIXPrimary_recocluster_globalX->size(), sPHENIXPrimary_recocluster_globalY->size(), sPHENIXPrimary_recocluster_globalZ->size()});
                        for (size_t iprim = 0; iprim < nprimary_reco; ++iprim)
                        {
                            if (!sPHENIXPrimary_charge || iprim >= sPHENIXPrimary_charge->size() || sPHENIXPrimary_charge->at(iprim) == 0.0)
                                continue;
                            const auto &reco_cl_x = sPHENIXPrimary_recocluster_globalX->at(iprim);
                            const auto &reco_cl_y = sPHENIXPrimary_recocluster_globalY->at(iprim);
                            const auto &reco_cl_z = sPHENIXPrimary_recocluster_globalZ->at(iprim);

                            TGraph *gr_reco_xy = new TGraph();
                            gr_reco_xy->SetName(Form("gr_sPHENIXPrimary_recocluster_xy_idx%zu", iprim));
                            gr_reco_xy->SetTitle(Form("sPHENIX primary idx=%zu reconstructed clusters (x-y)", iprim));

                            TGraph *gr_reco_zr = new TGraph();
                            gr_reco_zr->SetName(Form("gr_sPHENIXPrimary_recocluster_zr_idx%zu", iprim));
                            gr_reco_zr->SetTitle(Form("sPHENIX primary idx=%zu reconstructed clusters (z-r)", iprim));

                            const size_t nrecoclus = std::min({reco_cl_x.size(), reco_cl_y.size(), reco_cl_z.size()});
                            if (nrecoclus == 0)
                                continue;

                            for (size_t ic = 0; ic < nrecoclus; ++ic)
                            {
                                const double x = reco_cl_x[ic];
                                const double y = reco_cl_y[ic];
                                const double z = reco_cl_z[ic];
                                if (std::isnan(x) || std::isnan(y) || std::isnan(z))
                                    continue;

                                const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);

                                gr_reco_xy->SetPoint(gr_reco_xy->GetN(), x, y);
                                gr_reco_zr->SetPoint(gr_reco_zr->GetN(), z, r_signed);
                            }

                            if (gr_reco_xy->GetN() == 0 || gr_reco_zr->GetN() == 0)
                            {
                                delete gr_reco_xy;
                                delete gr_reco_zr;
                                continue;
                            }

                            display.grPrimaryRecoClustersXY.push_back(gr_reco_xy);
                            display.grPrimaryRecoClustersZR.push_back(gr_reco_zr);
                        }
                    }

                    out.crossingDisplays.push_back(display);
                    std::cout << "Prepared crossing display graphs for crossing key " << display.crossing_key << " (associated seeds = " << n_associated_silseeds << ", INTT clusters = " << n_intt_clusters_in_crossing << ")" << std::endl;
                }

                // Fill the same running histograms crossing-by-crossing.
                Tracklets_CombinatoricDR_FillEvent(tkldata, cfg, silicon_vtx, out, output_prefix);
                if (multiplicity_percentile_bin_index >= 0)
                {
                    Tracklets_CombinatoricDR_FillEvent(tkldata, cfg, silicon_vtx, multiplicity_percentile_outputs[multiplicity_percentile_bin_index], MultiplicityPercentileBinPrefix(output_prefix, multiplicity_percentile_bins[multiplicity_percentile_bin_index].percentile_bin));
                }
                ResetVec(tkldata);
            }
            else // if there are more than 1 vertices with this crossing, skip this event for simplicit
            {
                std::cout << "Skipping crossing " << pair.first << " with multiple vertices" << std::endl;
                continue;
            }
        }

        vertexCrossing_idx_map.clear();
        ResetVec(tkldata);
    }

    Tracklets_CombinatoricDR_Finalize(cfg, out, output_prefix);
    for (size_t i = 0; i < multiplicity_percentile_outputs.size(); ++i)
    {
        Tracklets_CombinatoricDR_Finalize(cfg, multiplicity_percentile_outputs[i], MultiplicityPercentileBinPrefix(output_prefix, multiplicity_percentile_bins[i].percentile_bin));
    }

    const bool write_ok = WriteTracklets_CombinatoricDR_Output(out, multiplicity_percentile_outputs, outfilename.Data());
    if (!write_ok)
    {
        std::cerr << "Failed to create output ROOT file: " << outfilename << std::endl;
        f->Close();
        delete f;
        return 1;
    }

    std::cout << "Wrote combinatoric histograms to " << outfilename << std::endl;
    ResetVec(tkldata);
    f->Close();
    delete f;

    return 0;
}
