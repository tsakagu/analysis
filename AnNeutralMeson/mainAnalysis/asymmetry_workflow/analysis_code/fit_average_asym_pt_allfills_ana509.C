#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

void fit_average_asym_pt_allfills_ana509(const std::string &production_date = "01312026",
                                         const std::string &selection_label = "MBD_INCOMPLETE",
                                         const std::string &inputfolder_polarization = "",
                                         const std::string &folder_asymmetry = "")
{
  const int nPtBins = 9;
  const float pTBins[nPtBins][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 10}, {10, 20}};
  const int nBeams = 2; // Yellow or blue
  const std::string beams[nBeams] = {"yellow", "blue"};
  const int nParticles = 2; // pi0 or eta
  const std::string particle[nParticles] = {"pi0", "eta"};
  const int nRegions = 2; // peak band or side_band
  const std::string regions[nRegions] = {"peak", "side"};
  //const std::string regions2[nRegions] = {"", "bkg"};
  const int nSpins = 2; // up or down
  const std::string spins[nSpins] = {"up", "down"};
  const int nSides = 2; // left or right
  const std::string sides[nSides] = {"left", "right"};
  const double etaThreshold = 0.35;
  const int nEtaRegions = 2;
  const std::string etaRegions[nEtaRegions] = {"low", "high"};
  static constexpr int nDirections = 2;
  const std::string directions = {"forward", "backward"};
  const int nMethods = 2; // 0: square root asymmetry   1: relative luminosity
  std::string methods[2] = {"sqrt", "rellum"};

  SetsPhenixStyle();

  const std::string inputfilename_template = folder_asymmetry + "/analysis_fills_asymmetry_nodir_complete_ana509_" + production_date + "_" + selection_label + "/analysis_";
  const std::string outputfilename = folder_asymmetry + "/analysis_asymmetry_nodir_complete_ana509_" + production_date + "_" + selection_label + "_averaged.root";

  TFile *outputfile = new TFile(outputfilename.c_str(), "RECREATE");
  outputfile->cd();

  TGraphErrors *graph_asym_average[nBeams][nParticles][nRegions][nPtBins][nMethods] = {0};
  TGraphErrors *graph_asym_fills[nBeams][nParticles][nRegions][nPtBins][nMethods][12] = {nullptr};
  TH1F *h_asym[nBeams][nParticles][nRegions][nPtBins][nMethods][12] = {nullptr};
  TH1F *h_asym_weighted[nBeams][nParticles][nRegions][nPtBins][nMethods][12] = {nullptr};

  float sum_weights[nBeams][nParticles][nRegions][nPtBins][nMethods][12] = {0};
  float sum_weighted_asymmetries[nBeams][nParticles][nRegions][nPtBins][nMethods][12] = {0};
  float sum_weighted_squares[nBeams][nParticles][nRegions][nPtBins][nMethods][12] = {0};

  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iR = 0; iR < nRegions; iR++)
    {
      for (int iP = 0; iP < nParticles; iP++)
      {
        for (int iPt = 0; iPt < nPtBins; iPt++)
        {
          for (int imethod = 0; imethod < nMethods; imethod++)
          {
            std::stringstream graph_AN_name;
            std::stringstream graph_AN_title;
            graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_pT_" << iPt << "_" << methods[imethod];
            graph_AN_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.")  << ", " << beams[iB] << " beam); #phi; " << (iR == 0 ? "Raw Asymmetry" : "Bkg Asymmetry");
            graph_asym_average[iB][iP][iR][iPt][imethod] = new TGraphErrors();
            graph_asym_average[iB][iP][iR][iPt][imethod]->SetName(graph_AN_name.str().c_str());
            graph_asym_average[iB][iP][iR][iPt][imethod]->SetTitle(graph_AN_title.str().c_str());
            for (int iphi = 0; iphi < 12; iphi++)
            {
              std::stringstream h_name;
              std::stringstream h_title;
              h_name << "h_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_pT_" << iPt << "_" << methods[imethod] << "_phi_" << iphi;
              h_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.") << ", " << beams[iB] << " beam) i_{#phi} = " << iphi << "; A_N; Counts";
              h_asym[iB][iP][iR][iPt][imethod][iphi] = new TH1F(h_name.str().c_str(), h_title.str().c_str(), 100, -1, 1);
              h_name << "_weighted";
              h_title.str("");
              h_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.") << ", " << beams[iB] << " beam) i_{#phi} = " << iphi << "; A_N; Weighted Counts";
              h_asym_weighted[iB][iP][iR][iPt][imethod][iphi] = new TH1F(h_name.str().c_str(), h_title.str().c_str(), 100, -1, 1);

              std::stringstream graph_fill_name;
              std::stringstream graph_fill_title;
              graph_fill_name << "graph_fill_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_pT_" << iPt << "_" << methods[imethod] << "_phi_" << iphi;
              graph_fill_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c  (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.")  << ", " << beams[iB] << " beam) i_{#phi} = " << iphi << "; Fill; A_N";
              graph_asym_fills[iB][iP][iR][iPt][imethod][iphi] = new TGraphErrors();
              graph_asym_fills[iB][iP][iR][iPt][imethod][iphi]->SetName(graph_fill_name.str().c_str());
              graph_asym_fills[iB][iP][iR][iPt][imethod][iphi]->SetTitle(graph_fill_title.str().c_str());
            }
          }
        }
      }
    }
  }

  TFile* f = TFile::Open((inputfolder_polarization + "/luminosity_fills_ana509.root").c_str());
  TGraph* graph_relative_luminosity_fill_blue_mbdns = (TGraph*)f->Get("graph_relative_luminosity_fill_blue_mbdns");
  int n = graph_relative_luminosity_fill_blue_mbdns->GetN();
  for (int i = 0; i < n; ++i)
  {
    int fillnumber = graph_relative_luminosity_fill_blue_mbdns->GetPointX(i);
    std::cout << "fill " << fillnumber << std::endl;
    std::string inputfilename = inputfilename_template + std::to_string(fillnumber) + ".root";
    TFile *inputfile = TFile::Open(inputfilename.c_str());

    if (!inputfile) {
      std::cout << "skip fill " << fillnumber << std::endl;
      continue;
    }

    for (int iB = 0; iB < nBeams; iB++)
    {
      for (int iR = 0; iR < nRegions; iR++)
      {
        for (int iP = 0; iP < nParticles; iP++)
        {
          for (int iPt = 0; iPt < nPtBins; iPt++)
          {
            for (int imethod = 0; imethod < nMethods; imethod++)
            {
              std::stringstream graph_AN_name;
              graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_pT_" << iPt << "_" << methods[imethod];
              TGraphErrors *graph_asym = (TGraphErrors*) inputfile->Get(graph_AN_name.str().c_str());
              for (int iphi = 0; iphi < 12; iphi++)
              {
                float asym = graph_asym->GetPointY(iphi);
                float asym_err = graph_asym->GetErrorY(iphi);
                if (asym == asym && asym_err == asym_err && asym_err > 0)
                {
                  h_asym[iB][iP][iR][iPt][imethod][iphi]->Fill(asym);
                  h_asym_weighted[iB][iP][iR][iPt][imethod][iphi]->Fill(asym, std::pow(asym_err, -2));
                  graph_asym_fills[iB][iP][iR][iPt][imethod][iphi]->SetPoint(i, fillnumber, asym);
                  graph_asym_fills[iB][iP][iR][iPt][imethod][iphi]->SetPointError(i, 0, asym_err);
                  sum_weights[iB][iP][iR][iPt][imethod][iphi] += std::pow(asym_err, -2);
                  sum_weighted_asymmetries[iB][iP][iR][iPt][imethod][iphi] += asym * std::pow(asym_err, -2);
                }
                else
                {
                  /*std::cout << "SCREAM empty bin :"
                            << "iB = " << iB << ", "
                            << "iP = " << iP << ", "
                            << "iR = " << iR << ", "
                            << "iPt = " << iPt << ", "
                            << "imethod = " << imethod << ", "
                            << "iphi = " << iphi << ", "
                            << "fill = " << i << std::endl;
                  */
                }
              }
            }
          }
        }
      }
    }
    inputfile->Close();
  }
  
  outputfile->cd();
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iR = 0; iR < nRegions; iR++)
    {
      for (int iP = 0; iP < nParticles; iP++)
      {
        for (int iPt = 0; iPt < nPtBins; iPt++)
        {
          for (int imethod = 0; imethod < nMethods; imethod++)
          {
            for (int iphi = 0; iphi < 12; iphi++)
            {
              float asym = sum_weighted_asymmetries[iB][iP][iR][iPt][imethod][iphi] / sum_weights[iB][iP][iR][iPt][imethod][iphi];
              float asym_err = std::sqrt(1. / sum_weights[iB][iP][iR][iPt][imethod][iphi]);
              if (asym == asym && asym_err == asym_err) {
                graph_asym_average[iB][iP][iR][iPt][imethod]->SetPoint(iphi,
                                                                             (float) iphi * (2 * M_PI / 12.) - (M_PI - M_PI / 12.),
                                                                             asym);
                graph_asym_average[iB][iP][iR][iPt][imethod]->SetPointError(iphi, 0, asym_err);
                graph_asym_fills[iB][iP][iR][iPt][imethod][iphi]->Write();
              }
              else { // NaN
                graph_asym_average[iB][iP][iR][iPt][imethod]->SetPoint(iphi,
                                                                             (float) iphi * (2 * M_PI / 12.) - (M_PI - M_PI / 12.),
                                                                             0);
                graph_asym_average[iB][iP][iR][iPt][imethod]->SetPointError(iphi, 0, 0);
                graph_asym_fills[iB][iP][iR][iPt][imethod][iphi]->Write();
              }
            }
            graph_asym_average[iB][iP][iR][iPt][imethod]->Write();
          }
        }
      }
    }
  }

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();
  delete outputfile;

  gSystem->Exit(0);
}
