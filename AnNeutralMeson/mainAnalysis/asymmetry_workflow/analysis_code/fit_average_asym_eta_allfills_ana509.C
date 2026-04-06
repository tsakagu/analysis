#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

void fit_average_asym_eta_allfills_ana509(const std::string &production_date = "01312026",
                                          const std::string &selection_label = "MBD_INCOMPLETE",
                                          const std::string &inputfolder_polarization = "",
                                          const std::string &folder_asymmetry = "")
{
  const int nEtaBins = 8;
  const float etaBins[nEtaBins][2] = {{-2.00, -1.05}, {-1.05, -0.86}, {-0.86, -0.61}, {-0.61, 0.0}, {0.0, 0.61}, {0.61, 0.86}, {0.86, 1.05}, {1.05, 2.0}};
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
  const int nMethods = 2; // 0: square root asymmetry   1: relative luminosity
  std::string methods[2] = {"sqrt", "rellum"};

  SetsPhenixStyle();

  const std::string inputfilename_template = folder_asymmetry + "/analysis_fills_asymmetry_eta_complete_ana509_" + production_date + "_" + selection_label + "/analysis_";
  const std::string outputfilename = folder_asymmetry + "/analysis_asymmetry_eta_complete_ana509_" + production_date + "_" + selection_label + "_averaged.root";

  TFile *outputfile = new TFile(outputfilename.c_str(), "RECREATE");
  outputfile->cd();

  TGraphErrors *graph_asym_average[nBeams][nParticles][nRegions][nEtaBins][nMethods] = {0};
  TGraphErrors *graph_asym_fills[nBeams][nParticles][nRegions][nEtaBins][nMethods][12] = {nullptr};
  TH1F *h_asym[nBeams][nParticles][nRegions][nEtaBins][nMethods][12] = {nullptr};
  TH1F *h_asym_weighted[nBeams][nParticles][nRegions][nEtaBins][nMethods][12] = {nullptr};

  float sum_weights[nBeams][nParticles][nRegions][nEtaBins][nMethods][12] = {0};
  float sum_weighted_asymmetries[nBeams][nParticles][nRegions][nEtaBins][nMethods][12] = {0};
  float sum_weighted_squares[nBeams][nParticles][nRegions][nEtaBins][nMethods][12] = {0};

  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iR = 0; iR < nRegions; iR++)
    {
      for (int iP = 0; iP < nParticles; iP++)
      {
        for (int ieta = 0; ieta < nEtaBins; ieta++)
        {
          for (int imethod = 0; imethod < nMethods; imethod++)
          {
            std::stringstream graph_AN_name;
            std::stringstream graph_AN_title;
            graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_eta_" << ieta << "_" << methods[imethod];
            graph_AN_title << "#eta #in [" << std::fixed << std::setprecision(1) << etaBins[ieta][0] << ", " << etaBins[ieta][1] << "] rad (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.")  << ", " << beams[iB] << " beam); #phi; " << (iR == 0 ? "Raw Asymmetry" : "Bkg Asymmetry");
            graph_asym_average[iB][iP][iR][ieta][imethod] = new TGraphErrors();
            graph_asym_average[iB][iP][iR][ieta][imethod]->SetName(graph_AN_name.str().c_str());
            graph_asym_average[iB][iP][iR][ieta][imethod]->SetTitle(graph_AN_title.str().c_str());
            for (int iphi = 0; iphi < 12; iphi++)
            {
              std::stringstream h_name;
              std::stringstream h_title;
              h_name << "h_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_eta_" << ieta << "_" << methods[imethod] << "_phi_" << iphi;
              h_title << "#eta #in [" << std::fixed << std::setprecision(1) << etaBins[ieta][0] << ", " << etaBins[ieta][1] << "] (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.") << ", " << beams[iB] << " beam) i_{#phi} = " << iphi << "; A_N; Counts";
              h_asym[iB][iP][iR][ieta][imethod][iphi] = new TH1F(h_name.str().c_str(), h_title.str().c_str(), 100, -1, 1);
              h_name << "_weighted";
              h_title.str("");
              h_title << "#eta #in [" << std::fixed << std::setprecision(1) <<  etaBins[ieta][0] << ", " <<  etaBins[ieta][1] << "]  (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.") << ", " << beams[iB] << " beam) i_{#phi} = " << iphi << "; A_N; Weighted Counts";
              h_asym_weighted[iB][iP][iR][ieta][imethod][iphi] = new TH1F(h_name.str().c_str(), h_title.str().c_str(), 100, -1, 1);
              
              std::stringstream graph_fill_name;
              std::stringstream graph_fill_title;
              graph_fill_name << "graph_fill_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_eta_" << ieta << "_" << methods[imethod] << "_phi_" << iphi;
              graph_fill_title << "#eta #in [" << std::fixed << std::setprecision(1) << etaBins[ieta][0] << ", " << etaBins[ieta][1] << "] rad  (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.")  << ", " << beams[iB] << " beam) i_{#phi} = " << iphi << "; Fill; A_N";
              graph_asym_fills[iB][iP][iR][ieta][imethod][iphi] = new TGraphErrors();
              graph_asym_fills[iB][iP][iR][ieta][imethod][iphi]->SetName(graph_fill_name.str().c_str());
              graph_asym_fills[iB][iP][iR][ieta][imethod][iphi]->SetTitle(graph_fill_title.str().c_str());
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
          for (int ieta = 0; ieta < nEtaBins; ieta++)
          {
            for (int imethod = 0; imethod < nMethods; imethod++)
            {
              std::stringstream graph_AN_name;
              graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_eta_" << ieta << "_" << methods[imethod];
              //std::cout << inputfile << " and " << graph_AN_name.str() << std::endl;
              TGraphErrors *graph_asym = (TGraphErrors*) inputfile->Get(graph_AN_name.str().c_str());
              for (int iphi = 0; iphi < 12; iphi++)
              {
                float asym = graph_asym->GetPointY(iphi);
                float asym_err = graph_asym->GetErrorY(iphi);
                if (asym == asym && asym_err == asym_err && asym_err > 0)
                {
                  h_asym[iB][iP][iR][ieta][imethod][iphi]->Fill(asym);
                  h_asym_weighted[iB][iP][iR][ieta][imethod][iphi]->Fill(asym, std::pow(asym_err, -2));
                  graph_asym_fills[iB][iP][iR][ieta][imethod][iphi]->SetPoint(i, fillnumber, asym);
                  graph_asym_fills[iB][iP][iR][ieta][imethod][iphi]->SetPointError(i, 0, asym_err);
                  sum_weights[iB][iP][iR][ieta][imethod][iphi] += std::pow(asym_err, -2);
                  sum_weighted_asymmetries[iB][iP][iR][ieta][imethod][iphi] += asym * std::pow(asym_err, -2);
                }
              }
            }
          }
        }
      }
    }
    inputfile->Close();
  }

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
          for (int ieta = 0; ieta < nEtaBins; ieta++)
          {
            for (int imethod = 0; imethod < nMethods; imethod++)
            {
              std::stringstream graph_AN_name;
              graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_eta_" << ieta << "_" << methods[imethod];
              TGraphErrors *graph_asym = (TGraphErrors*) inputfile->Get(graph_AN_name.str().c_str());
              for (int iphi = 0; iphi < 12; iphi++)
              {
                float asym = graph_asym->GetPointY(iphi);
                float asym_err = graph_asym->GetErrorY(iphi);
                if (asym == asym && asym_err == asym_err && asym_err != 0)
                {
                  float mean = sum_weighted_asymmetries[iB][iP][iR][ieta][imethod][iphi] / sum_weights[iB][iP][iR][ieta][imethod][iphi];
                  sum_weighted_squares[iB][iP][iR][ieta][imethod][iphi] += (asym - mean) * (asym - mean) * std::pow(asym_err, -2);
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
        for (int ieta = 0; ieta < nEtaBins; ieta++)
        {
          for (int imethod = 0; imethod < nMethods; imethod++)
          {
            for (int iphi = 0; iphi < 12; iphi++)
            {
              float asym = sum_weighted_asymmetries[iB][iP][iR][ieta][imethod][iphi] / sum_weights[iB][iP][iR][ieta][imethod][iphi];
              //float asym_err = std::sqrt(sum_weighted_squares[iB][iP][iR][ieta][imethod][iphi] / sum_weights[iB][iP][iR][ieta][imethod][iphi]);
              float asym_err = std::sqrt(1. / sum_weights[iB][iP][iR][ieta][imethod][iphi]);
              if (asym == asym && asym_err) {
                graph_asym_average[iB][iP][iR][ieta][imethod]->SetPoint(iphi,
                                                                        (float) iphi * (2 * M_PI / 12.) - (M_PI - M_PI / 12.),
                                                                        asym);
                graph_asym_average[iB][iP][iR][ieta][imethod]->SetPointError(iphi, 0, asym_err);
                graph_asym_fills[iB][iP][iR][ieta][imethod][iphi]->Write();
              }
              else { // NaN
                graph_asym_average[iB][iP][iR][ieta][imethod]->SetPoint(iphi,
                                                                        (float) iphi * (2 * M_PI / 12.) - (M_PI - M_PI / 12.),
                                                                        0);
                graph_asym_average[iB][iP][iR][ieta][imethod]->SetPointError(iphi, 0, 0);
                graph_asym_fills[iB][iP][iR][ieta][imethod][iphi]->Write();
              }
            }
            graph_asym_average[iB][iP][iR][ieta][imethod]->Write();
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
