#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

TGraphErrors* SqrtGraph(TH1* phi_up, TH1* phi_down, int iPt, double Pol);
double SqrtAsym(double NLup, double NLdown, double NRup, double NRdown);
double SqrtError(double NLup, double NLdown, double NRup, double NRdown,
                 double NLupErr, double NLdownErr, double NRupErr, double NRdownErr);

TGraphErrors* RelLumGraph(TH1* phi_up, TH1* phi_down, double RL, double Pol, bool print = false);
double RelLumAsym(double Nup, double Ndown, double RL);
double RelLumError(double Nup, double Ndown, double NupErr, double NdownErr, double RL);

void CreatePlotTGraphErrors(const std::string& plots_folder, TGraphErrors *graph_AN);

void fit_asym_pt_ana509(const std::string &production_date,
                        const std::string &selection_label,
                        const std::string &inputfolder_mass,
                        const std::string &inputfolder_analysis,
                        const std::string &outputfolder_asymmetry,
                        double RL_blue, double RL_yellow, double P_blue, double P_yellow,
                        int inputfill = 0);

bool printcanvas = false;

void fit_asym_pt_all_fills_ana509(const std::string &production_date = "01312026",
                                  const std::string &selection_label = "MBD_pt_l4_0mrad",
                                  bool separate_fills = true,
                                  const std::string &inputfolder_mass = "",
                                  const std::string &inputfolder_polarization = "",
                                  const std::string &inputfolder_analysis = "",
                                  const std::string &outputfolder_asymmetry = "")
{
  if (separate_fills)
  {
    // Open the file and get the TGraph
    TFile* f = TFile::Open((inputfolder_polarization + "/luminosity_fills_ana509.root").c_str());
    TGraph* graph_relative_luminosity_fill_blue_mbdns = (TGraph*)f->Get("graph_relative_luminosity_fill_blue_mbdns");
    TGraph* graph_relative_luminosity_fill_yellow_mbdns = (TGraph*)f->Get("graph_relative_luminosity_fill_yellow_mbdns");
    TGraph* graph_polarization_fill_blue_mbdns = (TGraph*)f->Get("graph_polarization_fill_blue_mbdns");
    TGraph* graph_polarization_fill_yellow_mbdns = (TGraph*)f->Get("graph_polarization_fill_yellow_mbdns");
    
    int n = graph_relative_luminosity_fill_blue_mbdns->GetN();
    for (int i = 0; i < n; ++i) {
      double fillnumber, RL_blue, RL_yellow, P_blue, P_yellow;
      graph_relative_luminosity_fill_blue_mbdns->GetPoint(i, fillnumber, RL_blue);
      graph_relative_luminosity_fill_yellow_mbdns->GetPoint(i, fillnumber, RL_yellow);
      graph_polarization_fill_blue_mbdns->GetPoint(i, fillnumber, P_blue);
      graph_polarization_fill_yellow_mbdns->GetPoint(i, fillnumber, P_yellow);
      std::cout << "fill: " << fillnumber << " (" << RL_blue << "," << RL_yellow << "," << P_blue << "," << P_yellow << ")" << std::endl;
      fit_asym_pt_ana509(production_date, selection_label, inputfolder_mass, inputfolder_analysis, outputfolder_asymmetry, RL_blue, RL_yellow, P_blue, P_yellow, fillnumber);
    }
    f->Close();
  }
  else
  {
    // Average values
    fit_asym_pt_ana509(production_date, selection_label, inputfolder_mass, inputfolder_analysis, outputfolder_asymmetry, 1.0012080359, 1.0013263725, 0.5082929939, 0.5167388276);
  }
  gSystem->Exit(0);
}

void fit_asym_pt_ana509(const std::string &production_date,
                        const std::string &selection_label,
                        const std::string &inputfolder_mass,
                        const std::string &inputfolder_analysis,
                        const std::string &outputfolder_asymmetry,
                        double RL_blue, double RL_yellow, double P_blue, double P_yellow,
                        int inputfill)
{
  std::string inputfilename = "";
  std::string outputfilename = "";

  if (inputfill == 0)
  {
    inputfilename = inputfolder_analysis + "/analysis_complete_ana509_" + production_date + "_" + selection_label + ".root";
    std::string outputfolder = outputfolder_asymmetry;
    gSystem->Exec(("mkdir -p " + outputfolder).c_str());
    outputfilename = outputfolder + "/analysis_asymmetry_complete_ana509_" + production_date + "_" + selection_label + ".root";
  }
  else
  {
    inputfilename = inputfolder_analysis + "/analysis_fills_complete_ana509_" + production_date + "_" + selection_label + "/analysis_" + std::to_string(inputfill) + ".root";
    std::string outputfolder = outputfolder_asymmetry + "/analysis_fills_asymmetry_nodir_complete_ana509_" + production_date + "_" + selection_label;
    gSystem->Exec(("mkdir -p " + outputfolder).c_str());
    outputfilename = outputfolder + "/analysis_" + std::to_string(inputfill) + ".root";
  }
  
  const int nPtBins = 9;
  const float pTBins[nPtBins][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 10}, {10, 20}};
    
  // Read input file
  TFile *inputfile = TFile::Open(inputfilename.c_str());

  if (!inputfile) return;
  
  const int nBeams = 2; // Yellow or blue
  const std::string beams[nBeams] = {"yellow", "blue"};
  const int nParticles = 2; // pi0 or eta
  const std::string particle[nParticles] = {"pi0", "eta"};
  const int nRegions = 2; // peak band or side_band
  const std::string regions[nRegions] = {"peak", "side"};
  const std::string regions2[nRegions] = {"", "bkg"};
  const int nSpins = 2; // up or down
  const std::string spins[nSpins] = {"up", "down"};
  const int nSides = 2; // left or right
  const std::string sides[nSides] = {"left", "right"};
  const double etaThreshold = 0.35;
  const int nEtaRegions = 2;
  const std::string etaRegions[nEtaRegions] = {"low", "high"};
  static constexpr int nDirections = 2;
  // (eta * beamDirection) > 0 or (eta * beamDirection < 0) -> forward or backward
  const std::string directions[nDirections] = {"forward", "backward"};
  TH1F *h_yield[nBeams][nParticles][nRegions][nPtBins][nSpins];
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iR = 0; iR < nRegions; iR++)
      {
        for (int iS = 0; iS < nSpins; iS++)
        {
          for (int iPt = 0; iPt < nPtBins; iPt++)
          {
            TH1F *h_yield_dir[2] = {nullptr};
            for (int iDir = 0; iDir < nDirections; iDir++)
            {
              std::stringstream h_yield_name;
              h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_"
                           << "pT_" << iPt << "_"
                           << directions[iDir] << "_"
                           << spins[iS];
              h_yield_dir[iDir] = (TH1F*)inputfile->Get(h_yield_name.str().c_str());
            }
            std::stringstream h_yield_name;
            h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_"
                         << "pT_" << iPt << "_" << spins[iS];
            const char *h_yield_title = h_yield_dir[0]->GetTitle();
            const int nBins = h_yield_dir[0]->GetXaxis()->GetNbins();
            const int xMin = h_yield_dir[0]->GetXaxis()->GetXmin();
            const int xMax = h_yield_dir[0]->GetXaxis()->GetXmax();
            h_yield[iB][iP][iR][iPt][iS] = new TH1F(h_yield_name.str().c_str(), h_yield_title, nBins, xMin, xMax);
            for (int iPhi = 1; iPhi <= 12; iPhi++)
            {
              h_yield[iB][iP][iR][iPt][iS]->SetBinContent(iPhi,
                                                          h_yield_dir[0]->GetBinContent(iPhi) + h_yield_dir[1]->GetBinContent(iPhi));
            }
          }
        }
      }
    }
  }
  
  std::string plots_folder = "plots_asymmetry/ana509_p022/analysis_" + production_date + "/pt_asymmetries/phi_asymmetries/" + selection_label + "/";
  gSystem->Exec(("mkdir -p " + plots_folder).c_str());

  TFile *outputfile = new TFile(outputfilename.c_str(), "RECREATE");
  
  const int nMethods = 2; // 0: square root asymmetry   1: relative luminosity
  std::string methods[2] = {"sqrt", "rellum"};

  outputfile->cd();
  // Compute the asymmetry
  for (int iB = 0; iB < nBeams; iB++)
  {
    double RL = (iB == 0 ? RL_yellow : RL_blue);

    // polarization
    double P = (iB == 0 ? P_yellow : P_blue);
    
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iPt = 0; iPt < nPtBins; iPt++)
      {
        for (int imethod = 0; imethod < nMethods; imethod++)
        {
          // raw asymmetry
          TH1F *phi_up = h_yield[iB][iP][0][iPt][0];
          TH1F *phi_down = h_yield[iB][iP][0][iPt][1];
          TGraphErrors *graph_AN = (imethod == 0 ? SqrtGraph(phi_up, phi_down, iPt, P) : RelLumGraph(phi_up, phi_down, RL, P));
          
          std::string h_name = h_yield[iB][iP][0][iPt][0]->GetName();
          
          std::stringstream graph_AN_name;
          graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_peak_pT_" << iPt << "_" << methods[imethod];
          graph_AN->SetName(graph_AN_name.str().c_str());
          
          std::string canvas_name = "canvas_" + h_name.substr(h_name.find("h_yield_") + 8) + "_" + methods[imethod];
          std::stringstream graph_title;
          graph_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c; #phi; raw asymmetry";
          graph_AN->SetTitle(graph_title.str().c_str());
          graph_AN->Write();
          delete graph_AN; // deleting the file should be enough
          
          // bkg asymmetry
          phi_up = h_yield[iB][iP][1][iPt][0];
          phi_down = h_yield[iB][iP][1][iPt][1];
          graph_AN = (imethod == 0 ? SqrtGraph(phi_up, phi_down, iPt, P) : RelLumGraph(phi_up, phi_down, RL, P));
          if (graph_AN->GetErrorY(1) == graph_AN->GetErrorY(1) && graph_AN->GetErrorY(1) < 0)
          {
            std::cout << "weird error: " << graph_AN->GetErrorY(1) << " at indices (" << iB << "," << iP << "," << 1 << "," << iPt << ")" << std::endl;
          }
            
          h_name = h_yield[iB][iP][1][iPt][0]->GetName();
          
          graph_AN_name.str("");
          graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_side_pT_" << iPt << "_" << methods[imethod];
          graph_AN->SetName(graph_AN_name.str().c_str());
          
          canvas_name = "canvas_" + h_name.substr(h_name.find("h_yield_") + 8) + "_" + methods[imethod];
          graph_title.str("");
          graph_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c; #phi; bkg asymmetry";
          graph_AN->SetTitle(graph_title.str().c_str());
          graph_AN->Write();
          delete graph_AN;
        }
      }
    }
  }
  
  inputfile->Close();
}

TGraphErrors* SqrtGraph(TH1* phi_up, TH1* phi_down, int iPt, double Pol)
{
  TGraphErrors* gr = new TGraphErrors();
  int nbins = phi_up->GetNbinsX();
  for (int i = 0; i < nbins; i++)
  {
    float phi = i*(2*M_PI/nbins) - (M_PI - M_PI/nbins);
    int phibinL = i; // N_Left
    int phibinR = (phibinL + (int)(nbins/2.)) % nbins; // N_Right
      
    double NLup = phi_up->GetBinContent(1+phibinL);
    double NLdown = phi_down->GetBinContent(1+phibinL);
    double NRup = phi_up->GetBinContent(1+phibinR);
    double NRdown = phi_down->GetBinContent(1+phibinR);
    double NLupErr = phi_up->GetBinError(1+phibinL); // Simply the square root of NLup
    double NLdownErr = phi_down->GetBinError(1+phibinL);
    double NRupErr = phi_up->GetBinError(1+phibinR);
    double NRdownErr = phi_down->GetBinError(1+phibinR);

    double asym = 1 * SqrtAsym(NLup, NLdown, NRup, NRdown) / Pol;
    double err = 1 * SqrtError(NLup, NLdown, NRup, NRdown, NLupErr, NLdownErr, NRupErr, NRdownErr) / Pol;

    // if (NLup <= 0 || NLdown <= 0 || NRup <= 0 || NRdown <= 0) {
    //   asym = std::numeric_limits<double>::quiet_NaN();
    //   err = std::numeric_limits<double>::quiet_NaN();
    // }
    
    //double err = std::sqrt(NLup + NLdown + NRup + NRdown) / Pol;

    /*if (err < 0)
    {
      std::cout << "(NLup, NLdown, NRup, NRdown, NLupErr, NLdownErr, NRupErr, NRdownErr, err) = ("
                << NLup << "," << NLdown << "," << NRup << "," << NRdown << "," << NLupErr << "," << NLdownErr << "," << NRupErr << "," << NRdownErr << "," << err << ")" << std::endl;
    }
    */
    
    // if ((asym == asym) && (err == err)) // i.e. asym and err are not NaN
    // {
    gr->SetPoint(i, phi, asym);
    gr->SetPointError(i, 0, err);
    // }
  }
  return gr;
}

TGraphErrors* RelLumGraph(TH1* phi_up, TH1* phi_down, double RL, double Pol, bool print) {
  TGraphErrors* gr = new TGraphErrors();
  int nbins = phi_up->GetNbinsX();
  for (int i = 0; i < nbins; i++)
  {
    float phi = i * (2*M_PI/nbins) - (M_PI - M_PI/nbins);
    int phibin = i;
    double Nup = phi_up->GetBinContent(1+phibin);
    double Ndown = phi_down->GetBinContent(1+phibin);
    double NupErr = phi_up->GetBinError(1+phibin);
    double NdownErr = phi_down->GetBinError(1+phibin);

    double asym = RelLumAsym(Nup, Ndown, RL) / Pol;
    double err = RelLumError(Nup, Ndown, NupErr, NdownErr, RL) / Pol;
    //double err = std::sqrt(Nup + Ndown) / Pol;

    // if (Nup <= 0 || Ndown <= 0) {
    //   asym = std::numeric_limits<double>::quiet_NaN();
    //   err = std::numeric_limits<double>::quiet_NaN();
    // }

    if (print)
    {
      std::cout << i << ":" << phi << "," << Pol << "," << RL << "," << Nup << "," << Ndown << "," << NupErr << "," << NdownErr << "," << asym << "," << err << std::endl;
    }

    /*if (err < 0)
    {
      std::cout << "(Nup, Ndown, NupErr, NdownErr, RL, Pol, err) = ("
                << Nup << "," << Ndown << "," << NupErr << "," << NdownErr << "," << RL << "," << Pol << "," << err << ")" << std::endl;
    }
    */
    
    // if ((asym == asym) && (err == err)) // i.e. asym and err are not NaN
    // { 
    gr->SetPoint(i, phi, asym);
    gr->SetPointError(i, 0, err);
    // }
  }
  return gr;
}


double SqrtAsym(double NLup, double NLdown, double NRup, double NRdown) {
  double numerator = sqrt(NLup*NRdown)-sqrt(NRup*NLdown);
  double denominator = sqrt(NLup*NRdown)+sqrt(NRup*NLdown);
  double asym = (denominator != 0 ? numerator / denominator : 1000);
  return asym;
}

double SqrtError(double NLup, double NLdown, double NRup, double NRdown,
                 double NLupErr, double NLdownErr, double NRupErr, double NRdownErr) {
  double t1 = sqrt(NLdown*NRup*NRdown/NLup)*NLupErr;
  double t2 = sqrt(NLup*NRup*NRdown/NLdown)*NLdownErr;
  double t3 = sqrt(NLup*NLdown*NRdown/NRup)*NRupErr;
  double t4 = sqrt(NLup*NLdown*NRup/NRdown)*NRdownErr;
  double denominator = pow(sqrt(NLup*NRdown)+sqrt(NRup*NLdown),2);
  double asym_err = sqrt(pow(t1,2)+pow(t2,2)+pow(t3,2)+pow(t4,2)) / denominator;
  return asym_err;
}

// Compute the asymmetry error with respect to the relative luminosity definition:
double RelLumAsym(double Nup, double Ndown, double RL)
{
  double numerator = Nup - RL * Ndown;
  double denominator = Nup + RL * Ndown;
  double  asym = (denominator != 0 ? numerator / denominator : 1000);
  return asym;
}

// Compute the asymmetry error with respect to the relative luminosity definition:
double RelLumError(double Nup, double Ndown, double NupErr, double NdownErr, double RL)
{
  // hypothesis: negligible uncertainty on the relative luminosity compared to the yields
  double t1 = 2 * RL * Ndown * NupErr;
  double t2 = 2 * RL * Nup * NdownErr;
  double denominator = pow(Nup + RL * Ndown,2);
  double asym_err = (denominator != 0 ? sqrt(pow(t1, 2) + pow(t2, 2)) / denominator : 0);
  return asym_err;
}

void CreatePlotTGraphErrors(const std::string& plots_folder, TGraphErrors *graph_AN)
{
  const char *characters = graph_AN->GetName();
  std::string graph_name(characters);
  std::string extension = graph_name.substr(graph_name.find_first_of("_")+1);
  std::string canvas_name = "canvas_" + extension;
  TCanvas *canvas = nullptr;
  if (printcanvas)
  {
    canvas = new TCanvas(canvas_name.c_str());
    canvas->cd();
    graph_AN->SetStats(0);
    //graph_AN->SetMinimum(-0.05);
    //graph_AN->SetMaximum(0.05);
    graph_AN->SetMarkerStyle(kPlus);
    graph_AN->SetLineStyle(0);
    graph_AN->Draw("ALP");
    canvas->SaveAs((plots_folder + "/" + canvas_name + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name + ".png").c_str());
    delete canvas;
  }
}
