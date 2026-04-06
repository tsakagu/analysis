#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

bool printcanvas = false;

void PlotOneAsym(std::string type, std::string plots_folder, std::string canvas_name, std::string graph_title, TGraphErrors* gr, float& epsilon, float& depsilon, float& offset, float& doffset);

void fit_sinusoid_asym_pt_ana509(const std::string &production_date = "01312026",
                                 const std::string &selection_label = "MBD_INCOMPLETE",
                                 const std::string &folder_asymmetry = "",
                                 const std::string &folder_csv_inputs = "",
                                 const std::string &folder_csv_ratios = "",
                                 bool gpr_ratio = false)
{
  std::string ratio_suffix = (gpr_ratio ? "_gpr" : "");
  std::string inputfilename = folder_asymmetry + "/analysis_asymmetry_nodir_complete_ana509_" + production_date + "_" + selection_label + "_averaged.root";
  std::string outputfilename = folder_asymmetry + "/analysis_asymmetry_nodir_complete_ana509_" + production_date + "_" + selection_label + "_averaged_fitted" + ratio_suffix + ".root";

  std::string plots_folder = "plots_asymmetry/ana509_p022/analysis_" + production_date + "/pt_asymmetries/" + selection_label + "/phi_asymmetries/runs";
  gSystem->Exec(("mkdir -p " + plots_folder).c_str());

  TFile *inputfile = TFile::Open(inputfilename.c_str());
  TFile *outputfile = new TFile(outputfilename.c_str(), "RECREATE");
  std::cout << "inputfile = " << inputfile << std::endl;

  const int nPtBins = 9;
  const float pTBins[nPtBins][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 10}, {10, 20}};

  float pTMeans[2][nPtBins] = {0};
  {
    std::string input_csv_name = folder_csv_inputs + "/input_pt_mean_" + production_date + "_" + selection_label + ".csv";
    std::cout << "csv name = " << input_csv_name << std::endl;
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::getline(line_pi0, entry, ',');
      float pTMean = std::stof(entry);
      pTMeans[0][iPt] = pTMean;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::getline(line_eta, entry, ',');
      float pTMean = std::stof(entry);
      pTMeans[1][iPt] = pTMean;
    }
  }
    
  float bkg_ratio[2][2][nPtBins] = {0};
  {
    std::string input_csv_name = folder_csv_ratios + "/pt_ratios_" + production_date + "_" + selection_label + ratio_suffix + ".csv";
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::getline(line_pi0, entry, ',');
      float pTRatio = std::stof(entry);
      bkg_ratio[0][0][iPt] = pTRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::getline(line_eta, entry, ',');
      float pTRatio = std::stof(entry);
      bkg_ratio[0][1][iPt] = pTRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::getline(line_pi0_err, entry, ',');
      float pTRatioErr = std::stof(entry);
      bkg_ratio[1][0][iPt] = pTRatioErr;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::getline(line_eta_err, entry, ',');
      float pTRatioErr = std::stof(entry);
      bkg_ratio[1][1][iPt] = pTRatioErr;
    }
  }
  
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
  outputfile->cd();

  TGraphErrors *graph_AN_raw[nBeams][nParticles][nMethods] = {nullptr};
  TGraphErrors *graph_AN_raw_offset[nBeams][nParticles][nMethods] = {nullptr};
  TGraphErrors *graph_AN_bkg[nBeams][nParticles][nMethods] = {nullptr};
  TGraphErrors *graph_AN_bkg_offset[nBeams][nParticles][nMethods] = {nullptr};
  TGraphErrors *graph_AN_peak[nBeams][nParticles][nMethods] = {nullptr};
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int imethod = 0; imethod < nMethods; imethod++)
      {
        graph_AN_raw[iB][iP][imethod]  = new TGraphErrors();
        std::stringstream graph_name;
        graph_name << "graph_AN_raw_" << beams[iB] << "_" << particle[iP] << "_" << methods[imethod];
        graph_AN_raw[iB][iP][imethod]->SetName(graph_name.str().c_str());

        graph_AN_raw_offset[iB][iP][imethod]  = new TGraphErrors();
        graph_name.str("");
        graph_name << "graph_AN_raw_offset_" << beams[iB] << "_" << particle[iP] << "_" << methods[imethod];
        graph_AN_raw_offset[iB][iP][imethod]->SetName(graph_name.str().c_str());

        graph_AN_bkg[iB][iP][imethod]  = new TGraphErrors();
        graph_name.str("");
        graph_name << "graph_AN_bkg_" << beams[iB] << "_" << particle[iP] << "_" << methods[imethod];
        graph_AN_bkg[iB][iP][imethod]->SetName(graph_name.str().c_str());

        graph_AN_bkg_offset[iB][iP][imethod]  = new TGraphErrors();
        graph_name.str("");
        graph_name << "graph_AN_bkg_offset_" << beams[iB] << "_" << particle[iP] << "_" << methods[imethod];
        graph_AN_bkg_offset[iB][iP][imethod]->SetName(graph_name.str().c_str());

        graph_AN_peak[iB][iP][imethod]  = new TGraphErrors();
        graph_name.str("");
        graph_name << "graph_AN_peak_" << beams[iB] << "_" << particle[iP] << "_" << methods[imethod];
        graph_AN_peak[iB][iP][imethod]->SetName(graph_name.str().c_str());

        for (int iPt = 0; iPt < nPtBins; iPt++)
        {
          float mean_pT = pTMeans[iP][iPt];
          graph_AN_raw[iB][iP][imethod]->SetPoint(iPt, mean_pT, 0);
          graph_AN_raw_offset[iB][iP][imethod]->SetPoint(iPt, mean_pT, 0);
          graph_AN_bkg[iB][iP][imethod]->SetPoint(iPt, mean_pT, 0);
          graph_AN_bkg_offset[iB][iP][imethod]->SetPoint(iPt, mean_pT, 0);
          graph_AN_peak[iB][iP][imethod]->SetPoint(iPt, mean_pT, 0);
        }
      }
    }
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
            std::stringstream graph_AN_title;
            graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_pT_" << iPt << "_" << methods[imethod];
            graph_AN_title << "P_{T} #in [" << (int) pTBins[iPt][0] << ", " << (int) pTBins[iPt][1] << "] GeV/c (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.")  << ", " << beams[iB] << " beam); #phi; " << (iR == 0 ? "Raw Asymmetry" : "Bkg Asymmetry");

            TGraphErrors *graph_AN;
            inputfile->GetObject(graph_AN_name.str().c_str(), graph_AN);
            float epsilon = 0;
            float depsilon = 0;
            float offset = 0;
            float doffset = 0;
            PlotOneAsym(methods[imethod].c_str(), plots_folder, graph_AN_name.str(), graph_AN_title.str(), graph_AN, epsilon, depsilon, offset, doffset);
            graph_AN->Write();

            if (iR == 0)
            {
              graph_AN_raw[iB][iP][imethod]->SetPointY(iPt, epsilon);
              graph_AN_raw[iB][iP][imethod]->SetPointError(iPt, 0, depsilon);
              graph_AN_raw_offset[iB][iP][imethod]->SetPointY(iPt, offset);
              graph_AN_raw_offset[iB][iP][imethod]->SetPointError(iPt, 0, doffset);
            }
            else
            {
              graph_AN_bkg[iB][iP][imethod]->SetPointY(iPt, epsilon);
              graph_AN_bkg[iB][iP][imethod]->SetPointError(iPt, 0, depsilon);
              graph_AN_bkg_offset[iB][iP][imethod]->SetPointY(iPt, offset);
              graph_AN_bkg_offset[iB][iP][imethod]->SetPointError(iPt, 0, doffset);
            }
          }
        }
      }
    }
  }


  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iPt = 0; iPt < nPtBins; iPt++)
      {
        for (int imethod = 0; imethod < nMethods; imethod++)
        {
          float R = bkg_ratio[0][iP][iPt];
          R = (R == 0 ? 1 : R);
          float AN_raw = graph_AN_raw[iB][iP][imethod]->GetPointY(iPt);
          float AN_raw_err = graph_AN_raw[iB][iP][imethod]->GetErrorY(iPt);
          float AN_bkg = graph_AN_bkg[iB][iP][imethod]->GetPointY(iPt);
          float AN_bkg_err = graph_AN_bkg[iB][iP][imethod]->GetErrorY(iPt);
          float R_err = bkg_ratio[1][iP][iPt];
          float Delta_R = 0;//(AN_raw - AN_bkg) / std::pow(1 - R, 2) * R_err;
          float Delta_AN_raw = 1 / (1 - R) * AN_raw_err;
          float Delta_AN_bkg = -R / (1 - R) * AN_bkg_err;
          float AN_peak_val = R < 1 ? (AN_raw - R * AN_bkg) / (1 - R) : 0;
          float AN_peak_err = R < 1 ? std::sqrt(std::pow(Delta_R, 2) + std::pow(Delta_AN_raw, 2) + std::pow(Delta_AN_bkg, 2)) : 0;
          graph_AN_peak[iB][iP][imethod]->SetPointY(iPt, AN_peak_val);
          graph_AN_peak[iB][iP][imethod]->SetPointError(iPt, 0, AN_peak_err);
        }
      }
    }
  }

  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int imethod = 0; imethod < nMethods; imethod++)
      {
        std::stringstream graph_AN_title;
        graph_AN_title << beams[iB] << " asymmetry; p_{T} [GeV/c]; Raw asymmetry";
        graph_AN_raw[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_raw[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_raw[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_raw[iB][iP][imethod]->Write();
        //CreatePlotTGraphErrors(plots_folder, graph_AN_raw[iB][iP][imethod]);

        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; p_{T} [GeV/c]; Raw constant offset";
        graph_AN_raw_offset[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_raw_offset[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_raw_offset[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_raw_offset[iB][iP][imethod]->Write();

        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; p_{T} [GeV/c]; Background asymmetry";
        graph_AN_bkg[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_bkg[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_bkg[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_bkg[iB][iP][imethod]->Write();
        //CreatePlotTGraphErrors(plots_folder, graph_AN_bkg[iB][iP][imethod]);

        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; p_{T} [GeV/c]; Background constant offset";
        graph_AN_bkg_offset[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_bkg_offset[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_bkg_offset[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_bkg_offset[iB][iP][imethod]->Write();

        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; p_{T} [GeV/c]; Corrected asymmetry";
        graph_AN_peak[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_peak[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_peak[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_peak[iB][iP][imethod]->Write();
        //CreatePlotTGraphErrors(plots_folder, graph_AN_peak[iB][iP][imethod]);
      }
    }
  }

  outputfile->Close();
  delete outputfile;

  gSystem->Exit(0);
}

void PlotOneAsym(std::string type, std::string plots_folder, std::string canvas_name, std::string graph_title, TGraphErrors* gr, float& epsilon, float& depsilon, float& offset, float& doffset) {
  TCanvas *canvas = nullptr;
  if (printcanvas)
  {
    canvas = new TCanvas(canvas_name.c_str());
    canvas->cd();
  }
  double fit_lower, fit_upper;
  if (type == "sqrt") {
    fit_lower = -M_PI/2.0;
    fit_upper = M_PI/2.0;
  }
  else {
    fit_lower = -M_PI;
    fit_upper = M_PI;
  }
  
  TF1* fit = new TF1("fitasym0", "[0]", fit_lower, fit_upper);
  //fit->SetParLimits(1,-0.2,0.2);
  fit->SetParName(0, "C_{0}");
  //fit->SetParName(1, "C_{0}");
  
  gr->SetTitle(graph_title.c_str());
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gr->GetXaxis()->SetLimits(fit_lower, fit_upper);
  if (printcanvas)
  {
    gr->Draw("ap0");
  }
  gr->Fit("fitasym0", "RQ");
  offset = fit->GetParameter(0);
  doffset = fit->GetParError(0);
  //std::cout << "offset = " << offset << std::endl;
  //std::cout << "doffset = " << doffset << std::endl;
  delete fit;

  fit = new TF1("fitasym", "-[0]*sin(x)+[1]", fit_lower, fit_upper);
  //fit->FixParameter(1, offset);
  //fit->SetParError(1, doffset);
  fit->FixParameter(1, 0);
  fit->SetParError(1, 0);
  //fit->SetParError(0, depsilon);

  gr->Fit("fitasym", "RQ");
  epsilon = fit->GetParameter(0);
  depsilon = fit->GetParError(0);
  offset = fit->GetParameter(1);
  //doffset = fit->GetParError(1);
  //gr->GetXaxis()->SetRangeUser(-M_PI/2, M_PI/2);
  //gr->GetYaxis()->SetRangeUser(-0.05, 0.05);
  if (printcanvas)
  {
    canvas->Draw();
    canvas->SaveAs((plots_folder + "/" + canvas_name + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name + ".png").c_str());
    delete canvas;
  }


}
