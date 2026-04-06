#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

bool printcanvas = false;

void PlotOneAsym(std::string type, std::string plots_folder, std::string canvas_name, std::string graph_title, TGraphErrors* gr, float& epsilon, float& depsilon, float& offset, float& doffset);

void fit_sinusoid_asym_xf_ana509(const std::string &production_date = "01312026",
                                 const std::string &selection_label = "MBD_pt_l4_0mrad",
                                 const std::string &folder_asymmetry = "",
                                 const std::string &folder_csv_inputs = "",
                                 const std::string &folder_csv_ratios = "",
                                 bool gpr_ratio = false)
{
  std::string ratio_suffix = (gpr_ratio ? "_gpr" : "");
  std::string inputfilename = folder_asymmetry + "/analysis_asymmetry_xf_complete_ana509_" + production_date + "_" + selection_label + "_averaged.root";
  std::string outputfilename = folder_asymmetry + "/analysis_asymmetry_xf_complete_ana509_" + production_date + "_" + selection_label + "_averaged_fitted" + ratio_suffix + ".root";

  std::string plots_folder = "plots_asymmetry/ana509_p022/analysis_" + production_date + "/xf_asymmetries/" + selection_label + "/phi_asymmetries/runs";
  gSystem->Exec(("mkdir -p " + plots_folder).c_str());

  TFile *inputfile = TFile::Open(inputfilename.c_str());
  TFile *outputfile = new TFile(outputfilename.c_str(), "RECREATE");
  std::cout << "inputfile = " << inputfile << std::endl;

  const int nXfBins = 8;
  const float xfBins[nXfBins][2] = {{-0.20, -0.048}, {-0.048, -0.035}, {-0.035, -0.022}, {-0.022, 0.0}, {0, 0.022}, {0.022, 0.035}, {0.035, 0.048}, {0.048, 0.20}};

  float xfMeans[2][nXfBins] = {0};
  {
    std::string input_csv_name = folder_csv_inputs + "/input_xf_mean_" + production_date + "_" + selection_label + ".csv";
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
    for (int ixf = 0; ixf < nXfBins; ixf++) {
      std::getline(line_pi0, entry, ',');
      float xfMean = std::stof(entry);
      xfMeans[0][ixf] = xfMean;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int ixf = 0; ixf < nXfBins; ixf++) {
      std::getline(line_eta, entry, ',');
      float xfMean = std::stof(entry);
      xfMeans[1][ixf] = xfMean;
    }
  }
    
  float bkg_ratio[2][2][nXfBins] = {0};
  {
    std::string input_csv_name = folder_csv_ratios + "/xf_ratios_" + production_date + "_" + selection_label + ratio_suffix + ".csv";
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
    for (int ixf = 0; ixf < nXfBins; ixf++) {
      std::getline(line_pi0, entry, ',');
      float xfRatio = std::stof(entry);
      bkg_ratio[0][0][ixf] = xfRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int ixf = 0; ixf < nXfBins; ixf++) {
      std::getline(line_eta, entry, ',');
      float xfRatio = std::stof(entry);
      bkg_ratio[0][1][ixf] = xfRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int ixf = 0; ixf < nXfBins; ixf++) {
      std::getline(line_pi0_err, entry, ',');
      float xfRatioErr = std::stof(entry);
      bkg_ratio[1][0][ixf] = xfRatioErr;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int ixf = 0; ixf < nXfBins; ixf++) {
      std::getline(line_eta_err, entry, ',');
      float xfRatioErr = std::stof(entry);
      bkg_ratio[1][1][ixf] = xfRatioErr;
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
        
        for (int ixf = 0; ixf < nXfBins; ixf++)
        {
          float mean_xf = xfMeans[iP][ixf];
          graph_AN_raw[iB][iP][imethod]->SetPoint(ixf, mean_xf, 0);
          graph_AN_raw_offset[iB][iP][imethod]->SetPoint(ixf, mean_xf, 0);
          graph_AN_bkg[iB][iP][imethod]->SetPoint(ixf, mean_xf, 0);
          graph_AN_bkg_offset[iB][iP][imethod]->SetPoint(ixf, mean_xf, 0);
          graph_AN_peak[iB][iP][imethod]->SetPoint(ixf, mean_xf, 0);
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
        for (int ixf = 0; ixf < nXfBins; ixf++)
        {
          for (int imethod = 0; imethod < nMethods; imethod++)
          {
            std::stringstream graph_AN_name;
            std::stringstream graph_AN_title;
            graph_AN_name << "graph_asym_" << beams[iB] << "_" << particle[iP] << "_" << regions[iR] << "_xf_" << ixf << "_" << methods[imethod];
            graph_AN_title << "x_F in [" << std::fixed << std::setprecision(1) << xfBins[ixf][0] << ", " << xfBins[ixf][1] << "] (" << (imethod == 0 ? "Geometric" : "Rel. Lumi.")  << ", " << beams[iB] << " beam); #phi; " << (iR == 0 ? "Raw Asymmetry" : "Bkg Asymmetry");
            
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
              graph_AN_raw[iB][iP][imethod]->SetPointY(ixf, epsilon);
              graph_AN_raw[iB][iP][imethod]->SetPointError(ixf, 0, depsilon);
              graph_AN_raw_offset[iB][iP][imethod]->SetPointY(ixf, offset);
              graph_AN_raw_offset[iB][iP][imethod]->SetPointError(ixf, 0, doffset);
            }
            else
            {
              graph_AN_bkg[iB][iP][imethod]->SetPointY(ixf, epsilon);
              graph_AN_bkg[iB][iP][imethod]->SetPointError(ixf, 0, depsilon);
              graph_AN_bkg_offset[iB][iP][imethod]->SetPointY(ixf, offset);
              graph_AN_bkg_offset[iB][iP][imethod]->SetPointError(ixf, 0, doffset);
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
      for (int ixf = 0; ixf < nXfBins; ixf++)
      {
        for (int imethod = 0; imethod < nMethods; imethod++)
        {
          float R = bkg_ratio[0][iP][ixf];
          R = (R == 0 ? 1 : R);
          float AN_raw = graph_AN_raw[iB][iP][imethod]->GetPointY(ixf);
          float AN_raw_err = graph_AN_raw[iB][iP][imethod]->GetErrorY(ixf);
          float AN_bkg = graph_AN_bkg[iB][iP][imethod]->GetPointY(ixf);
          float AN_bkg_err = graph_AN_bkg[iB][iP][imethod]->GetErrorY(ixf);
          float R_err = bkg_ratio[1][iP][ixf];
          float Delta_R = 0;//(AN_raw - AN_bkg) / std::pow(1 - R, 2) * R_err;
          float Delta_AN_raw = 1 / (1 - R) * AN_raw_err;
          float Delta_AN_bkg = -R / (1 - R) * AN_bkg_err;
          float AN_peak_val = R < 1 ? (AN_raw - R * AN_bkg) / (1 - R) : 0;
          float AN_peak_err = R < 1 ? std::sqrt(std::pow(Delta_R, 2) + std::pow(Delta_AN_raw, 2) + std::pow(Delta_AN_bkg, 2)) : 0;
          graph_AN_peak[iB][iP][imethod]->SetPointY(ixf, AN_peak_val);
          graph_AN_peak[iB][iP][imethod]->SetPointError(ixf, 0, AN_peak_err);
          if (iB == 1 && iP == 1 && ixf == 0 && imethod == 0) {
            std::cout << "(" << AN_raw << ", " << AN_raw_err << ", " << AN_bkg << ", " << AN_bkg_err << ", " << R << ", " << Delta_R << ", " << Delta_AN_raw << ", " << Delta_AN_bkg << ", " << AN_peak_val  << ", " << AN_peak_err << ")" << std::endl;
          }
        }
      }
    }
  }

  outputfile->cd();
  std::cout << "Write and draw" << std::endl;
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int imethod = 0; imethod < nMethods; imethod++)
      {
        std::stringstream graph_AN_title;
        graph_AN_title << beams[iB] << " asymmetry; x_F; Raw asymmetry";
        graph_AN_raw[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_raw[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_raw[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_raw[iB][iP][imethod]->Write();
        //CreatePlotTGraphErrors(plots_folder, graph_AN_raw[iB][iP][imethod]);
        
        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; x_F; Raw constant offset";
        graph_AN_raw_offset[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_raw_offset[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_raw_offset[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_raw_offset[iB][iP][imethod]->Write();
        
        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; x_F; Background asymmetry";
        graph_AN_bkg[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_bkg[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_bkg[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_bkg[iB][iP][imethod]->Write();
        //CreatePlotTGraphErrors(plots_folder, graph_AN_bkg[iB][iP][imethod]);
        
        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; x_F; Background constant offset";
        graph_AN_bkg_offset[iB][iP][imethod]->SetTitle(graph_AN_title.str().c_str());
        if (iP == 0) graph_AN_bkg_offset[iB][iP][imethod]->GetXaxis()->SetLimits(1,10);
        else graph_AN_bkg_offset[iB][iP][imethod]->GetXaxis()->SetLimits(2,20);
        graph_AN_bkg_offset[iB][iP][imethod]->Write();
        
        graph_AN_title.str("");
        graph_AN_title << beams[iB] << " asymmetry; x_F; Corrected asymmetry";
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

  fit = new TF1("fitasym", "-[0]*sin(x)", fit_lower, fit_upper);
  //fit->FixParameter(1, offset);
  //fit->SetParError(1, doffset);
  //fit->FixParameter(1, 0);
  //fit->SetParError(1, 0);
  //fit->SetParError(0, depsilon);

  gr->Fit("fitasym", "RQ");
  epsilon = fit->GetParameter(0);
  depsilon = fit->GetParError(0);
  //offset = fit->GetParameter(1);
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
