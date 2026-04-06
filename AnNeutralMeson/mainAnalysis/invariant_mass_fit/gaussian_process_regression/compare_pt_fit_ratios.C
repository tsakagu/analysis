#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

double IntegrateGraphInRange(const TGraph* g, double xmin, double xmax)
{
  if (!g) return 0.0;
  int n = g->GetN();

  const double* x = g->GetX();
  const double* y = g->GetY();

  double area = 0.0;

  for (int i = 0; i < n - 1; ++i) {
    double x1 = x[i],   y1 = y[i];
    double x2 = x[i+1], y2 = y[i+1];
    
    if (x2 <= xmin || x1 >= xmax) continue;
    
    // clip segment to [xmin, xmax]
    double xa = std::max(x1, xmin);
    double xb = std::min(x2, xmax);
    if (xb <= xa) continue;
    
    // linear interpolation at clipped boundaries
    double ya = y1 + (y2 - y1) * (xa - x1) / (x2 - x1);
    double yb = y1 + (y2 - y1) * (xb - x1) / (x2 - x1);
    
    area += 0.5 * (ya + yb) * (xb - xa);
  }
  return area;
}

void compare_pt_fit_ratios(const std::string& selection_label="MBD_pt_l3_0mrad",
                           const std::string& inputfolder_mass="",
                           const std::string& inputfolder_parametric_base="",
                           const std::string& inputfolder_gpr_base="",
                           const std::string& outputfolder_base="")
{
  const std::string production_date = "01312026";
  const int nParticles = 2;
  const int nRegions = 3; // left center right
  const int nPtBins = 9;
  const float pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 20};

  // Read mean pt values
  std::string input_csv_name = inputfolder_mass + "/csv_inputs/input_pt_mean_" + production_date + "_" + selection_label + ".csv";
  float ptMeans[2][nPtBins];
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
    float ptMean = std::stof(entry);
    ptMeans[0][iPt] = ptMean;
  }
  std::getline(input_csv, line);
  std::stringstream line_eta(line);
  for (int iPt = 0; iPt < nPtBins; iPt++) {
    std::getline(line_eta, entry, ',');
    float ptMean = std::stof(entry);
    ptMeans[1][iPt] = ptMean;
  }
  
  std::string inputfolder_parametric = inputfolder_parametric_base + "/" + selection_label + "/";
  std::string inputfolder_gpr = inputfolder_gpr_base + "/" + selection_label + "/";
  std::string outputfolder = outputfolder_base + "/" + selection_label + "/";

  gSystem->Exec(("mkdir -p " + outputfolder).c_str());

  TGraphErrors *pt_parametric_ratios_pi0 = new TGraphErrors();
  pt_parametric_ratios_pi0->SetName("pt_background_parametric_ratios_pi0");
  pt_parametric_ratios_pi0->SetTitle("#pi^{0} background ratio");
  pt_parametric_ratios_pi0->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_parametric_ratios_pi0->GetYaxis()->SetTitle("background fraction [%]");
  
  TGraphErrors *pt_parametric_ratios_eta = new TGraphErrors();
  pt_parametric_ratios_eta->SetTitle("#eta background ratio");
  pt_parametric_ratios_eta->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_parametric_ratios_eta->GetYaxis()->SetTitle("background ratio [%]");
  pt_parametric_ratios_eta->SetName("pt_background_parametric_ratios_eta");

  TGraphErrors *pt_gpr_ratios_pi0 = new TGraphErrors();
  pt_gpr_ratios_pi0->SetName("pt_background_gpr_ratios_pi0");
  pt_gpr_ratios_pi0->SetTitle("#pi^{0} background ratio");
  pt_gpr_ratios_pi0->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_gpr_ratios_pi0->GetYaxis()->SetTitle("background fraction [%]");
  
  TGraphErrors *pt_gpr_ratios_eta = new TGraphErrors();
  pt_gpr_ratios_eta->SetTitle("#eta background ratio");
  pt_gpr_ratios_eta->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_gpr_ratios_eta->GetYaxis()->SetTitle("background ratio [%]");
  pt_gpr_ratios_eta->SetName("pt_background_gpr_ratios_eta");
  
  // 3 sigma bands
  float band_bounds[nParticles][nRegions][2] =
    {{{0.080,0.199}, {0.030, 0.070}, {0.209, 0.249}},
     {{0.399,0.739}, {0.257, 0.371}, {0.767, 0.880}}};
  
  for (int iPt = 0; iPt < nPtBins; iPt++) { 
    std::string inputfilename_parametric = inputfolder_parametric + "/parametric_fit_mass_pt_" + std::to_string(iPt) + ".root";
    TFile *infile_parametric = TFile::Open(inputfilename_parametric.c_str());

    std::string inputfilename_gpr = inputfolder_gpr + "/gpr_fit_mass_pt_" + std::to_string(iPt) + ".root";
    TFile *infile_gpr = TFile::Open(inputfilename_gpr.c_str());

    if (!infile_parametric || infile_parametric->IsZombie() ||
        !infile_gpr || infile_gpr->IsZombie()) {
      pt_parametric_ratios_pi0->SetPoint(iPt, ptMeans[0][iPt], 100);
      pt_parametric_ratios_pi0->SetPointError(iPt, 0, 0);
      pt_parametric_ratios_eta->SetPoint(iPt, ptMeans[1][iPt], 100);
      pt_parametric_ratios_eta->SetPointError(iPt, 0, 0);
      pt_gpr_ratios_pi0->SetPoint(iPt, ptMeans[0][iPt], 100);
      pt_gpr_ratios_pi0->SetPointError(iPt, 0, 0);
      pt_gpr_ratios_eta->SetPoint(iPt, ptMeans[1][iPt], 100);
      pt_gpr_ratios_eta->SetPointError(iPt, 0, 0);
      continue;
    }

    // Extract info from parametric fit
    TH1F *h_pair_mass = (TH1F*)infile_parametric->Get("h_pair_mass");

    TF1 *fit_global_bkg_pi0 = (TF1*)infile_parametric->Get("fit_global_bkg_pi0");
    TF1 *fit_global_bkg_eta = (TF1*)infile_parametric->Get("fit_global_bkg_eta");
    TF1 *fit_global_pi0 = (TF1*)infile_parametric->Get("fit_global_pi0");
    TF1 *fit_global_eta = (TF1*)infile_parametric->Get("fit_global_eta");

    // Extract info from GPR fit
    // gpr_fit_mass_pt_${pt_index}.root
    TGraph *fit_gpr_bkg = (TGraph*)infile_gpr->Get(("graph_gpr_mass_fit_pt_" + std::to_string(iPt) + "_bkg").c_str());
    fit_gpr_bkg->SetLineWidth(3);
    fit_gpr_bkg->SetLineColor(kGreen+3);
    fit_gpr_bkg->SetLineStyle(kDashed);
    TGraph *fit_gpr_signal = (TGraph*)infile_gpr->Get(("graph_gpr_mass_fit_pt_" + std::to_string(iPt) + "_signal").c_str());
    fit_gpr_signal->SetLineWidth(3);
    fit_gpr_signal->SetLineColor(kCyan);
    fit_gpr_signal->SetLineStyle(kDashed);
    TGraph *fit_gpr_total = new TGraph();
    fit_gpr_total->SetName("fit_gpr_total");
    int nPoints = fit_gpr_bkg->GetN();
    for (int i=0; i < nPoints; i++) {
      double x = fit_gpr_bkg->GetPointX(i);
      double y1 = fit_gpr_bkg->GetPointY(i);
      double y2 = fit_gpr_signal->GetPointY(i);
      double y = y1 + y2;
      fit_gpr_total->SetPoint(i, x, y);
    }
    fit_gpr_total->SetLineWidth(3);
    fit_gpr_total->SetLineColor(kRed);
    fit_gpr_total->SetLineStyle(kSolid);

    // Compute background fractions for both methods
    float parametric_pi0_ratio = fit_global_bkg_pi0->Integral(band_bounds[0][0][0], band_bounds[0][0][1]) /
      (fit_global_pi0->Integral(band_bounds[0][0][0], band_bounds[0][0][1]) +
       fit_global_bkg_pi0->Integral(band_bounds[0][0][0], band_bounds[0][0][1]));
    float parametric_eta_ratio = fit_global_bkg_eta->Integral(band_bounds[1][0][0], band_bounds[1][0][1]) /
      (fit_global_eta->Integral(band_bounds[1][0][0], band_bounds[1][0][1]) +
       fit_global_bkg_eta->Integral(band_bounds[1][0][0], band_bounds[1][0][1]));
    float gpr_pi0_ratio = IntegrateGraphInRange(fit_gpr_bkg, band_bounds[0][0][0], band_bounds[0][0][1]) /
      (IntegrateGraphInRange(fit_gpr_signal, band_bounds[0][0][0], band_bounds[0][0][1]) +
       IntegrateGraphInRange(fit_gpr_bkg, band_bounds[0][0][0], band_bounds[0][0][1]));
    float gpr_eta_ratio = IntegrateGraphInRange(fit_gpr_bkg, band_bounds[1][0][0], band_bounds[1][0][1]) /
      (IntegrateGraphInRange(fit_gpr_signal, band_bounds[1][0][0], band_bounds[1][0][1]) +
       IntegrateGraphInRange(fit_gpr_bkg, band_bounds[1][0][0], band_bounds[1][0][1]));
    
    
    pt_parametric_ratios_pi0->SetPoint(iPt, ptMeans[0][iPt], parametric_pi0_ratio * 100);
    pt_parametric_ratios_pi0->SetPointError(iPt, 0, 0);
    pt_parametric_ratios_eta->SetPoint(iPt, ptMeans[1][iPt], parametric_eta_ratio * 100);
    pt_parametric_ratios_eta->SetPointError(iPt, 0, 0);
    pt_gpr_ratios_pi0->SetPoint(iPt, ptMeans[0][iPt], gpr_pi0_ratio * 100);
    pt_gpr_ratios_pi0->SetPointError(iPt, 0, 0);
    pt_gpr_ratios_eta->SetPoint(iPt, ptMeans[1][iPt], gpr_eta_ratio * 100);
    pt_gpr_ratios_eta->SetPointError(iPt, 0, 0);

    SetsPhenixStyle();
    std::stringstream canvas_name;
    canvas_name << "canvas_comparison_fit_mass_pt_" << iPt << "_" << selection_label;
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "c", 1600, 900);
    canvas->cd();

    h_pair_mass->SetLineColor(kBlack);
    h_pair_mass->SetStats(0);
    h_pair_mass->Draw();

    gPad->Update();
    double max_y = gPad->GetUymax();

    TBox *pi0PeakBand = new TBox(band_bounds[0][0][0], 0, band_bounds[0][0][1], max_y);
    pi0PeakBand->SetFillColorAlpha(kRed, 0.4);
    TBox *pi0SideBandLeft = new TBox(band_bounds[0][1][0], 0, band_bounds[0][1][1], max_y);
    pi0SideBandLeft->SetFillColorAlpha(kRed, 0.2);
    TBox *pi0SideBandRight = new TBox(band_bounds[0][2][0], 0, band_bounds[0][2][1], max_y);
    pi0SideBandRight->SetFillColorAlpha(kRed, 0.2);

    TBox *etaPeakBand = new TBox(band_bounds[1][0][0], 0, band_bounds[1][0][1], max_y);
    etaPeakBand->SetFillColorAlpha(kBlue, 0.4);
    TBox *etaSideBandLeft = new TBox(band_bounds[1][1][0], 0, band_bounds[1][1][1], max_y);
    etaSideBandLeft->SetFillColorAlpha(kBlue, 0.2);
    TBox *etaSideBandRight = new TBox(band_bounds[1][2][0], 0, band_bounds[1][2][1], max_y);
    etaSideBandRight->SetFillColorAlpha(kBlue, 0.2);

    fit_gpr_total->Draw("SAME");

    // Draw peak regions
    pi0PeakBand->Draw("same");
    pi0SideBandLeft->Draw("same");
    pi0SideBandRight->Draw("same");
    etaPeakBand->Draw("same");
    etaSideBandLeft->Draw("same");
    etaSideBandRight->Draw("same");
    
    // Draw the parametric signal and background
    fit_global_bkg_pi0->Draw("SAME");
    fit_global_bkg_eta->Draw("SAME");
    fit_global_pi0->Draw("SAME");
    fit_global_eta->Draw("SAME");

    // Draw the GPR signal and background
    fit_gpr_bkg->Draw("SAME");
    fit_gpr_signal->Draw("SAME");

    TPad *p = new TPad("p","p",0.,0.,1.,1.); p->SetFillStyle(0); p->Draw(); p->cd();
    TBox *whiteBox = new TBox(0.67, 0.72, 0.88, 0.9);
    whiteBox->Draw();
    canvas->cd();
    whiteBox->SetFillColorAlpha(kWhite, 1);
    std::stringstream stream;
    TLatex latex;
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.68, 0.85, "#font[72]{sPHENIX} Internal");
    latex.DrawLatex(0.68, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");
    stream.str("");
    stream << std::fixed << std::setprecision(0) << pTBins[iPt] << " < p_{T} [GeV] < " << pTBins[iPt+1];
    latex.DrawLatex(0.47,0.87, stream.str().c_str());
    
    canvas->Draw();
    canvas->SaveAs((outputfolder + "/" + canvas_name.str() + ".png").c_str());
    canvas->SaveAs((outputfolder + "/" + canvas_name.str() + ".pdf").c_str());
    delete canvas;
  }

  TFile *outfile_ratios = new TFile((outputfolder + "/pt_ratios.root").c_str(), "RECREATE");
  outfile_ratios->cd();
  pt_parametric_ratios_pi0->Write();
  pt_parametric_ratios_eta->Write();
  pt_gpr_ratios_pi0->Write();
  pt_gpr_ratios_eta->Write();
  outfile_ratios->Close();
  delete outfile_ratios;

  gSystem->Exit(0);
}
