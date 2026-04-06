#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

R__LOAD_LIBRARY(libCustomHistogramFit.so)
#include <CustomHistogramFit/HistogramUtils.h>
#include <CustomHistogramFit/FitFunctions.h>
#include <CustomHistogramFit/ConstrainedFitter.h>

void fit_argusmodified_2(const std::string &production_date = "01312026",
                         const std::string &selection_label = "MBD_INCOMPLETE",
                         const std::string &inputfolder_mass = "",
                         const std::string &inputfolder_analysis = "",
                         const std::string &folder_inputs = "",
                         const std::string &outputfolder_ratios = "",
                         const std::string &plots_folder_template = "")
{
  const int nPtBins = 9;
  const std::string pTBins[nPtBins + 1] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
  const std::string pTBounds[nPtBins][2] = {{"1","2"}, {"2","3"}, {"3", "4"}, {"4","5"}, {"5", "6"}, {"6","7"},{"7","8"},{"8","10"},{"10","20"}};

  // Read mean pT values
  std::string input_csv_name = folder_inputs + "/input_pt_mean_" + production_date + "_" + selection_label + ".csv";
  float pTMeans[2][nPtBins] = {0};
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

  // Histogram to store the background ratios

  TGraphErrors *pt_ratios_pi0 = new TGraphErrors();
  pt_ratios_pi0->SetName("pt_background_ratios_pi0");
  pt_ratios_pi0->SetTitle("#pi^{0} background ratio");
  pt_ratios_pi0->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_ratios_pi0->GetYaxis()->SetTitle("background fraction [%]");
  TGraphErrors *pt_ratios_eta = new TGraphErrors();
  pt_ratios_eta->SetTitle("#eta background ratio");
  pt_ratios_eta->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_ratios_eta->GetYaxis()->SetTitle("background ratio [%]");
  pt_ratios_eta->SetName("pt_background_ratios_eta");

  std::string plots_folder = plots_folder_template + "/plots_mass_argusmodified_2/analysis_ana509/ana509_" + production_date + "_sigma_3/" + selection_label + "/";
  
  gSystem->Exec(("mkdir -p " + plots_folder).c_str());
  double statSum = 0;

  for (int iPt = 0; iPt < nPtBins + 1; iPt++)
  {
    std::string filename = inputfolder_analysis + "/analysis_complete_ana509_" + production_date + "_" + selection_label + ".root";

    std::string hname = "h_pair_mass_pt_" + pTBins[iPt];
    if (iPt == nPtBins) hname = "h_pair_mass";
    gSystem->Exec(("mkdir -p " + plots_folder).c_str());
    SetsPhenixStyle();

    // Open file
    TFile *f = TFile::Open(filename.c_str());
    
    // Extract the unfitted histogram
    std::cout << "f = " << f << std::endl;
    TH1F *h_pair_mass = (TH1F*) f->Get(hname.c_str());
    std::cout << "hname = " << hname.c_str() << std::endl;
    std::cout << "h_pair_mass = " << h_pair_mass << std::endl;
    h_pair_mass->SetTitle("");
    h_pair_mass->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");

    statSum += h_pair_mass->GetEntries();
    std::ofstream statList;
    statList.open((plots_folder + "/statList.txt"), ios::out | ios::app);
    statList << "bin : p_T in " << pTBins[iPt] << " GeV/c: nentries = " << h_pair_mass->GetEntries() << std::endl;
    statList.close();
    
    // Rescale the histogram to make the fit easier
    double scaling_factor = h_pair_mass->Integral();
    h_pair_mass->Scale(1. / scaling_factor);


    // Determine when the histogram starts
    double thresh = HistogramUtils::GetMaxBinContentInRange(h_pair_mass, 0.005, 0.9) / 1000;
    double thresh_large = HistogramUtils::GetMaxBinContentInRange(h_pair_mass, 0.005, 0.9) / 10;
    std::cout << "thresholds (scaled): " << thresh << ", " << thresh_large << std::endl;
    const double offset = h_pair_mass->GetBinCenter(HistogramUtils::FindFirstOccupiedBin(h_pair_mass, thresh));
    const double offset_large = std::max(h_pair_mass->GetBinCenter(HistogramUtils::FindFirstOccupiedBin(h_pair_mass, thresh_large)), 0.10);

    if (true) {
      std::cout << "offsets = (" << offset << ", " << offset_large << ")\n";
    }
    
    // Define the limit between two fitting regions. There is some overlap
    double borderFitMin = 0.27;
    double borderFitMax = 0.35;
    double fitMax1 = 0.27;
    double fitMin2 = 0.35;

    // Prior background fit
    ConstrainedFitter *fitMethodBkg1 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::ARGUSMODIFIED});
    fitMethodBkg1->SetFitRange({{offset, offset_large}, {0.18, fitMax1}});
    if (iPt == 0) fitMethodBkg1->SetFitRange({{0.05, 0.11}, {0.20, 0.32}});
    fitMethodBkg1->SetInitialParams({0.002, offset, 0.15, 0.1, 0.1});
    fitMethodBkg1->SetParamBounds({{0.0001, 1.0}, {offset-0.001, offset+0.001}, {0.1, 0.3}, {0.001, 10}, {0.001, 10}});
    fitMethodBkg1->PerformFit();
    fitMethodBkg1->PrintResults();
    TF1 *fit_bkg_1 = fitMethodBkg1->GetFitFunction();

    ConstrainedFitter *fitMethodBkg2 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::POLY2});
    fitMethodBkg2->SetFitRange({{fitMin2, 0.45}, {0.75, 1.00}});
    fitMethodBkg2->SetInitialParams({0.5, 0.5, 0.5});
    fitMethodBkg2->SetParamBounds({{0.0, 10.0}, {0.0, 10.0}, {0.0, 10.0}});
    fitMethodBkg2->PerformFit();
    fitMethodBkg2->PrintResults();
    TF1 *fit_bkg_2 = fitMethodBkg2->GetFitFunction();

    // pi0 fit, fixed background
    ConstrainedFitter *fitMethodPi0 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::ARGUSMODIFIED, FitFunctions::FitModel::GAUSS});
    fitMethodPi0->SetFitRange({{offset, fitMax1}});
    fitMethodPi0->SetInitialParams({
        fit_bkg_1->GetParameter(0),
        offset,
        fit_bkg_1->GetParameter(2),
        fit_bkg_1->GetParameter(3),
        fit_bkg_1->GetParameter(4),
        0.01, 0.140, 0.02}); // Gaussian initial parameters
    fitMethodPi0->SetParamBounds({
        {fit_bkg_1->GetParameter(0) - 0.0001, fit_bkg_1->GetParameter(0) + 0.0001},
        {offset - 0.0001, offset + 0.0001},
        {fit_bkg_1->GetParameter(2) - 0.0001, fit_bkg_1->GetParameter(2) + 0.0001},
        {fit_bkg_1->GetParameter(3) - 0.0001, fit_bkg_1->GetParameter(2) + 0.0001},
        {fit_bkg_1->GetParameter(4) - 0.0001, fit_bkg_1->GetParameter(3) + 0.0001},
        {0.001, 1}, {0.135, 0.160}, {0.015, 0.025}}); // Gaussian parameter limits
    fitMethodPi0->PerformFit();
    fitMethodPi0->PrintResults();
    TF1 *fit_pi0 = fitMethodPi0->GetFitFunction();

    // eta fit, fixed background
    ConstrainedFitter *fitMethodEta = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::POLY2, FitFunctions::FitModel::GAUSS});
    fitMethodEta->SetFitRange({{fitMin2, 0.75}});
    fitMethodEta->SetInitialParams({
        fit_bkg_2->GetParameter(0),
        fit_bkg_2->GetParameter(1),
        fit_bkg_2->GetParameter(2),
        0.01, 0.580, 0.05}); // Gaussian initial parameters
    fitMethodEta->SetParamBounds({
        {fit_bkg_2->GetParameter(0) - 0.0001, fit_bkg_2->GetParameter(0) + 0.0001},
        {fit_bkg_2->GetParameter(1) - 0.0001, fit_bkg_2->GetParameter(1) + 0.0001},
        {fit_bkg_2->GetParameter(2) - 0.0001, fit_bkg_2->GetParameter(2) + 0.0001},
        {0.0001, 1}, {0.560, 0.590}, {0.045, 0.055}}); // Gaussian parameter limits
    fitMethodEta->PerformFit();
    fitMethodEta->PrintResults();
    TF1 *fit_eta = fitMethodEta->GetFitFunction();
    

    // bkg, pi0 and eta fit
    ConstrainedFitter *fitMethodGlobal1 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::ARGUSMODIFIED, FitFunctions::FitModel::GAUSS});
    fitMethodGlobal1->SetFitRange({{offset, fitMax1}});
    fitMethodGlobal1->SetInitialParams({
        fit_pi0->GetParameter(0),
        offset,
        fit_pi0->GetParameter(2),
        fit_pi0->GetParameter(3),
        fit_pi0->GetParameter(4),
        fit_pi0->GetParameter(5),
        fit_pi0->GetParameter(6),
        fit_pi0->GetParameter(7)
      }); 
    fitMethodGlobal1->SetParamBounds({
        {0.0001, 1.0}, {offset-0.001, offset+0.001}, {0.1, 0.3}, {0.001, 10}, {0.001, 10},
        {0.001, 1}, {0.135, 0.160}, {0.015, 0.025}, // Pi0 Gaussian parameter limits
      });
    
    fitMethodGlobal1->PerformFit();
    fitMethodGlobal1->PrintResults();
    TF1 *fit_global_1 = fitMethodGlobal1->GetFitFunction();

    ConstrainedFitter *fitMethodGlobal2 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::POLY2, FitFunctions::FitModel::GAUSS});
    fitMethodGlobal2->SetFitRange({{fitMin2, 1.0}});
    fitMethodGlobal2->SetInitialParams({
        fit_eta->GetParameter(0),
        fit_eta->GetParameter(1),
        fit_eta->GetParameter(2),
        fit_eta->GetParameter(3),
        fit_eta->GetParameter(4),
        fit_eta->GetParameter(5)
      }); 
    fitMethodGlobal2->SetParamBounds({
        {0, 10}, {0, 10}, {0, 10},
        {0.0001, 1}, {0.560, 0.590}, {0.045, 0.055} // eta Gaussian parameter limits
      });
    fitMethodGlobal2->PerformFit();
    fitMethodGlobal2->PrintResults();
    TF1 *fit_global_2 = fitMethodGlobal2->GetFitFunction();

    // Rescale
    fit_global_1->SetParameter(0, fit_global_1->GetParameter(0) * scaling_factor);
    fit_global_1->SetParameter(5, fit_global_1->GetParameter(5) * scaling_factor);

    fit_global_2->SetParameter(0, fit_global_2->GetParameter(0) * scaling_factor);
    fit_global_2->SetParameter(1, fit_global_2->GetParameter(1) * scaling_factor);
    fit_global_2->SetParameter(2, fit_global_2->GetParameter(2) * scaling_factor);
    fit_global_2->SetParameter(3, fit_global_2->GetParameter(3) * scaling_factor);

    // Extract the bkg, pi0 and eta contribution from the last fit
    ConstrainedFitter *fitMethodGlobalBkg1 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::ARGUSMODIFIED});
    fitMethodGlobalBkg1->SetInitialParams({
        fit_global_1->GetParameter(0),
        fit_global_1->GetParameter(1),
        fit_global_1->GetParameter(2),
        fit_global_1->GetParameter(3),
        fit_global_1->GetParameter(4)
      });
    TF1 *fit_global_bkg_1 = fitMethodGlobalBkg1->GetFitFunction();

    ConstrainedFitter *fitMethodGlobalBkg2 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::POLY2});
    fitMethodGlobalBkg2->SetInitialParams({
        fit_global_2->GetParameter(0),
        fit_global_2->GetParameter(1),
        fit_global_2->GetParameter(2)
      });
    TF1 *fit_global_bkg_2 = fitMethodGlobalBkg2->GetFitFunction();
    fit_global_bkg_2->SetName("fit_global_bkg_2");
    

    //auto lambda_scaled_global_bkg_2 = [scaling_factor, fit_global_bkg_2] (double *x, double *p) { return scaling_factor * fit_global_bkg_2->Eval(x[0]); };
    //TF1 *fit_scaled_global_bkg_2 = new TF1("fit_scaled_global_bkg_2", lambda_scaled_global_bkg_2);
    /*TF1 *fit_scaled_global_bkg_2 = new TF1("fit_scaled_global_bkg_2", "[0] * fit_global_bkg_2");
      fit_scaled_global_bkg_2->SetParameter(0, scaling_factor);*/

    ConstrainedFitter *fitMethodGlobalPi0 = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::GAUSS});
    fitMethodGlobalPi0->SetInitialParams({
        fit_global_1->GetParameter(5),
        fit_global_1->GetParameter(6),
        fit_global_1->GetParameter(7)
      });
    TF1 *fit_global_pi0 = fitMethodGlobalPi0->GetFitFunction();

    ConstrainedFitter *fitMethodGlobalEta = new ConstrainedFitter(h_pair_mass, {FitFunctions::FitModel::GAUSS});
    fitMethodGlobalEta->SetInitialParams({
        fit_global_2->GetParameter(3),
        fit_global_2->GetParameter(4),
        fit_global_2->GetParameter(5)
      });
    TF1 *fit_global_eta = fitMethodGlobalEta->GetFitFunction();


    // Define transparent bands for the peak regions and the side bands
    const int nParticles = 2;
    const int nRegions = 3;

    // 3 sigma
    float band_bounds[nParticles][nRegions][2] =
    {{{0.080,0.199}, {0.030, 0.070}, {0.209, 0.249}},
    {{0.399,0.739}, {0.257, 0.371}, {0.767, 0.880}}};

    // 2 sigma
    /*float band_bounds[nParticles][nRegions][2] =
    {{{0.100,0.180}, {0.030, 0.070}, {0.209, 0.249}},
    {{0.456,0.683}, {0.257, 0.371}, {0.767, 0.880}}};*/

    // 1.5 sigma
    /*float band_bounds[nParticles][nRegions][2] =
    {{{0.110,0.170}, {0.030, 0.070}, {0.209, 0.249}},
    {{0.484,0.654}, {0.257, 0.371}, {0.767, 0.880}}};*/

    h_pair_mass->Scale(scaling_factor);
    // Compute ratios
    int firstbin, secondbin;
    firstbin = h_pair_mass->FindBin(band_bounds[0][0][0]);
    secondbin = h_pair_mass->FindBin(band_bounds[0][0][1]);
    float pi0_ratio = fit_global_bkg_1->Integral(band_bounds[0][0][0], band_bounds[0][0][1]) /
      fit_global_1->Integral(band_bounds[0][0][0], band_bounds[0][0][1]);
      // h_pair_mass->Integral(firstbin, secondbin) * (secondbin - firstbin) /
      // (band_bounds[0][0][1] - band_bounds[0][0][0]);
    float pi0_ratio_1 = fit_global_bkg_1->Integral(band_bounds[0][0][0] + 0.002, band_bounds[0][0][1] + 0.002) /
      fit_global_1->Integral(band_bounds[0][0][0] + 0.002, band_bounds[0][0][1] + 0.002);
      // h_pair_mass->Integral(firstbin+1, secondbin+1) * (secondbin - firstbin) /
      // (band_bounds[0][0][1] - band_bounds[0][0][0]);
    float pi0_ratio_2 = fit_global_bkg_1->Integral(band_bounds[0][0][0] - 0.002, band_bounds[0][0][1] - 0.002) /
      fit_global_1->Integral(band_bounds[0][0][0] - 0.002, band_bounds[0][0][1] - 0.002);
      // h_pair_mass->Integral(firstbin-1, secondbin-1) * (secondbin - firstbin) /
      // (band_bounds[0][0][1] - band_bounds[0][0][0]);
    float pi0_ratio_err = (std::abs(pi0_ratio_1 - pi0_ratio) + std::abs(pi0_ratio_2 - pi0_ratio)) / 2;
    //fit_global_1->Integral(band_bounds[0][0][0], band_bounds[0][0][1]);
    float pi0_signal_range = h_pair_mass->Integral(firstbin, secondbin);

    firstbin = h_pair_mass->FindBin(band_bounds[1][0][0]);
    secondbin = h_pair_mass->FindBin(band_bounds[1][0][1]);
    float eta_ratio = fit_global_bkg_2->Integral(band_bounds[1][0][0], band_bounds[1][0][1]) /
      fit_global_2->Integral(band_bounds[1][0][0], band_bounds[1][0][1]);
      // h_pair_mass->Integral(firstbin, secondbin) * (secondbin - firstbin) /
      // (band_bounds[1][0][1] - band_bounds[1][0][0]);
    float eta_ratio_1 = fit_global_bkg_2->Integral(band_bounds[1][0][0] + 0.002, band_bounds[1][0][1] + 0.002) /
      fit_global_2->Integral(band_bounds[1][0][0] + 0.002, band_bounds[1][0][1] + 0.002);
      // h_pair_mass->Integral(firstbin + 1, secondbin + 1) * (secondbin - firstbin) /
      // (band_bounds[1][0][1] - band_bounds[1][0][0]);
    float eta_ratio_2 = fit_global_bkg_2->Integral(band_bounds[1][0][0] - 0.002, band_bounds[1][0][1] - 0.002) /
      fit_global_2->Integral(band_bounds[1][0][0] - 0.002, band_bounds[1][0][1] - 0.002);
      // h_pair_mass->Integral(firstbin + 1, secondbin + 1) * (secondbin - firstbin) /
      // (band_bounds[1][0][1] - band_bounds[1][0][0]);
    float eta_ratio_err = (std::abs(eta_ratio_1 - eta_ratio) + std::abs(eta_ratio_1 - eta_ratio)) / 2;
    //fit_global_2->Integral(band_bounds[1][0][0], band_bounds[1][0][1]);
    float eta_signal_range = h_pair_mass->Integral(firstbin, secondbin);
    

    // Access easily the total number of pi0 and eta
    if (iPt == 9) {
      float pi0_signal_stat = pi0_signal_range * (1 - pi0_ratio);
      float eta_signal_stat = eta_signal_range * (1 - eta_ratio);
      std::cout << "Total Number of entries = " << h_pair_mass->GetEntries() << std::endl;
      std::cout << "Total Number Integral = " << h_pair_mass->Integral() << std::endl;
      std::cout << "Total Number of pi0/eta mesons = " << pi0_signal_stat << "/" << eta_signal_stat << std::endl;
    }

    
    if (iPt < nPtBins) {
      //if (iPt > 1 && iPt < nPtBins - 1) {
      // Skip the first bins in pT (not included in the sPHENIX preliminary release)
      if (pi0_ratio > 0 && pi0_ratio < 1) {
        pt_ratios_pi0->SetPoint(iPt, pTMeans[0][iPt], pi0_ratio * 100);
        pt_ratios_pi0->SetPointError(iPt, 0, pi0_ratio_err * 100);
      }
      else {
        pt_ratios_pi0->SetPoint(iPt, pTMeans[0][iPt], 100);
        pt_ratios_pi0->SetPointError(iPt, 0, 0);
      }
      if (eta_ratio > 0 && eta_ratio < 1) {
        pt_ratios_eta->SetPoint(iPt, pTMeans[0][iPt], eta_ratio * 100);
        pt_ratios_eta->SetPointError(iPt, 0, eta_ratio_err * 100);
      }
      else {
        pt_ratios_eta->SetPoint(iPt, pTMeans[0][iPt], 100);
        pt_ratios_eta->SetPointError(iPt, 0, 0);
      }
    }

    std::ofstream ratioList;
    ratioList.open((plots_folder + "/ratioList.txt"), ios::out | ios::app);
    ratioList << "bin : p_T in " << pTBins[iPt] << " GeV/c: r_{pi0} = " << pi0_ratio << ", r_{eta} = " << eta_ratio << std::endl;
    ratioList.close();

    h_pair_mass->SetMinimum(0);
    h_pair_mass->GetYaxis()->SetTitle("counts / [2 MeV]");
    if (iPt == 8){
      double currentMax = h_pair_mass->GetMaximum();
      h_pair_mass->SetMaximum(3 * currentMax);
      h_pair_mass->SetMinimum(0);
    }
    if (iPt == 7){
      double currentMax = h_pair_mass->GetMaximum();
      h_pair_mass->SetMaximum(1.5 * currentMax);
      h_pair_mass->SetMinimum(0);
    }

    //if (iPt == 2) h_pair_mass->SetMaximum(800000);

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

    // Store the fit results
    std::string outfilename_root = plots_folder + "/parametric_fit_mass_pt_" + std::to_string(iPt) + ".root";
    TFile *outfile_root = new TFile(outfilename_root.c_str(), "RECREATE");
    outfile_root->cd();

    h_pair_mass->SetName("h_pair_mass");
    h_pair_mass->Write();

    fit_global_bkg_1->SetName("fit_global_bkg_pi0");
    fit_global_bkg_1->SetRange(offset, borderFitMin);
    fit_global_bkg_1->SetNpx(1000);
    fit_global_bkg_1->SetLineWidth(2);
    fit_global_bkg_1->SetLineColor(kGreen);
    fit_global_bkg_1->SetLineStyle(2);
    fit_global_bkg_1->Write();

    fit_global_bkg_2->SetName("fit_global_bkg_eta");
    fit_global_bkg_2->SetRange(borderFitMax, 1.0);
    fit_global_bkg_2->SetNpx(1000);
    fit_global_bkg_2->SetLineWidth(2);
    fit_global_bkg_2->SetLineColor(kGreen);
    fit_global_bkg_2->SetLineStyle(2);
    fit_global_bkg_2->Write();

    fit_global_pi0->SetName("fit_global_pi0");
    fit_global_pi0->SetRange(offset, 1.0);
    fit_global_pi0->SetNpx(1000);
    fit_global_pi0->SetLineWidth(2);
    fit_global_pi0->SetLineColor(kBlue);
    fit_global_pi0->SetLineStyle(kSolid);
    fit_global_pi0->Write();
    
    fit_global_eta->SetName("fit_global_eta");
    fit_global_eta->SetRange(offset, 1.0);
    fit_global_eta->SetNpx(1000);
    fit_global_eta->SetLineWidth(2);
    fit_global_eta->SetLineColor(kBlue);
    fit_global_eta->SetLineStyle(kSolid);
    fit_global_eta->Write();
    outfile_root->Close();
    delete outfile_root;
    
    // Print the histogram and the fits
    TFile *outfile = new TFile("fits.root", "RECREATE");
    std::string canvas_name = "canvas_pi0_eta_mass_pt_" + pTBins[iPt];
    TCanvas *canvas = new TCanvas(canvas_name.c_str(), "c", 1600, 900);
    
    // Draw original histogram
    h_pair_mass->SetLineColor(kBlack);
    h_pair_mass->SetStats(0);
    h_pair_mass->Draw("E");

    // Draw fits
    /*fit_bkg_1->SetRange(offset, borderFit);
    fit_bkg_1->SetNpx(1000);
    fit_bkg_1->SetLineWidth(3);
    fit_bkg_1->SetLineColor(kGreen);
    fit_bkg_1->SetLineStyle(kDashed);
    fit_bkg_1->Draw("same L");

    fit_bkg_2->SetRange(borderFit,1);
    fit_bkg_2->SetNpx(1000);
    fit_bkg_2->SetLineWidth(3);
    fit_bkg_2->SetLineColor(kGreen);
    fit_bkg_2->SetLineStyle(kDashed);
    fit_bkg_2->Draw("same L");

    fit_pi0->SetRange(offset, borderFit);
    fit_pi0->SetNpx(1000);
    fit_pi0->SetLineWidth(3);
    fit_pi0->SetLineColor(kBlue);
    fit_pi0->SetLineStyle(kDashed);
    fit_pi0->Draw("same L");

    fit_eta->SetRange(borderFit, 1.0);
    fit_eta->SetNpx(1000);
    fit_eta->SetLineWidth(3);
    fit_eta->SetLineColor(kYellow);
    fit_eta->SetLineStyle(kSolid);
    fit_eta->Draw("same L");*/

    fit_global_bkg_1->SetRange(offset, borderFitMin);
    fit_global_bkg_1->SetNpx(1000);
    fit_global_bkg_1->SetLineWidth(2);
    fit_global_bkg_1->SetLineColor(kGreen);
    fit_global_bkg_1->SetLineStyle(2);
    fit_global_bkg_1->Draw("same L");

    fit_global_bkg_2->SetRange(borderFitMax, 1.0);
    fit_global_bkg_2->SetNpx(1000);
    fit_global_bkg_2->SetLineWidth(2);
    fit_global_bkg_2->SetLineColor(kGreen);
    fit_global_bkg_2->SetLineStyle(2);
    fit_global_bkg_2->Draw("same L");

    fit_global_pi0->SetRange(offset, 1.0);
    fit_global_pi0->SetNpx(1000);
    fit_global_pi0->SetLineWidth(2);
    fit_global_pi0->SetLineColor(kBlue);
    fit_global_pi0->SetLineStyle(kSolid);
    fit_global_pi0->Draw("same L");

    fit_global_eta->SetRange(offset, 1.0);
    fit_global_eta->SetNpx(1000);
    fit_global_eta->SetLineWidth(2);
    fit_global_eta->SetLineColor(kBlue);
    fit_global_eta->SetLineStyle(kSolid);
    fit_global_eta->Draw("same L");
    
    fit_global_1->SetRange(offset, borderFitMin);
    fit_global_1->SetNpx(1000);
    fit_global_1->SetLineWidth(2);
    fit_global_1->SetLineColor(kRed);
    fit_global_1->SetLineStyle(kSolid);
    fit_global_1->Draw("same L");

    fit_global_2->SetRange(borderFitMax, 1.0);
    fit_global_2->SetNpx(1000);
    fit_global_2->SetLineWidth(2);
    fit_global_2->SetLineColor(kRed);
    fit_global_2->SetLineStyle(kSolid);
    fit_global_2->Draw("same L");

    // Draw peak regions
    pi0PeakBand->Draw("same");
    pi0SideBandLeft->Draw("same");
    pi0SideBandRight->Draw("same");
    etaPeakBand->Draw("same");
    etaSideBandLeft->Draw("same");
    etaSideBandRight->Draw("same");

    // Show relevant numerical information
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
    //std::cout << "text size: " << latex.GetTextSize() << std::endl;
    //latex.SetTextSize(18);
    if (iPt < nPtBins) {
      latex.DrawLatex(0.47,0.87, (pTBounds[iPt][0] + " < p_{T} [GeV] < " + pTBounds[iPt][1]).c_str());
    }
    /*latex.DrawLatex(0.5, 0.80,"#pi^{0}");
    stream << "mean: " << std::fixed << std::setprecision(5) << fit_global_pi0->GetParameter(1) << " #pm " << fit_global_1->GetParError(5);
    latex.DrawLatex(0.5, 0.75, stream.str().c_str());
    stream.str(""); stream << "#sigma: " << std::setprecision(5) << std::abs(fit_global_pi0->GetParameter(2)) << " #pm " << fit_global_1->GetParError(6);
    latex.DrawLatex(0.5, 0.70, stream.str().c_str());
    stream.str(""); stream << "r: " << std::setprecision(5) << pi0_ratio  << " #pm " << pi0_ratio_err;
    latex.DrawLatex(0.5, 0.65, stream.str().c_str());
    latex.DrawLatex(0.5, 0.55,"#eta");
    stream.str(""); stream << "mean: " << std::setprecision(5) << fit_global_eta->GetParameter(1) << " #pm " << fit_global_2->GetParError(5);
    latex.DrawLatex(0.5, 0.50, stream.str().c_str());
    stream.str(""); stream << "#sigma: " << std::setprecision(5) << std::abs(fit_global_eta->GetParameter(2)) << " #pm " << fit_global_2->GetParError(6);
    latex.DrawLatex(0.5, 0.45, stream.str().c_str());
    stream.str(""); stream << "r: " << std::setprecision(5) << eta_ratio << " #pm " << eta_ratio_err;
    latex.DrawLatex(0.5, 0.40, stream.str().c_str());*/

    // Add the full legend
    TLegend *full_legend = new TLegend(0.55, 0.5, 0.85, 0.70);
    full_legend->SetTextSize(0.05);
    full_legend->AddEntry(h_pair_mass, "Data");
    full_legend->AddEntry(fit_global_1, "Full Fit", "l");
    full_legend->AddEntry(fit_global_bkg_1, "Fitted Background", "l");
    full_legend->AddEntry(fit_global_pi0, "Fitted Signal", "l");
    full_legend->Draw();
    
    // Save histograms
    canvas->SaveAs((plots_folder + "/" + canvas_name + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name + ".png").c_str());

    // Close file
    outfile->cd();
    canvas->Write();
    outfile->Close();
    delete outfile;

    delete fitMethodBkg1;
    delete fitMethodBkg2;
    delete fitMethodPi0;
    delete fitMethodEta;
    delete fitMethodGlobal1;
    delete fitMethodGlobal2;
    delete fitMethodGlobalBkg1;
    delete fitMethodGlobalBkg2;
    delete fitMethodGlobalPi0;
    delete fitMethodGlobalEta;
  }

  // Draw background ratios
  std::string canvas_name = "canvas_pt_background_ratios";
  TCanvas *canvas = new TCanvas(canvas_name.c_str(), "c", 1600, 900);
  TMultiGraph *mg = new TMultiGraph();
  mg->GetXaxis()->SetTitle("p_{T} [GeV]");
  mg->GetYaxis()->SetTitle("Background Fraction [%]");
  pt_ratios_pi0->SetMarkerSize(2.5);
  pt_ratios_pi0->SetMarkerStyle(kFullCircle);
  pt_ratios_pi0->SetMarkerColor(kBlue);
  pt_ratios_eta->SetMarkerSize(2.5);
  pt_ratios_eta->SetMarkerStyle(kFullSquare);
  pt_ratios_eta->SetMarkerColor(kRed);
  mg->Add(pt_ratios_pi0);
  mg->Add(pt_ratios_eta);
  mg->SetMinimum(0);
  mg->SetMaximum(100);
  mg->GetXaxis()->SetRangeUser(2, 20);
  mg->Draw("AP");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextColor(kBlack);
  latex.DrawLatex(0.7, 0.75, "#font[72]{sPHENIX} Internal");
  latex.DrawLatex(0.7, 0.65, "p^{#uparrow}+p #sqrt{s} = 200 GeV");

  TLegend *legend = new TLegend(0.2, 0.5, 0.45, 0.7);
  //TLegend *legend = new TLegend(0.6, 0.2, 0.9, 0.4);
  legend->AddEntry(pt_ratios_pi0, "p^{#uparrow}p #rightarrow #pi^{0} X");
  legend->AddEntry(pt_ratios_eta, "p^{#uparrow}p #rightarrow #eta X");
  legend->SetBorderSize(1);
  legend->Draw();

  canvas->SaveAs((plots_folder + "/" + canvas_name + ".pdf").c_str());
  canvas->SaveAs((plots_folder + "/" + canvas_name + ".png").c_str());
  
  canvas_name = "canvas_pt_background_ratios_pi0";
  canvas = new TCanvas("canvas_pt_background_ratios_pi0", "c", 1600, 900);
  canvas->cd();
  pt_ratios_pi0->SetMinimum(0);
  pt_ratios_pi0->SetMaximum(100);
  //pt_ratios_pi0->SetLineWidth(0);
  pt_ratios_pi0->SetMarkerSize(1.5);
  pt_ratios_pi0->SetMarkerStyle(kFullCircle);
  pt_ratios_pi0->SetMarkerColor(kBlue);
  pt_ratios_pi0->Draw("AP");

  canvas->SaveAs((plots_folder + "/" + canvas_name + ".pdf").c_str());
  canvas->SaveAs((plots_folder + "/" + canvas_name + ".png").c_str());
  delete canvas;

  canvas_name = "canvas_pt_background_ratios_eta";
  canvas = new TCanvas("canvas_pt_background_ratios_eta", "c", 1600, 900);
  canvas->cd();
  pt_ratios_eta->SetMinimum(70);
  pt_ratios_eta->SetMaximum(100);
  //pt_ratios_eta->SetLineWidth(0);
  pt_ratios_eta->SetMarkerSize(1.5);
  pt_ratios_eta->SetMarkerStyle(kFullCircle);
  pt_ratios_eta->SetMarkerColor(kRed);
  pt_ratios_eta->Draw();

  std::ofstream ratio_csv;
  ratio_csv.open((outputfolder_ratios + "/" + production_date + "_" + selection_label + ".csv"));
  for (int iPt = 0; iPt < nPtBins - 1; iPt++) ratio_csv << pt_ratios_pi0->GetPointY(iPt) / 100. << ", ";
  ratio_csv << pt_ratios_pi0->GetPointY(nPtBins - 1) / 100. << "\n";
  for (int iPt = 0; iPt < nPtBins - 1; iPt++) ratio_csv << pt_ratios_eta->GetPointY(iPt) / 100. << ", ";
  ratio_csv << pt_ratios_eta->GetPointY(nPtBins - 1) / 100. << "\n";
  for (int iPt = 0; iPt < nPtBins - 1; iPt++) ratio_csv << 0  << ", ";
  ratio_csv << 0 << "\n";
  for (int iPt = 0; iPt < nPtBins - 1; iPt++) ratio_csv << 0 << ", ";
  ratio_csv << 0 << "\n";
  ratio_csv.close();

  canvas->SaveAs((plots_folder + "/" + canvas_name + ".pdf").c_str());
  canvas->SaveAs((plots_folder + "/" + canvas_name + ".png").c_str());
  delete canvas;
  
  gSystem->Exit(0);
}
