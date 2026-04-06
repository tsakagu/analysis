void read_average_bin(const std::string &production_date = "01312026",
                      const std::string &selection_label = "MBD_INCOMPLETE",
                      const std::string &inputfolder_mass = "",
                      const std::string &inputfolder_analysis = "",
                      const std::string &outputfolder_csv = "")
{
  std::string inputfilename = inputfolder_analysis + "/analysis_complete_ana509_" + production_date + "_" + selection_label + ".root";
  std::string output_csv_pt_name = outputfolder_csv + "/input_pt_mean_" + production_date + "_" + selection_label + ".csv";
  std::string output_csv_eta_name = outputfolder_csv + "/input_eta_mean_" + production_date + "_" + selection_label + ".csv";
  std::string output_csv_xf_name = outputfolder_csv + "/input_xf_mean_" + production_date + "_" + selection_label + ".csv";

  const int nPtBins = 9;
  const int nEtaBins = 8;
  const int nXfBins = 8;

  TFile *inputfile = TFile::Open(inputfilename.c_str());

  TH1F *h_average_pt_pi0 = (TH1F*) inputfile->Get("h_average_pt_pi0");
  TH1F *h_average_pt_eta = (TH1F*) inputfile->Get("h_average_pt_eta");
  TH1F *h_norm_pt_pi0 = (TH1F*) inputfile->Get("h_norm_pt_pi0");
  TH1F *h_norm_pt_eta = (TH1F*) inputfile->Get("h_norm_pt_eta");
  TH1F *h_average_eta_pi0 = (TH1F*) inputfile->Get("h_average_eta_pi0");
  TH1F *h_average_eta_eta = (TH1F*) inputfile->Get("h_average_eta_eta");
  TH1F *h_norm_eta_pi0 = (TH1F*) inputfile->Get("h_norm_eta_pi0");
  TH1F *h_norm_eta_eta = (TH1F*) inputfile->Get("h_norm_eta_eta");
  TH1F *h_average_xf_pi0 = (TH1F*) inputfile->Get("h_average_xf_pi0");
  TH1F *h_average_xf_eta = (TH1F*) inputfile->Get("h_average_xf_eta");
  TH1F *h_norm_xf_pi0 = (TH1F*) inputfile->Get("h_norm_xf_pi0");
  TH1F *h_norm_xf_eta = (TH1F*) inputfile->Get("h_norm_xf_eta");

  std::ofstream output_csv_pt;
  output_csv_pt.open(output_csv_pt_name);
  for (int i=1; i< nPtBins; i++) {
    double ratio = h_average_pt_pi0->GetBinContent(i) / h_norm_pt_pi0->GetBinContent(i);
    output_csv_pt << ratio << ", ";
  }
  output_csv_pt << h_average_pt_pi0->GetBinContent(nPtBins) / h_norm_pt_pi0->GetBinContent(nPtBins) << "\n";
  for (int i=1; i< nPtBins; i++) {
    double ratio = h_average_pt_eta->GetBinContent(i) / h_norm_pt_eta->GetBinContent(i);
    output_csv_pt << ratio << ", ";
  }
  output_csv_pt << h_average_pt_eta->GetBinContent(nPtBins) / h_norm_pt_eta->GetBinContent(nPtBins) << "\n";
  output_csv_pt.close();

  // Make the bins symmetric
  // std::ofstream output_csv_eta;
  // output_csv_eta.open(output_csv_eta_name);
  // for (int i=1; i< nEtaBins; i++) {
  //   double ratio = h_average_eta_pi0->GetBinContent(i) / h_norm_eta_pi0->GetBinContent(i);
  //   output_csv_eta << ratio << ", ";
  // }
  // output_csv_eta << (h_average_eta_pi0->GetBinContent(nEtaBins) / h_norm_eta_pi0->GetBinContent(nEtaBins)) << "\n";
  // for (int i=1; i< nEtaBins; i++) {
  //   double ratio = h_average_eta_eta->GetBinContent(i) / h_norm_eta_eta->GetBinContent(i);
  //   output_csv_eta << ratio << ", ";
  // }
  // output_csv_eta << (h_average_eta_eta->GetBinContent(nEtaBins) / h_norm_eta_eta->GetBinContent(nEtaBins)) << "\n";
  // output_csv_eta.close();

  // std::ofstream output_csv_xf;
  // output_csv_xf.open(output_csv_xf_name);
  // for (int i=1; i< nXfBins; i++) {
  //   double ratio = h_average_xf_pi0->GetBinContent(i) / h_norm_xf_pi0->GetBinContent(i);
  //   output_csv_xf << ratio << ", ";
  // }
  // output_csv_xf << (h_average_xf_pi0->GetBinContent(nXfBins) / h_norm_xf_pi0->GetBinContent(nXfBins)) << "\n";
  // for (int i=1; i< nXfBins; i++) {
  //   double ratio = h_average_xf_eta->GetBinContent(i) / h_norm_xf_eta->GetBinContent(i);
  //   output_csv_xf << ratio << ", ";
  // }
  // output_csv_xf << (h_average_xf_eta->GetBinContent(nXfBins) / h_norm_xf_eta->GetBinContent(nXfBins)) << "\n";
  // output_csv_xf.close();


  std::ofstream output_csv_eta;
  output_csv_eta.open(output_csv_eta_name);
  for (int i=1; i< nEtaBins; i++) {
    double ratio = (h_average_eta_pi0->GetBinContent(i) - h_average_eta_pi0->GetBinContent(nEtaBins+1-i)) / (h_norm_eta_pi0->GetBinContent(i) + h_norm_eta_pi0->GetBinContent(nEtaBins+1-i));
    output_csv_eta << ratio << ", ";
  }
  output_csv_eta << ((h_average_eta_pi0->GetBinContent(nEtaBins) - h_average_eta_pi0->GetBinContent(1)) / (h_norm_eta_pi0->GetBinContent(nEtaBins) + h_norm_eta_pi0->GetBinContent(1))) << "\n";
  for (int i=1; i< nEtaBins; i++) {
    double ratio = (h_average_eta_eta->GetBinContent(i) - h_average_eta_eta->GetBinContent(nEtaBins+1-i)) / (h_norm_eta_eta->GetBinContent(i) + h_norm_eta_eta->GetBinContent(nEtaBins+1-i));
    output_csv_eta << ratio << ", ";
  }
  output_csv_eta << ((h_average_eta_eta->GetBinContent(nEtaBins) - h_average_eta_eta->GetBinContent(1)) / (h_norm_eta_eta->GetBinContent(nEtaBins) + h_norm_eta_eta->GetBinContent(1))) << "\n";
  output_csv_eta.close();

  std::ofstream output_csv_xf;
  output_csv_xf.open(output_csv_xf_name);
  for (int i=1; i< nXfBins; i++) {
    double ratio = (h_average_xf_pi0->GetBinContent(i) - h_average_xf_pi0->GetBinContent(nXfBins+1-i)) / (h_norm_xf_pi0->GetBinContent(i) + h_norm_xf_pi0->GetBinContent(nXfBins+1-i));
    output_csv_xf << ratio << ", ";
  }
  output_csv_xf << ((h_average_xf_pi0->GetBinContent(nXfBins) - h_average_xf_pi0->GetBinContent(1)) / (h_norm_xf_pi0->GetBinContent(nXfBins) + h_norm_xf_pi0->GetBinContent(1))) << "\n";
  for (int i=1; i< nXfBins; i++) {
    double ratio = (h_average_xf_eta->GetBinContent(i) - h_average_xf_eta->GetBinContent(nXfBins+1-i)) / (h_norm_xf_eta->GetBinContent(i) + h_norm_xf_eta->GetBinContent(nXfBins+1-i));
    output_csv_xf << ratio << ", ";
  }
  output_csv_xf << ((h_average_xf_eta->GetBinContent(nXfBins) - h_average_xf_eta->GetBinContent(1)) / (h_norm_xf_eta->GetBinContent(nXfBins) + h_norm_xf_eta->GetBinContent(1))) << "\n";
  output_csv_xf.close();
  
  gSystem->Exit(0);
}
