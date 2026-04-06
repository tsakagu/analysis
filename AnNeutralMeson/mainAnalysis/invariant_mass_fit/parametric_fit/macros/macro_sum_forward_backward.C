void macro_sum_forward_backward(const std::string &production_date = "01312026",
                                const std::string &selection_label = "MBD_INCOMPLETE",
                                const std::string &inputfolder_analysis = "")
{
  std::string filename = inputfolder_analysis + "/analysis_complete_ana509_" + production_date + "_" + selection_label + ".root";
  TFile *infile = TFile::Open(filename.c_str());
  std::string outfilename = inputfolder_analysis + "/analysis_complete_ana509_" + production_date + "_" + selection_label + "_symmetric.root";
  TFile *outfile = new TFile(outfilename.c_str(), "RECREATE");
  outfile->cd();

  const int nEtaBins = 8;
  const int nXfBins = 8;
  for (int iEta = 0; iEta < nEtaBins/2; iEta++) {
    std::string hname_forward = "h_pair_mass_eta_" + std::to_string(iEta);
    std::string hname_backward = "h_pair_mass_eta_" + std::to_string(nEtaBins-1-iEta);
    std::string hname_sym = "h_pair_mass_eta_" + std::to_string(iEta) + "_sym";
    TH1F *h_pair_mass_sym = new TH1F(hname_sym.c_str(),
                                 ";M_{#gamma#gamma} [GeV]; Counts",
                                 500, 0.0, 1.0);

    TH1F *h_pair_mass_forward = (TH1F*)infile->Get(hname_forward.c_str());
    TH1F *h_pair_mass_backward = (TH1F*)infile->Get(hname_backward.c_str());
    
    const int nBins = 500;
    for (int iBin = 1; iBin <= 500; iBin++) {
      h_pair_mass_sym->SetBinContent(iBin,
                                     h_pair_mass_forward->GetBinContent(iBin) +
                                     h_pair_mass_backward->GetBinContent(iBin));
    }
    h_pair_mass_sym->Write();
  }
  for (int iXf = 0; iXf < nXfBins/2; iXf++) {
    std::string hname_forward = "h_pair_mass_xf_" + std::to_string(iXf);
    std::string hname_backward = "h_pair_mass_xf_" + std::to_string(nXfBins-1-iXf);
    std::string hname_sym = "h_pair_mass_xf_" + std::to_string(iXf) + "_sym";
    TH1F *h_pair_mass_sym = new TH1F(hname_sym.c_str(),
                                 ";M_{#gamma#gamma} [GeV]; Counts",
                                 500, 0.0, 1.0);

    TH1F *h_pair_mass_forward = (TH1F*)infile->Get(hname_forward.c_str());
    TH1F *h_pair_mass_backward = (TH1F*)infile->Get(hname_backward.c_str());
    
    const int nBins = 500;
    for (int iBin = 1; iBin <= 500; iBin++) {
      h_pair_mass_sym->SetBinContent(iBin,
                                     h_pair_mass_forward->GetBinContent(iBin) +
                                     h_pair_mass_backward->GetBinContent(iBin));
    }
    h_pair_mass_sym->Write();
  }
  outfile->Close();
  delete outfile;
  gSystem->Exit(0);
}
