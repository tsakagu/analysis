void fit_symmetric_graph(TGraphErrors* (&graph_double),
                         const TGraphErrors* graph_oneside)
{
  graph_double = new TGraphErrors();
  std::stringstream graphname;
  graphname << graph_oneside->GetName() << "_double";
  graph_double->SetName(graphname.str().c_str());
  const int N = graph_oneside->GetN() * 2;
  for (int i = 0; i < N/2; i++) {
    double x = graph_oneside->GetPointX(i);
    double y = graph_oneside->GetPointY(i);
    double ey = graph_oneside->GetErrorY(i);
    graph_double->SetPoint(i, x, y); graph_double->SetPoint(N-1-i, -x, y);
    graph_double->SetPointError(i, 0, ey); graph_double->SetPointError(N-1-i, 0, ey);
  }
}

void write_csv_ratios(const std::string &selection_label = "MBD_pt_l3_0mrad",
                      const std::string &production_date = "01312026",
                      const std::string &inputfolder_gpr_folder = "",
                      const std::string &outputfolder = "csv_ratios")
{ 
  const std::string inputfolder = inputfolder_gpr_folder + "/result_comparison/analysis_ana509/01312026/" + selection_label;
  
  const int nPtBins = 9;
  const int nEtaBins = 8;
  const int nXfBins = 8;

  // pT asymmetries
  {
    std::string inputfilename = inputfolder + "/pt_ratios.root";
    TFile *inputfile_ratios = TFile::Open(inputfilename.c_str());
    
    std::stringstream csv_parametric_name;
    csv_parametric_name << outputfolder << "/pt_ratios_" << production_date << "_" << selection_label << ".csv";
    std::stringstream csv_gpr_name;
    csv_gpr_name << outputfolder << "/pt_ratios_" << production_date << "_" << selection_label << "_gpr.csv";

    TGraphErrors *pt_parametric_ratios_pi0 = (TGraphErrors*) inputfile_ratios->Get("pt_background_parametric_ratios_pi0");
    TGraphErrors *pt_parametric_ratios_eta = (TGraphErrors*) inputfile_ratios->Get("pt_background_parametric_ratios_eta");
    TGraphErrors *pt_gpr_ratios_pi0 = (TGraphErrors*) inputfile_ratios->Get("pt_background_gpr_ratios_pi0");
    TGraphErrors *pt_gpr_ratios_eta = (TGraphErrors*) inputfile_ratios->Get("pt_background_gpr_ratios_eta");

    std::ofstream fout_parametric(csv_parametric_name.str());
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      double ratio = (pt_parametric_ratios_pi0->GetPointY(iPt) / 100);
      fout_parametric << (ratio == ratio ? ratio : 1) << ", ";
    }
    fout_parametric << (pt_parametric_ratios_pi0->GetPointY(nPtBins-1) / 100) << "\n";
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      double ratio = (pt_parametric_ratios_eta->GetPointY(iPt) / 100);
      fout_parametric << (ratio == ratio ? ratio : 1) << ", ";
    }
    fout_parametric << (pt_parametric_ratios_eta->GetPointY(nPtBins-1) / 100) << "\n";
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      fout_parametric << 0 << ", ";
    }
    fout_parametric << 0 << "\n";
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      fout_parametric << 0 << ", ";
    }
    fout_parametric << 0 << "\n";
    fout_parametric.close();

    std::ofstream fout_gpr(csv_gpr_name.str());
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      double ratio = (pt_gpr_ratios_pi0->GetPointY(iPt) / 100);
      fout_gpr << (ratio == ratio ? ratio : 1) << ", ";
    }
    fout_gpr << (pt_gpr_ratios_pi0->GetPointY(nPtBins-1) / 100) << "\n";
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      double ratio = (pt_gpr_ratios_eta->GetPointY(iPt) / 100);
      fout_gpr << (ratio == ratio ? ratio : 1) << ", ";
    }
    fout_gpr << (pt_gpr_ratios_eta->GetPointY(nPtBins-1) / 100) << "\n";
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      fout_gpr << 0 << ", ";
    }
    fout_gpr << 0 << "\n";
    for (int iPt = 0; iPt < nPtBins - 1; iPt++) {
      fout_gpr << 0 << ", ";
    }
    fout_gpr << 0 << "\n";
    fout_gpr.close();
  }

  // xF asymmetries
  {
    std::string inputfilename = inputfolder + "/xf_ratios.root";
    TFile *inputfile_ratios = TFile::Open(inputfilename.c_str());
    
    std::stringstream csv_parametric_name;
    csv_parametric_name << outputfolder << "/xf_ratios_" << production_date << "_" << selection_label << ".csv";
    std::stringstream csv_gpr_name;
    csv_gpr_name << outputfolder << "/xf_ratios_" << production_date << "_" << selection_label << "_gpr.csv";

    TGraphErrors *xf_parametric_ratios_pi0_single = (TGraphErrors*) inputfile_ratios->Get("xf_background_parametric_ratios_pi0");
    TGraphErrors *xf_parametric_ratios_eta_single = (TGraphErrors*) inputfile_ratios->Get("xf_background_parametric_ratios_eta");
    TGraphErrors *xf_gpr_ratios_pi0_single = (TGraphErrors*) inputfile_ratios->Get("xf_background_gpr_ratios_pi0");
    TGraphErrors *xf_gpr_ratios_eta_single = (TGraphErrors*) inputfile_ratios->Get("xf_background_gpr_ratios_eta");

    TGraphErrors *xf_parametric_ratios_pi0 = nullptr; fit_symmetric_graph(xf_parametric_ratios_pi0, xf_parametric_ratios_pi0_single);
    TGraphErrors *xf_parametric_ratios_eta = nullptr; fit_symmetric_graph(xf_parametric_ratios_eta, xf_parametric_ratios_eta_single);
    TGraphErrors *xf_gpr_ratios_pi0 = nullptr; fit_symmetric_graph(xf_gpr_ratios_pi0, xf_gpr_ratios_pi0_single);
    TGraphErrors *xf_gpr_ratios_eta = nullptr; fit_symmetric_graph(xf_gpr_ratios_eta, xf_gpr_ratios_eta_single);

    std::ofstream fout_parametric(csv_parametric_name.str());
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_parametric << (xf_parametric_ratios_pi0->GetPointY(iXf) / 100) << ", ";
    }
    fout_parametric << (xf_parametric_ratios_pi0->GetPointY(nXfBins-1) / 100) << "\n";
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_parametric << (xf_parametric_ratios_eta->GetPointY(iXf) / 100) << ", ";
    }
    fout_parametric << (xf_parametric_ratios_eta->GetPointY(nXfBins-1) / 100) << "\n";
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_parametric << 0 << ", ";
    }
    fout_parametric << 0 << "\n";
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_parametric << 0 << ", ";
    }
    fout_parametric << 0 << "\n";
    fout_parametric.close();

    std::ofstream fout_gpr(csv_gpr_name.str());
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_gpr << (xf_gpr_ratios_pi0->GetPointY(iXf) / 100) << ", ";
    }
    fout_gpr << (xf_gpr_ratios_pi0->GetPointY(nXfBins-1) / 100) << "\n";
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_gpr << (xf_gpr_ratios_eta->GetPointY(iXf) / 100) << ", ";
    }
    fout_gpr << (xf_gpr_ratios_eta->GetPointY(nXfBins-1) / 100) << "\n";
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_gpr << 0 << ", ";
    }
    fout_gpr << 0 << "\n";
    for (int iXf = 0; iXf < nXfBins - 1; iXf++) {
      fout_gpr << 0 << ", ";
    }
    fout_gpr << 0 << "\n";
    fout_gpr.close();
  }

  // eta asymmetries
  {
    std::string inputfilename = inputfolder + "/eta_ratios.root";
    TFile *inputfile_ratios = TFile::Open(inputfilename.c_str());
    
    std::stringstream csv_parametric_name;
    csv_parametric_name << outputfolder << "/eta_ratios_" << production_date << "_" << selection_label << ".csv";
    std::stringstream csv_gpr_name;
    csv_gpr_name << outputfolder << "/eta_ratios_" << production_date << "_" << selection_label << "_gpr.csv";

    TGraphErrors *eta_parametric_ratios_pi0_single = (TGraphErrors*) inputfile_ratios->Get("eta_background_parametric_ratios_pi0");
    TGraphErrors *eta_parametric_ratios_eta_single = (TGraphErrors*) inputfile_ratios->Get("eta_background_parametric_ratios_eta");
    TGraphErrors *eta_gpr_ratios_pi0_single = (TGraphErrors*) inputfile_ratios->Get("eta_background_gpr_ratios_pi0");
    TGraphErrors *eta_gpr_ratios_eta_single = (TGraphErrors*) inputfile_ratios->Get("eta_background_gpr_ratios_eta");

    TGraphErrors *eta_parametric_ratios_pi0 = nullptr; fit_symmetric_graph(eta_parametric_ratios_pi0, eta_parametric_ratios_pi0_single);
    TGraphErrors *eta_parametric_ratios_eta = nullptr; fit_symmetric_graph(eta_parametric_ratios_eta, eta_parametric_ratios_eta_single);
    TGraphErrors *eta_gpr_ratios_pi0 = nullptr; fit_symmetric_graph(eta_gpr_ratios_pi0, eta_gpr_ratios_pi0_single);
    TGraphErrors *eta_gpr_ratios_eta = nullptr; fit_symmetric_graph(eta_gpr_ratios_eta, eta_gpr_ratios_eta_single);

    std::ofstream fout_parametric(csv_parametric_name.str());
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_parametric << (eta_parametric_ratios_pi0->GetPointY(iEta) / 100) << ", ";
    }
    fout_parametric << (eta_parametric_ratios_pi0->GetPointY(nEtaBins-1) / 100) << "\n";
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_parametric << (eta_parametric_ratios_eta->GetPointY(iEta) / 100) << ", ";
    }
    fout_parametric << (eta_parametric_ratios_eta->GetPointY(nEtaBins-1) / 100) << "\n";
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_parametric << 0 << ", ";
    }
    fout_parametric << 0 << "\n";
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_parametric << 0 << ", ";
    }
    fout_parametric << 0 << "\n";
    fout_parametric.close();

    std::ofstream fout_gpr(csv_gpr_name.str());
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_gpr << (eta_gpr_ratios_pi0->GetPointY(iEta) / 100) << ", ";
    }
    fout_gpr << (eta_gpr_ratios_pi0->GetPointY(nEtaBins-1) / 100) << "\n";
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_gpr << (eta_gpr_ratios_eta->GetPointY(iEta) / 100) << ", ";
    }
    fout_gpr << (eta_gpr_ratios_eta->GetPointY(nEtaBins-1) / 100) << "\n";
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_gpr << 0 << ", ";
    }
    fout_gpr << 0 << "\n";
    for (int iEta = 0; iEta < nEtaBins - 1; iEta++) {
      fout_gpr << 0 << ", ";
    }
    fout_gpr << 0 << "\n";
    fout_gpr.close();
  }

  gSystem->Exit(0);
}
