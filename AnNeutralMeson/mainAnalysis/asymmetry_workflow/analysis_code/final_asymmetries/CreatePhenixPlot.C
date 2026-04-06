#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

double GetAbsMaxInRange(TGraphErrors *graph, double xMin, double xMax)
{
  if (!graph ) return -1; // Error handling

  double maxVal = -1;
  int n = graph->GetN();

  for (int i = 0; i < n; i++)
  {
    double x, y, ey, yVal;
    graph->GetPoint(i, x, y);
    ey = graph->GetErrorY(i);
    yVal = std::max(std::abs(y-ey), std::abs(y+ey));
    if (x >= xMin && x <= xMax)
    {
      if (yVal > maxVal)
      {
        maxVal = yVal;
      }
    }
  }
  return maxVal;
}

void parse_theory_file(const std::string& filename,
                       std::vector<double>& x,
                       std::vector<double>& y,
                       std::vector<double>& ey)
{
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::string line;
  // Skip first line  (labels)
  std::getline(fin, line);
  // Read individual entries
  while (std::getline(fin, line)) {
    // Skip empty lines
    if (line.empty()) continue;

    // Replace '&' with spaces so stringstream can parse numbers
    std::replace(line.begin(), line.end(), '&', ' ');

    std::stringstream ss(line);
    double xv, yv, eyv;
    if (!(ss >> xv >> yv >> eyv)) {
      std::cerr << "Skipping malformed line: " << line << std::endl;
      continue;
    }

    x.push_back(xv);
    y.push_back(yv);
    ey.push_back(eyv);
  }
}

void parse_file(const std::string& filename,
                std::vector<double>& x,
                std::vector<double>& y,
                std::vector<double>& estat,
                std::vector<double>& esyst,
                std::vector<double>& exsyst,
                const double width = 0.15
                )
{  
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::string line;
  while (std::getline(fin, line)) {
    // Skip empty lines
    if (line.empty()) continue;

    // Replace '&' with spaces so stringstream can parse numbers
    std::replace(line.begin(), line.end(), '&', ' ');

    std::stringstream ss(line);
    double xv, yv, statv, systv;
    if (!(ss >> xv >> yv >> statv >> systv)) {
      std::cerr << "Skipping malformed line: " << line << std::endl;
      continue;
    }

    x.push_back(xv);
    y.push_back(yv);
    estat.push_back(statv);
    esyst.push_back(systv);

    // Constant half-width for the systematic boxes in x.
    // Choose something visually reasonable for your spacing.
    exsyst.push_back(width);
  }
}

void parse_phenix_file(const std::string& filename,
                       std::vector<double>& x,
                       std::vector<double>& y,
                       std::vector<double>& estat,
                       std::vector<double>& esyst,
                       std::vector<double>& exsyst,
                       const double width = 0.15
                )
{  
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::string line;
  // Skip first lines
  for (int i = 0; i < 8; i++) std::getline(fin, line);
  while (std::getline(fin, line)) {
    // Skip empty lines
    if (line.empty()) continue;

    // Replace ',' with spaces so stringstream can parse numbers
    std::replace(line.begin(), line.end(), ',', ' ');

    std::stringstream ss(line);
    double xv, yv, statv, systv, dummy;
    if (!(ss >> xv >> dummy >> dummy >> yv >> statv >> dummy >> dummy >> dummy >> dummy >> dummy >> systv)) {
      std::cerr << "Skipping malformed line: " << line << std::endl;
      continue;
    }

    x.push_back(xv);
    y.push_back(yv);
    estat.push_back(statv);
    esyst.push_back(systv);

    // Constant half-width for the systematic boxes in x.
    // Choose something visually reasonable for your spacing.
    exsyst.push_back(width);
  }
}

void CreatePhenixPlot()
{
  std::string plots_folder_pt = "pt_asymmetries/PHENIX/";

  gSystem->Exec(("mkdir -p " + plots_folder_pt).c_str());

  const std::string inputfolder = "asym_values/";
  const std::string inputfolder_phenix = "PHENIX21_AN/";
  const int nParticles = 2;
  const std::string particles[nParticles] = {"pi0", "eta"};

  float pTBounds[nParticles][2] = {{-0.026, 0.026}, {-0.09, 0.08}};
  
  for (int iP = 0; iP < nParticles; iP++) {
    std::stringstream inputfilename;
    inputfilename << inputfolder << "asym_summary_" << particles[iP] << "_pT" << ".txt";
    std::stringstream canvas_name;
    canvas_name << "canvas_asym_" << particles[iP] << "_pT";
    std::vector<double> x, y, estat, esyst, exsyst;
    parse_file(inputfilename.str(), x, y, estat, esyst, exsyst, 0.20);

    // Parse phenix plots
    std::stringstream th_filename;
    th_filename << inputfolder_phenix << "/phenix_asym_" << particles[iP] << ".csv";
    std::vector<double> x_ph, y_ph, estat_ph, esyst_ph, exsyst_ph;
    parse_phenix_file(th_filename.str(), x_ph, y_ph, estat_ph, esyst_ph, exsyst_ph, 0.20);

    const int n = x.size();
    auto gSyst = new TGraphErrors(n);
    auto gStat = new TGraphErrors(n);
    auto gSyst_ph = new TGraphErrors();
    auto gStat_ph = new TGraphErrors();
    int i_th = 0;
    for (int i = (iP == 0 ? 0 : 1); i < (iP == 0 ? n-1 : n); ++i) {
      gSyst->SetPoint(i, x[i], y[i]);
      gSyst->SetPointError(i, exsyst[i], esyst[i]);

      gStat->SetPoint(i, x[i], y[i]);
      gStat->SetPointError(i, 0.0, estat[i]);
    }
    const int n_ph = x_ph.size();
    for (int i = 0; i < n_ph; ++i) {
      gSyst_ph->SetPoint(i, x_ph[i], y_ph[i]);
      gSyst_ph->SetPointError(i, exsyst[0], esyst_ph[i]);

      gStat_ph->SetPoint(i, x_ph[i], y_ph[i]);
      gStat_ph->SetPointError(i, 0.0, estat_ph[i]);
    }

    SetsPhenixStyle();

    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "", 1600, 900);
    gStat->SetTitle(";p_{T} [GeV];A_{N}");

    double y_bound = 0;
    if (iP == 0){
      y_bound = GetAbsMaxInRange(gStat, 0.5, 10);
    } else {
      y_bound = GetAbsMaxInRange(gStat, 1.5, 12);
    }

    int color = (iP == 0 ? kRed : kBlue);

    gStat->SetLineWidth(2);
    gStat->SetLineColor(color);
    gStat->SetMarkerColor(color);
    gStat->SetMarkerStyle(kFullSquare);
    gStat->SetMarkerSize(2.4);
    // gStat->SetMinimum(-2 * y_bound);
    // gStat->SetMaximum(2 * y_bound);
    gStat->SetMinimum(pTBounds[iP][0]);
    gStat->SetMaximum(pTBounds[iP][1]);
    gStat->Draw("AP E1");

    // Style systematic boxes (transparent fill)
    if (iP == 0) {
      gSyst->SetFillColorAlpha(kPink-9, 0.45); // alpha controls transparency
      gSyst->SetLineColor(kPink-9);
    } else {
      gSyst->SetFillColorAlpha(kAzure-4, 0.45); // alpha controls transparency
      gSyst->SetLineColor(kAzure-4);
    }
    gSyst->SetMarkerStyle(0);
    gSyst->SetMarkerColor(color);
    gSyst->Draw("P E2 SAME");

    std::cout << "PHENIX stats:" << std::endl;
    gStat_ph->Print();
    std::cout << "sPHENIX stats:" << std::endl;
    gStat->Print();


    
    gSyst_ph->Print();
      
    gStat_ph->SetLineWidth(2);
    gStat_ph->SetLineColor(kGray+2);
    gStat_ph->SetMarkerColor(kGray+2);
    gStat_ph->SetMarkerStyle(59);
    gStat_ph->SetMarkerSize(2.6);
    gStat_ph->SetMinimum(-2 * y_bound);
    gStat_ph->SetMaximum(2 * y_bound);
    gStat_ph->Draw("P E1 SAME");

    gSyst_ph->SetMarkerStyle(0);
    gSyst_ph->SetMarkerColor(kGray+2);
    gSyst_ph->SetFillColorAlpha(kGray, 0.45); // alpha controls transparency
    gSyst_ph->SetLineColor(kGray);
    gSyst_ph->Draw("P E2 SAME");

    if (iP == 0) gStat->GetXaxis()->SetLimits(1.0, 20.0);
    else gStat->GetXaxis()->SetLimits(2.0, 20.0);

    TLegend *legend = new TLegend(0.2, 0.25, 0.5, 0.35);
    legend->SetFillColor(kOrange);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gStat, "sPHENIX 2024, |#eta| < 2.0");
    legend->AddEntry(gStat_ph, "PHENIX PRD 103 (2021) 052009, |#eta| < 0.35");
    legend->Draw();

    gPad->Modified();
    gPad->Update();
    double min_x = gPad->GetUxmin();
    double max_x = gPad->GetUxmax();

    // Add "sPHENIX internal"
    TPad *p = new TPad("p","p",0.,0.,1.,1.); p->SetFillStyle(0); p->Draw(); p->cd();
    TBox *whiteBox = new TBox(0.17, 0.72, 0.46, 0.90);
    whiteBox->Draw();
    canvas->cd();
    whiteBox->SetFillColorAlpha(kWhite, 1);
    std::stringstream stream;
    stream.str("");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.22, 0.85, "#font[72]{sPHENIX} Internal");
    latex.DrawLatex(0.22, 0.79, "p^{#uparrow}+p #sqrt{s} = 200 GeV");
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.19, 0.73, "7% polarization scale uncertainty not shown");
    latex.SetTextSize(0.07);

    if (iP == 0) {
      latex.DrawLatex(0.57, 0.84, "p^{#uparrow}+p #rightarrow #pi^{0} X");
    } else if (iP == 1) {
      latex.DrawLatex(0.57, 0.84, "p^{#uparrow}+p #rightarrow #eta X");
    }
    
    TLine *tline = new TLine();
    tline->SetLineWidth(2);
    tline->SetLineColor(kBlack);
    tline->SetLineStyle(kDashed);
    tline->DrawLine(min_x, 0, max_x, 0);

    canvas->Update();
    canvas->Draw();
    canvas->SaveAs((plots_folder_pt +  "/" + canvas_name.str() + ".png").c_str());
    canvas->SaveAs((plots_folder_pt +  "/" + canvas_name.str() + ".pdf").c_str());
  }
  gSystem->Exit(0);
}
