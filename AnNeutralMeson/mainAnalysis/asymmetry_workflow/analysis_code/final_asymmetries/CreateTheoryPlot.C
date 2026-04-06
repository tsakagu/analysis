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

void CreateTheoryPlot()
{
  std::string plots_folder_pt = "pt_asymmetries/THEORY/";
  std::string plots_folder_eta = "eta_asymmetries/THEORY/";
  std::string plots_folder_xf = "xf_asymmetries/THEORY/";

  gSystem->Exec(("mkdir -p " + plots_folder_pt).c_str());
  gSystem->Exec(("mkdir -p " + plots_folder_eta).c_str());
  gSystem->Exec(("mkdir -p " + plots_folder_xf).c_str());

  
  const std::string inputfolder = "asym_values/";
  const std::string inputfolder_theory_qgq = "sPHENIX26_AN/";
  const int nParticles = 2;
  const std::string particles[nParticles] = {"pi0", "eta"};
  const int nDirections = 2;
  const std::string directions[nDirections] = {"forward", "backward"};
  for (int iP = 0; iP < nParticles; iP++) {
    for (int iDir = 0; iDir < nDirections; iDir++) {
      std::stringstream inputfilename;
      inputfilename << inputfolder << "asym_summary_" << particles[iP] << "_" << directions[iDir] << "_pT" << ".txt";
      std::stringstream canvas_name;
      canvas_name << "canvas_asym_" << particles[iP] << "_" << directions[iDir] << "_pT";
      std::vector<double> x, y, estat, esyst, exsyst;
      parse_file(inputfilename.str(), x, y, estat, esyst, exsyst, 0.12);

      // Parse theory plots
      std::stringstream th_filename;
      th_filename << inputfolder_theory_qgq << "/" << particles[iP] << "_AN_vs_pT-" << directions[iDir] << "-200.dat";
      std::vector<double> x_th, y_th, ey_th;
      parse_theory_file(th_filename.str(), x_th, y_th, ey_th);

      const int n = x.size();
      auto gSyst = new TGraphErrors(n);
      auto gStat = new TGraphErrors(n);
      auto gTh_qgq = new TGraphErrors(n-1);
      int i_th = 0;
      for (int i = 0; i < n; ++i) {
        gSyst->SetPoint(i, x[i], y[i]);
        gSyst->SetPointError(i, exsyst[i], esyst[i]);
        
        gStat->SetPoint(i, x[i], y[i]);
        gStat->SetPointError(i, 0.0, estat[i]);

        if (iP == 0 && x[i] < 10 ||
            iP == 1 && x[i] > 2) {
          gTh_qgq->SetPoint(i_th, x_th[i_th], y_th[i_th]);
          gTh_qgq->SetPointError(i_th, 0, ey_th[i_th]);
          i_th++;
        }
      }
      gStat->Print();
      gSyst->Print();
      gTh_qgq->Print();
      SetsPhenixStyle();

      TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "", 1600, 900);
      gStat->SetTitle(";p_{T} GeV;A_{N}");

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
      gStat->SetMarkerStyle(kFullCircle);
      gStat->SetMarkerSize(2.2);
      gStat->SetMinimum(-2 * y_bound);
      gStat->SetMaximum(2 * y_bound);
      gStat->Draw("AP E1");
      
      // Style systematic boxes (transparent fill)
      if (iP == 0) {
        gSyst->SetFillColorAlpha(kPink-9, 0.35); // alpha controls transparency
        gSyst->SetLineColor(kPink-9);
      } else {
        gSyst->SetFillColorAlpha(kAzure-4, 0.35); // alpha controls transparency
        gSyst->SetLineColor(kAzure-4);
      }
      gSyst->SetMarkerStyle(0);
      gSyst->Draw("P E2 SAME");

      // Theory plots
      gTh_qgq->SetLineWidth(3);
      gTh_qgq->SetLineStyle(5);
      gTh_qgq->SetLineColor(kGreen);
      gTh_qgq->SetMarkerColor(kGreen);
      gTh_qgq->SetMarkerSize(0.0);
      gTh_qgq->Draw("SAME");
      
      if (iP == 0) gStat->GetXaxis()->SetLimits(1.0, 10.0);
      else gStat->GetXaxis()->SetLimits(2.0, 13.0);

      TLegend *legend = new TLegend(0.2, 0.25, 0.4, 0.35);
      legend->SetTextSize(0.04);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->AddEntry(gTh_qgq, "qgq contribution");
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
      latex.DrawLatex(0.22, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");
      latex.SetTextSize(0.03);
      latex.DrawLatex(0.19, 0.68, "7% polarization scale uncertainty not shown");
      latex.SetTextSize(0.05);

      if (iP == 0 && iDir == 0) {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #pi^{0} X,   x_{F}>0");
      } else if (iP == 1 && iDir == 0) {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #eta X,   x_{F}>0");
      } else if (iP == 0 && iDir == 1) {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #pi^{0} X,   x_{F}<0");
      } else if (iP == 1 && iDir == 1) {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #eta X,   x_{F}<0");
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
    {
      std::stringstream inputfilename;
      inputfilename << inputfolder << "asym_summary_" << particles[iP] << "_eta" << ".txt";
      std::stringstream canvas_name;
      canvas_name << "canvas_asym_" << particles[iP] << "_eta";
      std::vector<double> x, y, estat, esyst, exsyst;

      parse_file(inputfilename.str(), x, y, estat, esyst, exsyst, 0.05);

      // Parse theory plots
      std::stringstream th_filename;
      th_filename << inputfolder_theory_qgq << "/" << particles[iP] << "_AN_vs_eta-200.dat";
      std::vector<double> x_th, y_th, ey_th;
      parse_theory_file(th_filename.str(), x_th, y_th, ey_th);

      const int n = x.size();
      auto gSyst = new TGraphErrors(n);
      auto gStat = new TGraphErrors(n);
      auto gTh_qgq = new TGraphErrors(n);

      for (int i = 0; i < n; ++i) {
        gSyst->SetPoint(i, x[i], y[i]);
        gSyst->SetPointError(i, exsyst[i], esyst[i]);
        
        gStat->SetPoint(i, x[i], y[i]);
        gStat->SetPointError(i, 0.0, estat[i]);

        gTh_qgq->SetPoint(i, x_th[i], y_th[i]);
        gTh_qgq->SetPointError(i, 0, ey_th[i]);
      }
      gStat->Print();
      gSyst->Print();
      SetsPhenixStyle();

      TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "", 1600, 900);
      gStat->SetTitle(";#eta;A_{N}");

      double y_bound = 0;
      y_bound = GetAbsMaxInRange(gStat, -3, 3);

      int color = (iP == 0 ? kRed : kBlue);

      gStat->SetLineWidth(2);
      gStat->SetLineColor(color);
      gStat->SetMarkerColor(color);
      gStat->SetMarkerStyle(kFullCircle);
      gStat->SetMarkerSize(2.2);
      gStat->SetMinimum(-2 * y_bound);
      gStat->SetMaximum(2 * y_bound);
      gStat->Draw("AP E1");

      // Theory plots
      gTh_qgq->SetLineWidth(3);
      gTh_qgq->SetLineStyle(5);
      gTh_qgq->SetLineColor(kGreen);
      gTh_qgq->SetMarkerColor(kGreen);
      gTh_qgq->SetMarkerSize(0.0);
      gTh_qgq->Draw("SAME");

       // Style systematic boxes (transparent fill)
      if (iP == 0) {
        gSyst->SetFillColorAlpha(kPink-9, 0.35); // alpha controls transparency
        gSyst->SetLineColor(kPink-9);
      } else {
        gSyst->SetFillColorAlpha(kAzure-4, 0.35); // alpha controls transparency
        gSyst->SetLineColor(kAzure-4);
      }
      gSyst->SetMarkerStyle(0);
      gSyst->Draw("P E2 SAME");

      TLegend *legend = new TLegend(0.2, 0.25, 0.4, 0.35);
      legend->SetTextSize(0.04);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->AddEntry(gTh_qgq, "qgq contribution");
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
      latex.DrawLatex(0.22, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");
      latex.SetTextSize(0.03);
      latex.DrawLatex(0.19, 0.68, "7% polarization scale uncertainty not shown");
      latex.SetTextSize(0.05);

      if (iP == 0) {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #pi^{0} X");
      } else {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #eta X");
      }

      TLine *tline = new TLine();
      tline->SetLineWidth(2);
      tline->SetLineColor(kBlack);
      tline->SetLineStyle(kDashed);
      tline->DrawLine(min_x, 0, max_x, 0);

      canvas->Update();
      canvas->Draw();
      canvas->SaveAs((plots_folder_eta +  "/" + canvas_name.str() + ".png").c_str());
      canvas->SaveAs((plots_folder_eta +  "/" + canvas_name.str() + ".pdf").c_str());
    }
    {
      std::stringstream inputfilename;
      inputfilename << inputfolder << "asym_summary_" << particles[iP] << "_xf" << ".txt";
      std::stringstream canvas_name;
      canvas_name << "canvas_asym_" << particles[iP] << "_xf";
      std::vector<double> x, y, estat, esyst, exsyst;

      parse_file(inputfilename.str(), x, y, estat, esyst, exsyst, 0.002);

      // Parse theory plots
      std::stringstream th_filename;
      th_filename << inputfolder_theory_qgq << "/" << particles[iP] << "_AN_vs_xF-200.dat";
      std::vector<double> x_th, y_th, ey_th;
      parse_theory_file(th_filename.str(), x_th, y_th, ey_th);

      const int n = x.size();
      auto gSyst = new TGraphErrors(n);
      auto gStat = new TGraphErrors(n);
      auto gTh_qgq = new TGraphErrors(n);

      int ind_th = 0;
      for (int i = 0; i < n; ++i) {
        gSyst->SetPoint(i, x[i], y[i]);
        gSyst->SetPointError(i, exsyst[i], esyst[i]);
        
        gStat->SetPoint(i, x[i], y[i]);
        gStat->SetPointError(i, 0.0, estat[i]);

        gTh_qgq->SetPoint(i, x_th[i], y_th[i]);
        gTh_qgq->SetPointError(i, 0, ey_th[i]);
      }
      gStat->Print();
      gSyst->Print();
      SetsPhenixStyle();

      TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "", 1600, 900);
      gStat->SetTitle(";x_{F};A_{N}");

      double y_bound = 0;
      y_bound = GetAbsMaxInRange(gStat, -0.3, 0.3);

      int color = (iP == 0 ? kRed : kBlue);

      gStat->SetLineWidth(2);
      gStat->SetLineColor(color);
      gStat->SetMarkerColor(color);
      gStat->SetMarkerStyle(kFullCircle);
      gStat->SetMarkerSize(2.2);
      gStat->SetMinimum(-2 * y_bound);
      gStat->SetMaximum(2 * y_bound);
      gStat->Draw("AP E1");

       // Style systematic boxes (transparent fill)
      if (iP == 0) {
        gSyst->SetFillColorAlpha(kPink-9, 0.35); // alpha controls transparency
        gSyst->SetLineColor(kPink-9);
      } else {
        gSyst->SetFillColorAlpha(kAzure-4, 0.35); // alpha controls transparency
        gSyst->SetLineColor(kAzure-4);
      }
      gSyst->SetMarkerStyle(0);
      gSyst->Draw("P E2 SAME");

      // Theory plots
      gTh_qgq->SetLineWidth(3);
      gTh_qgq->SetLineStyle(5);
      gTh_qgq->SetLineColor(kGreen);
      gTh_qgq->SetMarkerColor(kGreen);
      gTh_qgq->SetMarkerSize(0.0);
      gTh_qgq->Draw("SAME");

      TLegend *legend = new TLegend(0.2, 0.25, 0.4, 0.35);
      legend->SetTextSize(0.04);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->AddEntry(gTh_qgq, "qgq contribution");
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
      latex.DrawLatex(0.22, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");
      latex.SetTextSize(0.03);
      latex.DrawLatex(0.19, 0.68, "7% polarization scale uncertainty not shown");
      latex.SetTextSize(0.05);

      if (iP == 0) {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #pi^{0} X");
      } else {
        latex.DrawLatex(0.5, 0.85, "p^{#uparrow}+p #rightarrow #eta X");
      }

      TLine *tline = new TLine();
      tline->SetLineWidth(2);
      tline->SetLineColor(kBlack);
      tline->SetLineStyle(kDashed);
      tline->DrawLine(min_x, 0, max_x, 0);

      canvas->Update();
      canvas->Draw();
      canvas->SaveAs((plots_folder_xf +  "/" + canvas_name.str() + ".png").c_str());
      canvas->SaveAs((plots_folder_xf +  "/" + canvas_name.str() + ".pdf").c_str());
    }
  }
  gSystem->Exit(0);
}
