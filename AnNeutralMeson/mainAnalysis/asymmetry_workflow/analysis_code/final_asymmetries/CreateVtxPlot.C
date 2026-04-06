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

void CreateVtxPlot()
{
  std::string plots_folder_xf = "xf_asymmetries_low_vtx/SYSTEMATIC/";

  gSystem->Exec(("mkdir -p " + plots_folder_xf).c_str());

  
  const std::string inputfolder = "asym_values/";
  const int nParticles = 2;
  const std::string particles[nParticles] = {"pi0", "eta"};
  const int nDirections = 2;
  const std::string directions[nDirections] = {"forward", "backward"};
  for (int iP = 0; iP < nParticles; iP++) {
    {
      std::stringstream inputfilename_low;
      inputfilename_low << inputfolder << "asym_summary_" << particles[iP] << "_xf_low_vtx" << ".txt";
      std::stringstream inputfilename_high;
      inputfilename_high << inputfolder << "asym_summary_" << particles[iP] << "_xf_high_vtx" << ".txt";
      std::stringstream canvas_name;
      canvas_name << "canvas_asym_" << particles[iP] << "_xf_vs_vertex";
      std::vector<double> x_low, y_low, estat_low, esyst_low, exsyst_low;
      std::vector<double> x_high, y_high, estat_high, esyst_high, exsyst_high;

      parse_file(inputfilename_low.str(), x_low, y_low, estat_low, esyst_low, exsyst_low, 0.002);
      parse_file(inputfilename_high.str(), x_high, y_high, estat_high, esyst_high, exsyst_high, 0.002);

      const int n = x_low.size();
      auto gSyst_low = new TGraphErrors(n);
      auto gStat_low = new TGraphErrors(n);
      auto gSyst_high = new TGraphErrors(n);
      auto gStat_high = new TGraphErrors(n);

      for (int i = 0; i < n; ++i) {
        gSyst_low->SetPoint(i, x_low[i], y_low[i]);
        gSyst_low->SetPointError(i, exsyst_low[i], esyst_low[i]);
        
        gStat_low->SetPoint(i, x_low[i], y_low[i]);
        gStat_low->SetPointError(i, 0.0, estat_low[i]);

        gSyst_high->SetPoint(i, x_high[i], y_high[i]);
        gSyst_high->SetPointError(i, exsyst_high[i], esyst_high[i]);
        
        gStat_high->SetPoint(i, x_high[i], y_high[i]);
        gStat_high->SetPointError(i, 0.0, estat_high[i]);
      }
      SetsPhenixStyle();

      TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "", 1600, 900);
      gStat_low->SetTitle(";x_{F};A_{N}");

      double y_bound = 0;
      y_bound = GetAbsMaxInRange(gStat_low, -0.3, 0.3);

      int color = (iP == 0 ? kRed : kBlue);

      gStat_low->SetLineWidth(2);
      gStat_low->SetLineColor(color);
      gStat_low->SetMarkerColor(color);
      gStat_low->SetMarkerStyle(kFullCircle);
      gStat_low->SetMarkerSize(2.2);
      gStat_low->SetMinimum(-2 * y_bound);
      gStat_low->SetMaximum(2 * y_bound);
      gStat_low->Draw("AP E1");

      gStat_high->SetLineWidth(2);
      gStat_high->SetLineColor(kViolet);
      gStat_high->SetMarkerColor(kViolet);
      gStat_high->SetMarkerStyle(kFullCircle);
      gStat_high->SetMarkerSize(2.2);
      gStat_high->SetMinimum(-2 * y_bound);
      gStat_high->SetMaximum(2 * y_bound);
      gStat_high->Draw("P E1 SAME");

       // Style systematic boxes (transparent fill)
      if (iP == 0) {
        gSyst_low->SetFillColorAlpha(kPink-9, 0.35); // alpha controls transparency
        gSyst_low->SetLineColor(kPink-9);
        gSyst_high->SetFillColorAlpha(kViolet-9, 0.35); // alpha controls transparency
        gSyst_high->SetLineColor(kViolet-9);
      } else {
        gSyst_low->SetFillColorAlpha(kAzure-4, 0.35); // alpha controls transparency
        gSyst_low->SetLineColor(kAzure-4);
      }
      gSyst_low->SetMarkerStyle(0);
      gSyst_low->Draw("P E2 SAME");
      gSyst_high->SetMarkerStyle(0);
      gSyst_high->Draw("P E2 SAME");

      TLegend *legend = new TLegend(0.25, 0.25, 0.5, 0.4);
      legend->SetTextSize(0.05);
      legend->AddEntry(gStat_low, "|z_{vtx}| < 30 cm");
      legend->AddEntry(gStat_high, "|z_{vtx}| > 30 cm");
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
