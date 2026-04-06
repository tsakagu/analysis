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

void CreateForwardPlot()
{
  std::string plots_folder_pt = "pt_asymmetries/SYSTEMATIC/";

  gSystem->Exec(("mkdir -p " + plots_folder_pt).c_str());
  
  const std::string inputfolder = "asym_values/";

  float pTBounds[2] = {-0.010, 0.020};
  
  std::stringstream inputfilename_low;
  inputfilename_low << inputfolder << "asym_summary_pi0_pT_low_xf" << ".txt";
  std::stringstream inputfilename_high;
  inputfilename_high << inputfolder << "asym_summary_pi0_pT_high_xf" << ".txt";
  std::stringstream canvas_name;
  canvas_name << "canvas_asym_pi0_pT_vs_xF";
  std::vector<double> x_low, y_low, estat_low, esyst_low, exsyst_low;
  std::vector<double> x_high, y_high, estat_high, esyst_high, exsyst_high;

  parse_file(inputfilename_low.str(), x_low, y_low, estat_low, esyst_low, exsyst_low, 0.12);
  parse_file(inputfilename_high.str(), x_high, y_high, estat_high, esyst_high, exsyst_high, 0.12);
  
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
  }

  // Skip first bin (too little stat for large xF)
  for (int i = 1; i < n; ++i) {
    gSyst_high->SetPoint(i-1, x_high[i]+0.1, y_high[i]);
    gSyst_high->SetPointError(i-1, exsyst_high[i], esyst_high[i]);

    gStat_high->SetPoint(i-1, x_high[i]+0.1, y_high[i]);
    gStat_high->SetPointError(i-1, 0.0, estat_high[i]);
  }
  
  SetsPhenixStyle();

  TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "", 1600, 900);
  gStat_low->SetTitle(";p_{T} [GeV];A_{N}");

  double y_bound = 0;
  y_bound = GetAbsMaxInRange(gStat_high, 1.5, 10);
  
  gStat_low->SetLineWidth(2);
  gStat_low->SetLineColor(kGray+2);
  gStat_low->SetMarkerColor(kGray+2);
  gStat_low->SetMarkerStyle(59);
  gStat_low->SetMarkerSize(2.6);
  // gStat_low->SetMinimum(-2 * y_bound);
  // gStat_low->SetMaximum(2 * y_bound);
  gStat_low->SetMinimum(pTBounds[0]);
  gStat_low->SetMaximum(pTBounds[1]);
  gStat_low->Draw("AP E1");

  gStat_high->SetLineWidth(2);
  gStat_high->SetLineColor(kRed);
  gStat_high->SetMarkerColor(kRed);
  gStat_high->SetMarkerStyle(kFullSquare);
  gStat_high->SetMarkerSize(2.4);
  // gStat_high->SetMinimum(-2 * y_bound);
  // gStat_high->SetMaximum(2 * y_bound);
  gStat_high->SetMinimum(pTBounds[0]);
  gStat_high->SetMaximum(pTBounds[1]);
  gStat_high->Draw("P E1 SAME");

  // Style systematic boxes (transparent fill)
  gSyst_low->SetFillColorAlpha(kGray, 0.65);
  gSyst_low->SetLineColor(kGray);
  gSyst_high->SetFillColorAlpha(kPink-9, 0.45);
  gSyst_high->SetLineColor(kPink-9);

  gSyst_low->SetMarkerStyle(0);
  gSyst_low->Draw("P E2 SAME");
  gSyst_high->SetMarkerStyle(0);
  gSyst_high->Draw("P E2 SAME");

  gStat_low->GetXaxis()->SetLimits(1.0, 10.0);
  
  TLegend *legend = new TLegend(0.22, 0.54, 0.38, 0.68);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(gStat_high, "x_{F} > 0.035");
  legend->AddEntry(gStat_low, "x_{F} < 0.035");
  legend->Draw();

  // TLegend *legend = new TLegend(0.62, 0.75, 0.88, 0.88);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);

  gPad->Modified();
  gPad->Update();
  double min_x = gPad->GetUxmin();
  double max_x = gPad->GetUxmax();

  // Add "sPHENIX internal"
  TPad *p = new TPad("p","p",0.,0.,1.,1.); p->SetFillStyle(0); p->Draw(); p->cd();
  TBox *whiteBox = new TBox(0.185, 0.72, 0.49, 0.90);
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
  latex.DrawLatex(0.57, 0.84, "p^{#uparrow}+p #rightarrow #pi^{0} X");

  TLine *tline = new TLine();
  tline->SetLineWidth(2);
  tline->SetLineColor(kGray+2);
  tline->SetLineStyle(kDashed);
  tline->DrawLine(min_x, 0, max_x, 0);

  canvas->Update();
  canvas->Draw();
  canvas->SaveAs((plots_folder_pt +  "/" + canvas_name.str() + ".png").c_str());
  canvas->SaveAs((plots_folder_pt +  "/" + canvas_name.str() + ".pdf").c_str());
  canvas->SaveAs((plots_folder_pt +  "/" + canvas_name.str() + ".C").c_str());
  gSystem->Exit(0);
}
