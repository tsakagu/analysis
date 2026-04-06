#include <any>
#include <cmath>

using namespace std;

bool saveFinalPlots = true;
bool printMarkdownTables = true;
bool printLatexTables = true;

//Custom binning scheme for LF analysis
float lower_bin_bounds[] = {0.5, 0.8, 1.1, 1.4, 1.8, 2.2, 3, 4};
const unsigned int n_variable_bins = sizeof(lower_bin_bounds)/sizeof(lower_bin_bounds[0]) - 1; 

template <typename T>
string to_string_with_precision(const T a_value, const int n = 0)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

//Initialise histograms
TH1F makeHisto(int nBins, float min, float max, string type, string xAxisTitle, string unit, int precision, string yAxisTitle = "Geo. Accept.")
{
  string histo_name = type + "_" + xAxisTitle;
  TH1F myHisto(histo_name.c_str(), histo_name.c_str(), nBins, min, max);
  
  if (unit != "") xAxisTitle += " [" + unit +  "]";  
  myHisto.GetXaxis()->SetTitle(xAxisTitle.c_str());

  float binWidth = (float) (max - min)/nBins;
  if (unit != "") yAxisTitle += " / " + to_string_with_precision(binWidth, precision) + " " + unit;
  myHisto.GetYaxis()->SetTitle(yAxisTitle.c_str());
  
  return myHisto;
}

TH1F makeHisto(int nBins, float* min, string type, string xAxisTitle, string unit, int precision, string yAxisTitle = "Geo. Accept.")
{
  string histo_name = type + "_" + xAxisTitle;
  TH1F myHisto(histo_name.c_str(), histo_name.c_str(), nBins, min);

  if (unit != "") xAxisTitle += " [" + unit +  "]";
  myHisto.GetXaxis()->SetTitle(xAxisTitle.c_str());
  myHisto.GetYaxis()->SetTitle(yAxisTitle.c_str());

  return myHisto;
}

//What branches to read into the analysis
struct variableMap
{
  float floats[7];
  int ints[2];

  //any and variant are a pain here with TBranch 
  map<string, float> float_vars{
    {"true_mother_mass", floats[0]},
    {"true_mother_p", floats[1]},
    {"true_mother_pT", floats[2]},
    {"true_mother_eta", floats[3]},
    {"true_mother_phi", floats[4]},
    {"true_primary_vertex_z", floats[5]},
    {"true_secondary_vertex_z", floats[6]},
  };
    
  map<string, int> int_vars{
    {"reco_track_1_silicon_seeds", ints[0]},
    {"reco_track_2_silicon_seeds", ints[1]},
  };
};

//What histograms to build for each ratio check
class histos
{
  public:
    TH1F Lambda0_all;
    TH1F K_S0_all;
    TH1F Lambda0_inGeo;
    TH1F K_S0_inGeo;
    TH1F ratio;
    TH1F inv_ratio;
    bool variable_bins = false;
    int nBins = 15;

  histos(bool constructor_variable_bins = false, float range = 1., string variable = "", string unit = "")
  {
    variable_bins = constructor_variable_bins;
    if (constructor_variable_bins)
    {
      Lambda0_all   = makeHisto(n_variable_bins, lower_bin_bounds, "Lambda0_all", variable.c_str(), unit.c_str(), 1);
      K_S0_all      = makeHisto(n_variable_bins, lower_bin_bounds, "K_S0_all", variable.c_str(), unit.c_str(), 1);
      Lambda0_inGeo = makeHisto(n_variable_bins, lower_bin_bounds, "Lambda0_inGeo", variable.c_str(), unit.c_str(), 1);
      K_S0_inGeo    = makeHisto(n_variable_bins, lower_bin_bounds, "K_S0_inGeo", variable.c_str(), unit.c_str(), 1);
      ratio         = makeHisto(n_variable_bins, lower_bin_bounds, "ratio", variable.c_str(), unit.c_str(), 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
      inv_ratio     = makeHisto(n_variable_bins, lower_bin_bounds, "inv_ratio", variable.c_str(), unit.c_str(), 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
    }
    else
    {
      Lambda0_all   = makeHisto(nBins, -1*range, range, "Lambda0_all", variable.c_str(), unit.c_str(), 1);
      K_S0_all      = makeHisto(nBins, -1*range, range, "K_S0_all", variable.c_str(), unit.c_str(), 1);
      Lambda0_inGeo = makeHisto(nBins, -1*range, range, "Lambda0_inGeo", variable.c_str(), unit.c_str(), 1);
      K_S0_inGeo    = makeHisto(nBins, -1*range, range, "K_S0_inGeo", variable.c_str(), unit.c_str(), 1);
      ratio         = makeHisto(nBins, -1*range, range, "ratio", variable.c_str(), unit.c_str(), 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
      inv_ratio     = makeHisto(nBins, -1*range, range, "inv_ratio", variable.c_str(), unit.c_str(), 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
    }
  }
};

histos pT(true, 0, "pT", "GeV");
histos eta(false, 1.1, "#eta", "");
histos phi(false, M_PI, "#phi", "");
histos rap(false, 1.0, "y", "");

template <typename T>
void savePlots(T myPlot, string plotName, bool logY = false, float yMin = 0, float yMax = 1)
{
  TGaxis::SetMaxDigits(3);
  string plotPath = "plots/";
  string makeDirectory = "mkdir -p " + plotPath;
  system(makeDirectory.c_str());

  TCanvas *c1  = new TCanvas("myCanvas", "myCanvas",800,800);

  myPlot.GetYaxis()->SetRangeUser(yMin, yMax);

  if (strncmp(typeid(myPlot).name(), "4TH2F", 5) == 0)
  {
    myPlot.Draw("COLZ");
  }
  else
  {
    myPlot.Sumw2();
    if (logY) gPad->SetLogy();
    myPlot.Draw("PE1");
  }

  TPaveText *pt;
  pt = new TPaveText(0.15,0.9,0.95,1., "NDC");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  TText *pt_LaTex = pt->AddText("#it{#bf{sPHENIX}} Simulation, #sqrt{s} = 200 GeV, pp");
  pt->SetBorderSize(0);
  pt->Draw();
  gPad->Modified();

  string extensions[] = {".C", ".pdf", ".png", ".root"};
  for (auto extension : extensions)
  {
    string output = plotPath + "/" + plotName + extension;
    c1->SaveAs(output.c_str());
  }
}

void processData(string type = "Kshort2pipi")
{
  string fileName = type == "Kshort2pipi" ? "files/outputHFTrackEff_Kshort2pipi_CR_2_mode_pTref_3p0_extra_high_pT_points.root"
                                          : "files/outputHFTrackEff_Lambda2ppi_CR_2_mode_pTref_3p0_extra_high_pT_points.root";
  TFile *file = new TFile(fileName.c_str());
  TTree* data = (TTree*)file->Get("HFTrackEfficiency");

  variableMap inputMap;
  for (auto &[branch, var] : inputMap.float_vars) data->SetBranchAddress(branch.c_str(), &var);
  for (auto &[branch, var] : inputMap.int_vars)   data->SetBranchAddress(branch.c_str(), &var);

  int tmp = 0;
  int barWidth = 50;

  //Geometric Acceptance
  int num_entries = data->GetEntries();
  for (int  l = 0; l < num_entries; ++l)
  {
    if (tmp != (int)100*l/num_entries)
    {
      tmp = (int)100*l/num_entries;
      if ((tmp%1)  == 0)
      {
        cout << "[";
        int pos = barWidth * tmp/100;
        for (int i = 0; i < barWidth; ++i)
        {
          if (i < pos) cout << "=";
          else if (i == pos) cout << ">";
          else cout << " ";
        }
        cout << "] " << tmp << " %\r";
        cout.flush();
      }
    }

    data->GetEntry(l);

    int sign = inputMap.float_vars["true_secondary_vertex_z"] > inputMap.float_vars["true_primary_vertex_z"] ? 1 : -1;
    float pz = sqrt(pow(inputMap.float_vars["true_mother_p"], 2) - pow(inputMap.float_vars["true_mother_pT"], 2));
    float E = sqrt(pow(inputMap.float_vars["true_mother_p"], 2) + pow(inputMap.float_vars["true_mother_mass"], 2));
    float rapidity = sign*0.5*log((E + pz)/(E - pz));

    if (abs(rapidity) > 1.) continue;
 
    if (type == "Kshort2pipi")
    {
      pT.K_S0_all.Fill(inputMap.float_vars["true_mother_pT"]);
      eta.K_S0_all.Fill(inputMap.float_vars["true_mother_eta"]);
      rap.K_S0_all.Fill(rapidity);
      phi.K_S0_all.Fill(inputMap.float_vars["true_mother_phi"]);
    }
    else
    {
      pT.Lambda0_all.Fill(inputMap.float_vars["true_mother_pT"]);
      eta.Lambda0_all.Fill(inputMap.float_vars["true_mother_eta"]);
      rap.Lambda0_all.Fill(rapidity);
      phi.Lambda0_all.Fill(inputMap.float_vars["true_mother_phi"]);
    }

    bool accepted = min(inputMap.int_vars["reco_track_1_silicon_seeds"], inputMap.int_vars["reco_track_2_silicon_seeds"]) > 0;

    if (accepted)
    {
      if (type == "Kshort2pipi")
      {
        pT.K_S0_inGeo.Fill(inputMap.float_vars["true_mother_pT"]);
        eta.K_S0_inGeo.Fill(inputMap.float_vars["true_mother_eta"]);
        rap.K_S0_inGeo.Fill(rapidity);
        phi.K_S0_inGeo.Fill(inputMap.float_vars["true_mother_phi"]);
      }
      else
      {
        pT.Lambda0_inGeo.Fill(inputMap.float_vars["true_mother_pT"]);
        eta.Lambda0_inGeo.Fill(inputMap.float_vars["true_mother_eta"]);
        rap.Lambda0_inGeo.Fill(rapidity);
        phi.Lambda0_inGeo.Fill(inputMap.float_vars["true_mother_phi"]);
      }
    }
  }

  cout << "[";
  int pos = barWidth * tmp;
  for (int i = 0; i < barWidth; ++i)
  {
    if (i < pos) cout << "=";
    else if (i == pos) cout << ">";
    else cout << " ";
  }
  cout << "] 100 %\r";
  cout.flush();
  cout<<endl;
}

void makeRatiosAndSave(histos &histoSet, string trailer)
{
  string K_S0_saveName = "KS0_geometric_acceptance_ratio_" + trailer;
  string Lambda0_saveName = "Lambda0_geometric_acceptance_ratio_" + trailer;
  string ratio_saveName = "Lambda0_to_KS0_geometric_acceptance_ratio_" + trailer;
  string inv_ratio_saveName = "K_S0_to_Lambda0_geometric_acceptance_ratio_" + trailer;

  histoSet.Lambda0_all.Sumw2();
  histoSet.K_S0_all.Sumw2();

  histoSet.Lambda0_inGeo.Sumw2();
  histoSet.K_S0_inGeo.Sumw2();

  histoSet.Lambda0_inGeo.Divide(&histoSet.Lambda0_all);
  histoSet.K_S0_inGeo.Divide(&histoSet.K_S0_all);
  savePlots(histoSet.K_S0_inGeo, K_S0_saveName.c_str(), false, 0, 0.035);
  savePlots(histoSet.Lambda0_inGeo, Lambda0_saveName.c_str(), false, 0, 0.035);

  histoSet.ratio = histoSet.Lambda0_inGeo;
  histoSet.ratio.Divide(&histoSet.K_S0_inGeo);

  histoSet.inv_ratio = histoSet.K_S0_inGeo;
  histoSet.inv_ratio.Divide(&histoSet.Lambda0_inGeo);

  savePlots(histoSet.ratio, ratio_saveName.c_str(), false, 0, 1);
  savePlots(histoSet.inv_ratio, inv_ratio_saveName.c_str(), false, 0, 3);
}

void makeTable(string var, histos theseHistos, string type = "Markdown")
{
  cout << endl;
 
  string startline = type == "Markdown" ? "| " : "";
  string separator = type == "Markdown" ? " | " : " & ";
  string endline  = type == "Markdown" ? " | " : " \\\\ ";

  int nBins = theseHistos.variable_bins ? n_variable_bins : theseHistos.nBins;

  if (type != "Markdown") cout << "\\begin{tabular}{@{}ccccc@{}}" << endl << "\\toprule[1pt]" << endl;
  
  cout << startline << var << separator << "$K_S^0$" << separator << "$\\Lambda^0$" << separator << "$\\Lambda^0 / K_S^0$" << separator << "$K_S^0 / \\Lambda^0$" << endline << endl;

  if (type == "Markdown") cout << "|:--:|:--:|:--:|:--:|:--:|" << endl;
  else cout << "\\midrule[0.2pt]" << endl;

  for (int i = 1; i <= nBins; ++i)
  {
    string low = to_string_with_precision(theseHistos.ratio.GetXaxis()->GetBinLowEdge(i), 2);
    string high = to_string_with_precision(theseHistos.ratio.GetXaxis()->GetBinUpEdge(i), 2);

    string Ks0_content = to_string_with_precision(theseHistos.K_S0_inGeo.GetBinContent(i), 4);
    string Ks0_error = to_string_with_precision(theseHistos.K_S0_inGeo.GetBinError(i), 4);

    string Lambda_content = to_string_with_precision(theseHistos.Lambda0_inGeo.GetBinContent(i), 4);
    string Lambda_error = to_string_with_precision(theseHistos.Lambda0_inGeo.GetBinError(i), 4);

    string content = to_string_with_precision(theseHistos.ratio.GetBinContent(i), 2);
    string error = to_string_with_precision(theseHistos.ratio.GetBinError(i), 2);

    string inv_content = to_string_with_precision(theseHistos.inv_ratio.GetBinContent(i), 2);
    string inv_error = to_string_with_precision(theseHistos.inv_ratio.GetBinError(i), 2);
    cout << startline << low << " $\\rightarrow$ " << high << separator << Ks0_content << " $\\pm$ " << Ks0_error << separator << Lambda_content << " $\\pm$ " << Lambda_error << separator << content << " $\\pm$ " << error << separator << inv_content << " $\\pm$ " << inv_error << endline << endl;
  }

  if (type != "Markdown") cout << "\\bottomrule[1pt]" << endl << "\\end{tabular}" << endl;
}

void calcGeomAccept()
{
  cout << "Processing K-short data" << endl;
  processData();
  cout << "Processing Lambda0 data" << endl;
  processData("Lambda02ppi");

  if (saveFinalPlots)
  {
    makeRatiosAndSave(pT, "pT");
    makeRatiosAndSave(eta, "eta");
    makeRatiosAndSave(phi, "phi");
    makeRatiosAndSave(rap, "rap");
  }

  if (printMarkdownTables)
  {
    makeTable("$p_{T}$ [GeV]", pT, "Markdown");
    makeTable("$\\eta$", eta, "Markdown");
    makeTable("$y$", rap, "Markdown");
    makeTable("$\\phi$", phi, "Markdown");
  }

  if (printLatexTables)
  {
    makeTable("$p_{T}$ [GeV]", pT, "LaTex");
    makeTable("$\\eta$", eta, "LaTex");
    makeTable("$y$", rap, "LaTex");
    makeTable("$\\phi$", phi, "LaTex");
  }
}
