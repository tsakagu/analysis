#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

double GetAbsMaxInRangeSimple(TGraph *graph, double xMin, double xMax)
{
  if (!graph ) return -1; // Error handling

  double maxVal = -1;
  int n = graph->GetN();

  for (int i = 0; i < n; i++)
  {
    double x, y, yVal;
    graph->GetPoint(i, x, y);
    yVal = std::abs(y);
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

void trim_graph(TGraphErrors* graph_original,
                TGraphErrors* (&graph_trimmed),
                const float y_limit)
{
  graph_trimmed = new TGraphErrors();
  std::stringstream graphname; graphname << graph_original->GetName() << "_trimmed";
  graph_trimmed->SetName(graphname.str().c_str());
  graph_trimmed->SetTitle(graph_original->GetTitle());
  int ind = 0;
  for (int i = 0; i < graph_original->GetN(); i++) {
    double y = graph_original->GetPointY(i);
    double ey = graph_original->GetErrorY(i); 
    double yVal = std::max(std::abs(y-ey), std::abs(y+ey));
    if (yVal < y_limit) {
      graph_trimmed->SetPoint(ind, graph_original->GetPointX(i), graph_original->GetPointY(i));
      graph_trimmed->SetPointError(ind, 0, graph_original->GetErrorY(i));
      ind++;
    }
  }
}

void cut_graph(TGraphErrors* graph_original,
               TGraphErrors* (&graph_cut),
               const int iMin,
               const int iMax)
{
  graph_cut = new TGraphErrors();
  std::stringstream graphname; graphname << graph_original->GetName() << "_cut";
  graph_cut->SetName(graphname.str().c_str());
  graph_cut->SetTitle(graph_original->GetTitle());
  int ind = 0;
  for (int i = std::max(0, iMin) ; i <= std::min(iMax, graph_original->GetN()-1); i++) {
    double y = graph_original->GetPointY(i);
    double ey = graph_original->GetErrorY(i); 
    double yVal = std::max(std::abs(y-ey), std::abs(y+ey));
    graph_cut->SetPoint(ind, graph_original->GetPointX(i), graph_original->GetPointY(i));
    graph_cut->SetPointError(ind, 0, graph_original->GetErrorY(i));
    ind++;
  }
}

void shift_graph(TGraphErrors *graph, const float shift)
{
  const int nPoints = graph->GetN();
  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    double x = graph->GetPointX(iPoint);
    graph->SetPointX(iPoint, x + shift);
  }
}

void weighted_average(double &average_asym,
                      double &average_asym_err,
                      const double value_1,
                      const double value_2,
                      const double weight_1,
                      const double weight_2)
{
  average_asym = (value_1 * weight_1 + value_2 * weight_2) / (weight_1 + weight_2);
  average_asym_err = std::sqrt(1. / (weight_1 + weight_2));
}

void compute_t_statistic(TGraph* (&graph_t_test),
                         const std::string &graph_t_test_name,
                         TGraphErrors* graph_1,
                         TGraphErrors* graph_2
                         )
{
  graph_t_test = new TGraph();
  graph_t_test->SetName(graph_t_test_name.c_str());
  int nPoints = graph_1->GetN();
  if (graph_2->GetN() != nPoints) {
    std::cout << "Error, Averaging two histograms with unequal number of bins" << std::endl;
    gSystem->Exit(0);
  }

  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    double x1 = graph_1->GetPointX(iPoint);
    double x2 = graph_2->GetPointX(iPoint);
    double val1 = graph_1->GetPointY(iPoint);
    double val2 = graph_2->GetPointY(iPoint);
    double e1 = graph_1->GetErrorY(iPoint); double w1 = 1. / e1 / e1;
    double e2 = graph_2->GetErrorY(iPoint); double w2 = 1. / e2 / e2;
    double xMean, xErr; weighted_average(xMean, xErr, x1, x2, w1, w2);
    double t_stat = (val1 - val2) / std::sqrt(e1 * e1 + e2 * e2);
    graph_t_test->SetPoint(iPoint, xMean, t_stat);
  }
}

void extract_average_asymmetry(TGraphErrors* (&graph_average),
                               const std::string &graph_average_name,
                               const TGraphErrors *graph_1,
                               const TGraphErrors *graph_2)
{
  graph_average = new TGraphErrors();
  graph_average->SetName(graph_average_name.c_str());

  if (!graph_1 || !graph_2) {
    std::cout << "Error: Missing Graph" << std::endl;
    gSystem->Exit(0);
  }
  
  const int nPoints = graph_1->GetN();
  if (graph_2->GetN() != nPoints) {
    std::cout << "Error, Averaging two histograms with unequal number of bins" << std::endl;
    gSystem->Exit(0);
  }

  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    double x1 = graph_1->GetPointX(iPoint);
    double x2 = graph_2->GetPointX(iPoint);
    double val1 = graph_1->GetPointY(iPoint);
    double val2 = graph_2->GetPointY(iPoint);
    double e1 = graph_1->GetErrorY(iPoint); double w1 = 1. / e1 / e1;
    double e2 = graph_2->GetErrorY(iPoint); double w2 = 1. / e2 / e2;
    double xMean, xErr; weighted_average(xMean, xErr, x1, x2, w1, w2);
    double average_asym, average_asym_err; weighted_average(average_asym, average_asym_err, val1, val2, w1, w2);
    //std::cout << iPoint << ": " << x1 << ", " << x2 << ", " << val1 << ", " << val2 << ", " << e1 << ", " << e2 <<  ", " << xMean << ", " << average_asym << ", " << average_asym_err << std::endl;

    graph_average->SetPointX(iPoint, xMean);
    graph_average->SetPointY(iPoint, average_asym);
    graph_average->SetPointError(iPoint, 0, average_asym_err);
  }
}

void final_average_asym_eta_ana509(bool gpr_ratio = false,
                                   const std::string& special_tag = "")
{

  SetsPhenixStyle();
  std::string plots_folder = "eta_asymmetries" + special_tag + "/FINAL/";
  gSystem->Exec(("mkdir -p " + plots_folder).c_str());
  
  const std::string inputfolder = "/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/analysis_asymmetry/analysis_ana509";
  
  const int nBeams = 2;
  const std::string beams[nBeams] = {"blue", "yellow"};
  const int nParticles = 2;
  const std::string particles[nParticles] = {"pi0", "eta"};
  const int nMethods = 2;
  const std::string methods[nMethods] = {"sqrt", "rellum"};
  const int nTriggers = 2;
  const std::string trigger_labels[nTriggers] = {"MBD_pt_l", "efficiency_30_pt_g"};
  const std::string trigger_shorts[nTriggers] = {"MBD", "efficiency_30"};
  const int nPtThresholds = nParticles;
  const std::string pt_thresholds[nPtThresholds] = {"3", "4"};
  const int nCrossingAngles = 2;
  const std::string crossing_angles[nCrossingAngles] = {"0mrad", "15mrad"};

  // Store final graph in output file
  std::string results_folder = "final_results";
  std::string ratio_suffix = (gpr_ratio ? "_gpr" : "");
  gSystem->Exec(("mkdir -p " + results_folder).c_str());
  TFile *outputfile = new TFile(("full_average_eta" + special_tag + ratio_suffix + ".root").c_str(), "RECREATE"); 

  TFile *inputfile[nTriggers][nCrossingAngles][nParticles] = {nullptr};
  TGraphErrors *graph_AN[nTriggers][nCrossingAngles][nBeams][nParticles][nMethods] = {nullptr};


  // Load graphs
  for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
    for (int iCrossingAngle = 0; iCrossingAngle < nCrossingAngles; iCrossingAngle++) {
      for (int iParticle = 0; iParticle < nParticles; iParticle++) {
        std::string inputfilename = inputfolder + "/analysis_asymmetry_eta_complete_ana509_01312026" + special_tag + "_" + trigger_labels[iTrigger] + pt_thresholds[iParticle] + "_" + crossing_angles[iCrossingAngle] + "_averaged_fitted" + ratio_suffix + ".root";
        inputfile[iTrigger][iCrossingAngle][iParticle] = TFile::Open(inputfilename.c_str());
        if (!inputfile[iTrigger][iCrossingAngle][iParticle]) {
          std::cerr << "Error: Cannot open file " << inputfilename << std::endl;
          gSystem->Exit(0);
        }
        for (int iBeam = 0; iBeam < nBeams; iBeam++) {
          for (int iMethod = 0; iMethod < nMethods; iMethod++) {
            std::stringstream graphname;
            graphname << "graph_AN_peak_" << beams[iBeam] << "_" << particles[iParticle] << "_" << methods[iMethod];
            graph_AN[iTrigger][iCrossingAngle][iBeam][iParticle][iMethod] = (TGraphErrors*) inputfile[iTrigger][iCrossingAngle][iParticle]->Get(graphname.str().c_str());
            if (!graph_AN[iTrigger][iCrossingAngle][iBeam][iParticle][iMethod]) {
              std::cout << "Error. Can't open " << graphname.str() << " in file " << inputfilename << std::endl;
              gSystem->Exit(0);
            }
          }
        }
      }
    }
  }

  outputfile->cd();
  // New averaged graphs
  TGraphErrors *graph_AN_triggerA[nCrossingAngles][nBeams][nParticles][nMethods] = {nullptr}; // trigger-averaged
  TGraphErrors *graph_AN_angleA[nTriggers][nBeams][nParticles][nMethods] = {nullptr}; // angle-averaged
  TGraphErrors *graph_AN_triggerA_beamA[nCrossingAngles][nParticles][nMethods] = {nullptr}; // trigger-averaged and beam-averaged
  TGraphErrors *graph_AN_triggerA_angleA[nBeams][nParticles][nMethods] = {nullptr}; // trigger-averaged and angle-averaged
  TGraphErrors *graph_AN_angleA_beamA[nTriggers][nParticles][nMethods] = {nullptr}; // angle-averaged and beam-averaged
  TGraphErrors *graph_AN_final[nParticles][nMethods] = {nullptr}; // trigger-, beam- and angle-averaged

  // Average -> Separate Crossing Angles
  // Average over triggers
  for (int iCrossingAngle = 0; iCrossingAngle < nCrossingAngles; iCrossingAngle++) {
    for (int iBeam = 0; iBeam < nBeams; iBeam++) {
      for (int iMethod = 0; iMethod < nMethods; iMethod++) {
        for (int iParticle = 0; iParticle < nParticles; iParticle++) {
          std::stringstream graphname;
          graphname << "graph_AN_peak_" << beams[iBeam] << "_" << particles[iParticle] << "_" << methods[iMethod] << "_" << crossing_angles[iCrossingAngle];
          extract_average_asymmetry(
            graph_AN_triggerA[iCrossingAngle][iBeam][iParticle][iMethod],
            graphname.str(),
            graph_AN[0][iCrossingAngle][iBeam][iParticle][iMethod],
            graph_AN[1][iCrossingAngle][iBeam][iParticle][iMethod]);
          graph_AN_triggerA[iCrossingAngle][iBeam][iParticle][iMethod]->Write();
        }
      }
    }
  }
  // Then average over beam 
  for (int iCrossingAngle = 0; iCrossingAngle < nCrossingAngles; iCrossingAngle++) {
    for (int iMethod = 0; iMethod < nMethods; iMethod++) {
      for (int iParticle = 0; iParticle < nParticles; iParticle++) {
        std::stringstream graphname;
        graphname << "graph_AN_peak_" << particles[iParticle] << "_" << methods[iMethod] << "_" << crossing_angles[iCrossingAngle];
        extract_average_asymmetry(
          graph_AN_triggerA_beamA[iCrossingAngle][iParticle][iMethod],
          graphname.str(),
          graph_AN_triggerA[iCrossingAngle][0][iParticle][iMethod],
          graph_AN_triggerA[iCrossingAngle][1][iParticle][iMethod]);
        graph_AN_triggerA_beamA[iCrossingAngle][iParticle][iMethod]->Write();
      }
    }
  }

  // Average -> Separate Triggers
  // Average over crossing angle
  for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
    for (int iBeam = 0; iBeam < nBeams; iBeam++) {
      for (int iMethod = 0; iMethod < nMethods; iMethod++) {
        for (int iParticle = 0; iParticle < nParticles; iParticle++) {
          std::stringstream graphname;
          graphname << "graph_AN_peak_" << beams[iBeam] << "_" << particles[iParticle] << "_" << methods[iMethod] << "_" << trigger_shorts[iTrigger];
          extract_average_asymmetry(
            graph_AN_angleA[iTrigger][iBeam][iParticle][iMethod],
            graphname.str(),
            graph_AN[iTrigger][0][iBeam][iParticle][iMethod],
            graph_AN[iTrigger][1][iBeam][iParticle][iMethod]);
          graph_AN_angleA[iTrigger][iBeam][iParticle][iMethod]->Write();
        }
      }
    }
  }
  // Then average over beam 
  for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
    for (int iMethod = 0; iMethod < nMethods; iMethod++) {
      for (int iParticle = 0; iParticle < nParticles; iParticle++) {
        std::stringstream graphname;
        graphname << "graph_AN_peak_" << particles[iParticle] << "_" << methods[iMethod] << "_" << trigger_shorts[iTrigger];
        extract_average_asymmetry(
          graph_AN_angleA_beamA[iTrigger][iParticle][iMethod],
          graphname.str(),
          graph_AN_angleA[iTrigger][0][iParticle][iMethod],
          graph_AN_angleA[iTrigger][1][iParticle][iMethod]);
        graph_AN_angleA_beamA[iTrigger][iParticle][iMethod]->Write();
      }
    }
  }

  // Average -> Separate Beams
  // Average over triggers and crossing angle
  for (int iCrossingAngle = 0; iCrossingAngle < nCrossingAngles; iCrossingAngle++) {
    for (int iBeam = 0; iBeam < nBeams; iBeam++) {
      for (int iMethod = 0; iMethod < nMethods; iMethod++) {
        for (int iParticle = 0; iParticle < nParticles; iParticle++) {
          std::stringstream graphname;
          graphname << "graph_AN_peak_" << beams[iBeam] << "_" << particles[iParticle] << "_" << methods[iMethod];
          extract_average_asymmetry(
            graph_AN_triggerA_angleA[iBeam][iParticle][iMethod],
            graphname.str(),
            graph_AN_triggerA[0][iBeam][iParticle][iMethod],
            graph_AN_triggerA[1][iBeam][iParticle][iMethod]);
          graph_AN_triggerA_angleA[iBeam][iParticle][iMethod]->Write();
        }
      }
    }
  }
  

  // Full average
  // Average over trigger, crossing angle, and beam
  for (int iMethod = 0; iMethod < nMethods; iMethod++) {
    for (int iParticle = 0; iParticle < nParticles; iParticle++) {
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iParticle] << "_" << methods[iMethod] << "_final";
      extract_average_asymmetry(
        graph_AN_final[iParticle][iMethod],
        graphname.str(),
        graph_AN_triggerA_angleA[0][iParticle][iMethod],
        graph_AN_triggerA_angleA[1][iParticle][iMethod]);
      graph_AN_final[iParticle][iMethod]->Write();
      if (iMethod == 0) {
        std::cout << "eta final : iP = " << iParticle << std::endl;
        graph_AN_final[iParticle][0]->Print();
      }
    }
  }

  // Save final asymmetry plots
  // Plot the two graphs
  for (int iP = 0; iP < nParticles; iP++) {
    std::stringstream graphname;
    graphname << "graph_AN_eta_peak_" << particles[iP];
    std::stringstream canvas_name;
    canvas_name << "canvas_" <<  graphname.str();
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "c");// 1600, 900);
    canvas->cd();
    TMultiGraph *mg = new TMultiGraph();
    std::stringstream graph_title; graph_title << (iP == 0 ? "#pi^{0}" : "#eta") << " Asymmetry";
    mg->SetTitle(graph_title.str().c_str());

    graph_AN_final[iP][0]->SetStats(0);
    graph_AN_final[iP][0]->SetLineColor(kBlue);
    graph_AN_final[iP][0]->SetMarkerColor(kBlue);
    graph_AN_final[iP][1]->SetStats(0);
    graph_AN_final[iP][1]->SetLineColor(kRed);
    graph_AN_final[iP][1]->SetMarkerColor(kRed);
    shift_graph(graph_AN_final[iP][1], 0.05);

    mg->Add(graph_AN_final[iP][0]);
    mg->Add(graph_AN_final[iP][1]);

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("A_{N}");

    double y_bound = GetAbsMaxInRange(graph_AN_final[iP][1], -3.0, 3.0);
    mg->SetMinimum(-2 * y_bound);
    mg->SetMaximum(2 * y_bound);

    canvas->Draw();

    // Add "sPHENIX internal"
    TLatex latex;
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.18, 0.85, "#font[72]{sPHENIX} Internal");
    latex.DrawLatex(0.18, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");

    TLegend *legend = new TLegend(0.18, 0.2, 0.4, 0.4);
    legend->SetTextSize(0.05);
    legend->AddEntry(graph_AN_final[iP][0], "Geometric");
    legend->AddEntry(graph_AN_final[iP][1], "Rel. Lumi");
    legend->Draw();

    gPad->Update();
    gPad->Modified();

    // Show a dashed line at y = 0 (no asymmetry)
    TLine *tline = new TLine();
    tline->SetLineWidth(2);
    tline->SetLineColor(kBlack);
    tline->SetLineStyle(kDashed);
    tline->DrawLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);

    // Save image
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".png").c_str());
    delete canvas;
  }

  // Save final asymmetry plots
  // Plot the two graphs -> Beam Comparison
  for (int iP = 0; iP < nParticles; iP++) {
    std::stringstream graphname;
    graphname << "graph_AN_eta_peak_beams_" << particles[iP];
    std::stringstream canvas_name;
    canvas_name << "canvas_" <<  graphname.str();
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "c");// 1600, 900);
    canvas->cd();
    TMultiGraph *mg = new TMultiGraph();
    std::stringstream graph_title; graph_title << (iP == 0 ? "#pi^{0}" : "#eta") << " Asymmetry";
    mg->SetTitle(graph_title.str().c_str());

    graph_AN_triggerA_angleA[0][iP][0]->SetStats(0);
    graph_AN_triggerA_angleA[0][iP][0]->SetLineColor(kOcean);
    graph_AN_triggerA_angleA[0][iP][0]->SetMarkerColor(kOcean);
    graph_AN_triggerA_angleA[1][iP][0]->SetStats(0);
    graph_AN_triggerA_angleA[1][iP][0]->SetLineColor(kOrange-3);
    graph_AN_triggerA_angleA[1][iP][0]->SetMarkerColor(kOrange-3);
    shift_graph(graph_AN_triggerA_angleA[1][iP][0], 0.05);

    mg->Add(graph_AN_triggerA_angleA[0][iP][0]);
    mg->Add(graph_AN_triggerA_angleA[1][iP][0]);

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("A_{N}");

    double y_bound = GetAbsMaxInRange(graph_AN_triggerA_angleA[1][iP][0], -3.0, 3.0);
    mg->SetMinimum(-2 * y_bound);
    mg->SetMaximum(2 * y_bound);

    canvas->Draw();

    // Add "sPHENIX internal"
    TLatex latex;
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.18, 0.85, "#font[72]{sPHENIX} Internal");
    latex.DrawLatex(0.18, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");

    TLegend *legend = new TLegend(0.18, 0.2, 0.4, 0.4);
    legend->SetTextSize(0.05);
    legend->AddEntry(graph_AN_triggerA_angleA[0][iP][0], "Blue Beam");
    legend->AddEntry(graph_AN_triggerA_angleA[1][iP][0], "Yellow Beam");
    legend->Draw();

    gPad->Update();
    gPad->Modified();

    // Show a dashed line at y = 0 (no asymmetry)
    TLine *tline = new TLine();
    tline->SetLineWidth(2);
    tline->SetLineColor(kBlack);
    tline->SetLineStyle(kDashed);
    tline->DrawLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);

    // Save image
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".png").c_str());
    delete canvas;

    // Compute the t-test
    {
      std::stringstream graphname_t_test; graphname_t_test << graphname.str() << "_test";
      TGraph *t_test = nullptr;
      compute_t_statistic(t_test, graphname_t_test.str(), graph_AN_triggerA_angleA[0][iP][0], graph_AN_triggerA_angleA[1][iP][0]);
      t_test->GetXaxis()->SetTitle("#eta");
      t_test->GetYaxis()->SetTitle("t Statistic");
      std::stringstream canvas_name_t_test; canvas_name_t_test << canvas_name.str() << "_test";
      TCanvas *canvas_t_test = new TCanvas(canvas_name_t_test.str().c_str());
      canvas_t_test->cd();
      t_test->SetLineColor(kBlack);
      t_test->SetMarkerColor(kBlack);
      t_test->SetMarkerStyle(kFullCircle);

      t_test->Draw("AP");

      y_bound = GetAbsMaxInRangeSimple(t_test, -3.0, 3.0);
      t_test->SetMinimum(-y_bound);
      t_test->SetMaximum(y_bound);

      gPad->Update();
      gPad->Modified();
      double min_x = gPad->GetUxmin();
      double max_x = gPad->GetUxmax();

      // Show lines at 0 and +- 1 sigma
      tline->SetLineStyle(kSolid);
      tline->DrawLine(min_x, 0, max_x, 0);
      tline->SetLineStyle(kDashed);
      tline->DrawLine(min_x, -1, max_x, -1);
      tline->DrawLine(min_x, 1, max_x, 1);

      canvas_t_test->Draw();
      canvas_t_test->SaveAs((plots_folder + "/" + canvas_name_t_test.str() + ".pdf").c_str());
      canvas_t_test->SaveAs((plots_folder + "/" + canvas_name_t_test.str() + ".png").c_str());

      delete canvas_t_test;
    }
  }

  // Save final asymmetry plots
  // Plot the two graphs -> Crossing Angle Comparison
  for (int iP = 0; iP < nParticles; iP++) {
    std::stringstream graphname;
    graphname << "graph_AN_eta_peak_angles_" << particles[iP];
    std::stringstream canvas_name;
    canvas_name << "canvas_" <<  graphname.str();
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "c");// 1600, 900);
    canvas->cd();
    TMultiGraph *mg = new TMultiGraph();
    std::stringstream graph_title; graph_title << (iP == 0 ? "#pi^{0}" : "#eta") << " Asymmetry";
    mg->SetTitle(graph_title.str().c_str());

    graph_AN_triggerA_beamA[0][iP][0]->SetStats(0);
    graph_AN_triggerA_beamA[0][iP][0]->SetLineColor(kBlue);
    graph_AN_triggerA_beamA[0][iP][0]->SetMarkerColor(kBlue);
    graph_AN_triggerA_beamA[1][iP][0]->SetStats(0);
    graph_AN_triggerA_beamA[1][iP][0]->SetLineColor(kRed);
    graph_AN_triggerA_beamA[1][iP][0]->SetMarkerColor(kRed);
    //shift_graph(graph_AN_triggerA_beamA[1][iP][0], 0.002);

    mg->Add(graph_AN_triggerA_beamA[0][iP][0]);
    mg->Add(graph_AN_triggerA_beamA[1][iP][0]);

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("A_{N}");

    double y_bound = GetAbsMaxInRange(graph_AN_triggerA_beamA[1][iP][0], -3.0, 3.0);
    mg->SetMinimum(-2 * y_bound);
    mg->SetMaximum(2 * y_bound);

    canvas->Draw();

    // Add "sPHENIX internal"
    TLatex latex;
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.18, 0.85, "#font[72]{sPHENIX} Internal");
    latex.DrawLatex(0.18, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");

    TLegend *legend = new TLegend(0.18, 0.2, 0.4, 0.4);
    legend->SetTextSize(0.05);
    legend->AddEntry(graph_AN_triggerA_beamA[0][iP][0], "0 mrad");
    legend->AddEntry(graph_AN_triggerA_beamA[1][iP][0], "1.5 mrad");
    legend->Draw();

    gPad->Update();
    gPad->Modified();

    // Show a dashed line at y = 0 (no asymmetry)
    TLine *tline = new TLine();
    tline->SetLineWidth(2);
    tline->SetLineColor(kBlack);
    tline->SetLineStyle(kDashed);
    tline->DrawLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);

    // Save image
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".png").c_str());
    delete canvas;

    // Compute the t-test
    {
      std::stringstream graphname_t_test; graphname_t_test << graphname.str() << "_test";
      TGraph *t_test = nullptr;
      compute_t_statistic(t_test, graphname_t_test.str(), graph_AN_triggerA_beamA[0][iP][0], graph_AN_triggerA_beamA[1][iP][0]);
      t_test->GetXaxis()->SetTitle("#eta");
      t_test->GetYaxis()->SetTitle("t Statistic");
      std::stringstream canvas_name_t_test; canvas_name_t_test << canvas_name.str() << "_test";
      TCanvas *canvas_t_test = new TCanvas(canvas_name_t_test.str().c_str());
      canvas_t_test->cd();
      t_test->SetLineColor(kBlack);
      t_test->SetMarkerColor(kBlack);
      t_test->SetMarkerStyle(kFullCircle);

      t_test->Draw("AP");
      y_bound = GetAbsMaxInRangeSimple(t_test, -3.0, 3.0);
      t_test->SetMinimum(-y_bound);
      t_test->SetMaximum(y_bound);

      gPad->Update();
      gPad->Modified();
      double min_x = gPad->GetUxmin();
      double max_x = gPad->GetUxmax();

      // Show lines at 0 and +- 1 sigma
      tline->SetLineStyle(kSolid);
      tline->DrawLine(min_x, 0, max_x, 0);
      tline->SetLineStyle(kDashed);
      tline->DrawLine(min_x, -1, max_x, -1);
      tline->DrawLine(min_x, 1, max_x, 1);

      canvas_t_test->Draw();
      canvas_t_test->SaveAs((plots_folder + "/" + canvas_name_t_test.str() + ".pdf").c_str());
      canvas_t_test->SaveAs((plots_folder + "/" + canvas_name_t_test.str() + ".png").c_str());

      delete canvas_t_test;
    }
  }

  // Save final asymmetry plots
  // Plot the two graphs -> Trigger Comparison
  for (int iP = 0; iP < nParticles; iP++) {
    std::stringstream graphname;
    graphname << "graph_AN_eta_peak_triggers_" << particles[iP];
    std::stringstream canvas_name;
    canvas_name << "canvas_" <<  graphname.str();
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(), "c");// 1600, 900);
    canvas->cd();

    // Trim the MBD graph where statistics are small
    double y_bound = GetAbsMaxInRange(graph_AN_angleA_beamA[1][iP][0], -3.0, 3.0);
    TGraphErrors *graph_mbd_trimmed = nullptr;
    //trim_graph(graph_AN_angleA_beamA[0][iP][0], graph_mbd_trimmed, 1.5 * y_bound);
    cut_graph(graph_AN_angleA_beamA[0][iP][0], graph_mbd_trimmed, 1, 6);
    
    TMultiGraph *mg = new TMultiGraph();
    std::stringstream graph_title; graph_title << (iP == 0 ? "#pi^{0}" : "#eta") << " Asymmetry";
    mg->SetTitle(graph_title.str().c_str());

    graph_mbd_trimmed->SetStats(0);
    graph_mbd_trimmed->SetLineColor(kBlue);
    graph_mbd_trimmed->SetMarkerColor(kBlue);
    graph_AN_angleA_beamA[1][iP][0]->SetStats(0);
    graph_AN_angleA_beamA[1][iP][0]->SetLineColor(kRed);
    graph_AN_angleA_beamA[1][iP][0]->SetMarkerColor(kRed);
    //shift_graph(graph_AN_angleA_beamA[1][iP][0], 0.002);

    mg->Add(graph_mbd_trimmed);
    mg->Add(graph_AN_angleA_beamA[1][iP][0]);

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("A_{N}");

    mg->SetMinimum(-2 * y_bound);
    mg->SetMaximum(2 * y_bound);

    canvas->Draw();

    // Add "sPHENIX internal"
    TLatex latex;
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.18, 0.85, "#font[72]{sPHENIX} Internal");
    latex.DrawLatex(0.18, 0.75, "p^{#uparrow}+p #sqrt{s} = 200 GeV");

    TLegend *legend = new TLegend(0.18, 0.2, 0.4, 0.4);
    legend->SetTextSize(0.05);
    legend->AddEntry(graph_mbd_trimmed, "MBD Trigger");
    legend->AddEntry(graph_AN_angleA_beamA[1][iP][0], "Photon Trigger");
    legend->Draw();

    gPad->Update();
    gPad->Modified();

    // Show a dashed line at y = 0 (no asymmetry)
    TLine *tline = new TLine();
    tline->SetLineWidth(2);
    tline->SetLineColor(kBlack);
    tline->SetLineStyle(kDashed);
    tline->DrawLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);

    // Save image
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
    canvas->SaveAs((plots_folder + "/" + canvas_name.str() + ".png").c_str());
    delete canvas;
  }

  outputfile->Close();
  delete outputfile;
  gSystem->Exit(0);
  
  gSystem->Exit(0);
}
