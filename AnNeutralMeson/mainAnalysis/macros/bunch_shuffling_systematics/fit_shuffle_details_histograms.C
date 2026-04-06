#include "/sphenix/u/virgilemahaut/style/sPhenixStyle_Greg.C"

void fit_shuffle_details_histograms()
{
  SetsPhenixStyle();

  std::string plots_folder = "figures_details_runs_0312";
  gSystem->Exec(("mkdir -p " + plots_folder).c_str());

  std::string inputfilename = "shuffle_1_100000_runs_0312.root";

  TFile *inputfile = TFile::Open(inputfilename.c_str());

  const int nBeams = 2;
  const std::string beams[nBeams] = {"yellow", "blue"};
  const int nConfigs = 4;
  const std::string configs[nConfigs] = {"0mrad_mbd", "0mrad_photon", "15mrad_mbd", "15mrad_photon"};
  const int nParticles = 2;
  const std::string particles[nParticles] = {"pi0", "eta"};
  const int nDirections = 2;
  const std::string directions[nDirections] = {"forward", "backward"};
  const int nPtBins = 9;
  const double pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 20};
  const double pTMeans[nParticles][nPtBins] = {
    {1.665, 2.385, 3.440, 4.406, 5.399, 6.398, 7.397, 11.715},
    {1.665, 2.349, 3.356, 4.435, 5.418, 6.417, 7.424, 8.750, 11.715}
  };
  const int nEtaBins = 8;
  const double etaBins[nEtaBins + 1] = {-2.00, -1.05, -0.86, -0.61, 0.00, 0.61, 0.86, 1.05, 2.00};
  const double etaMeans[nParticles][nEtaBins] = {
    {-1.201, -0.957, -0.744, -0.342, 0.342, 0.744, 0.957, 1.201},
    {-1.192, -0.951, -0.741, -0.344, 0.344, 0.741, 0.951, 1.192}
  };
  const int nXfBins = 8;
  const double xfBins[nXfBins + 1] = {-0.200, -0.048, -0.035, -0.022, 0.00, 0.022, 0.035, 0.048, 0.200};
  const double xfMeans[nParticles][nXfBins] = {
    {-0.0606, -0.0410, -0.0284, -0.0130, 0.0130, 0.0284, 0.0410, 0.0606}, 
    {-0.0641, -0.0414, -0.0282, -0.0132, 0.0132, 0.0282, 0.0414, 0.0641}
  };
  

  TGraphErrors *gr_final_pt_mean[nParticles][nDirections] = {nullptr};
  TGraphErrors *gr_final_pt_sigma[nParticles][nDirections] = {nullptr};
  TGraphErrors *gr_final_eta_mean[nParticles] = {nullptr};
  TGraphErrors *gr_final_eta_sigma[nParticles] = {nullptr};
  TGraphErrors *gr_final_xf_mean[nParticles] = {nullptr};
  TGraphErrors *gr_final_xf_sigma[nParticles] = {nullptr};

  for (int iB = 0; iB < nBeams; iB++) {
    for (int iC = 0; iC < nConfigs; iC++) {
      for (int iP = 0; iP < nParticles; iP++) {
        std::string local_plots_folder = plots_folder + "/" + beams[iB] + "/" + configs[iC] + "/pT/";
        gSystem->Exec(("mkdir -p " + local_plots_folder).c_str());
        {
          std::stringstream graphname_template;
          graphname_template << "graph_shuffle_pt_nodir_"
                             << beams[iB] << "_"
                             << particles[iP] << "_"
                             << configs[iC] << "_";
          std::string graphname_mean = graphname_template.str() + "mean";
          std::string graphname_sigma = graphname_template.str() + "sigma";
          TGraphErrors *graph_mean = new TGraphErrors();
          graph_mean->SetName(graphname_mean.c_str());
          graph_mean->SetTitle(";p_{T} [GeV];#mu");
          TGraphErrors *graph_sigma = new TGraphErrors();
          graph_sigma->SetName(graphname_sigma.c_str());
          graph_sigma->SetTitle(";p_{T} [GeV];#sigma");
          for (int iPt = 0; iPt < nPtBins; iPt++) {
            // h_AN_xf_shuffle_eta_xf_3_sqrt
            std::stringstream hname;
            hname << "h_AN_pt_nodir_shuffle_"
                  << beams[iB] << "_"
                  << particles[iP] << "_pt_"
                  << iPt << "_"
                  << configs[iC] << "_sqrt";
            TH1F *h_AN = (TH1F*)inputfile->Get(hname.str().c_str());
            if (!h_AN) {
              std::cerr << "Error. Could not open histogram " << hname.str() << " in file " << inputfilename << "." << std::endl;
              continue;
            }
            TF1 *fit_amplitude = new TF1("fit_amplitude", "gaus", -4, 4);
            h_AN->Fit(fit_amplitude);

            graph_mean->SetPoint(iPt, pTMeans[iP][iPt], fit_amplitude->GetParameter(1));
            graph_mean->SetPointError(iPt, 0, fit_amplitude->GetParError(1));
            graph_sigma->SetPoint(iPt, pTMeans[iP][iPt], fit_amplitude->GetParameter(2));
            graph_sigma->SetPointError(iPt, 0, fit_amplitude->GetParError(2));

            // Plot histogram
            {
              std::stringstream canvas_name;
              canvas_name << "canvas_" << hname.str();
              TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
              canvas->cd();
              fit_amplitude->SetLineColor(kRed);
              h_AN->Draw("HIST");
              fit_amplitude->Draw("SAME");

              std::stringstream stream;
              TLatex latex;
              latex.SetNDC();
              latex.SetTextColor(kBlue);
              std::cout << "text size: " << latex.GetTextSize() << std::endl;
              latex.SetTextSize(0.04);
              stream << "#mu = " << std::setprecision(2) << fit_amplitude->GetParameter(1) << " #pm " << fit_amplitude->GetParError(1);
              latex.DrawLatex(0.43, 0.25, stream.str().c_str());
              stream.str(""); stream << std::setprecision(2) << "#sigma = " << fit_amplitude->GetParameter(2) << " #pm " << fit_amplitude->GetParError(2);
              latex.DrawLatex(0.43,0.20, stream.str().c_str());
              latex.SetTextColor(kBlack);
              latex.SetTextSize(0.05);
              stream.str(""); stream << std::fixed << std::setprecision(0) << pTBins[iPt] << " < p_{T} < " << pTBins[iPt+1] << " GeV";
              latex.DrawLatex(0.2, 0.80, stream.str().c_str());

              canvas->Draw();
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
              delete canvas;
            }
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_mean->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_mean" << std::endl;
            graph_mean->Print();
            graph_mean->SetLineColor(color);
            graph_mean->SetMarkerColor(color);
            graph_mean->SetMarkerStyle(kFullCircle);
            graph_mean->Draw("AP");

            graph_mean->SetMinimum(-0.05);
            graph_mean->SetMaximum(0.05);

            if (iP == 0) graph_mean->GetXaxis()->SetRangeUser(1.0, 10.0);
            else graph_mean->GetXaxis()->SetRangeUser(2.0, 20.0);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 0, max_x, 0);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_sigma->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_sigma" << std::endl;
            graph_sigma->Print();
            graph_sigma->SetLineColor(color);
            graph_sigma->SetMarkerColor(color);
            graph_sigma->SetMarkerStyle(kFullCircle);
            graph_sigma->Draw("AP");

            graph_sigma->SetMinimum(0.95);
            graph_sigma->SetMaximum(1.05);

            if (iP == 0) graph_sigma->GetXaxis()->SetRangeUser(1.0, 10.0);
            else graph_sigma->GetXaxis()->SetRangeUser(2.0, 20.0);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 1, max_x, 1);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
        }
        for (int iDir = 0; iDir < nDirections; iDir++) {
          std::stringstream graphname_template;
          graphname_template << "graph_shuffle_pt_"
                             << beams[iB] << "_"
                             << particles[iP] << "_"
                             << directions[iDir] << "_"
                             << configs[iC] << "_";
          std::string graphname_mean = graphname_template.str() + "mean";
          std::string graphname_sigma = graphname_template.str() + "sigma";
          TGraphErrors *graph_mean = new TGraphErrors();
          graph_mean->SetName(graphname_mean.c_str());
          graph_mean->SetTitle(";p_{T} [GeV];#mu");
          TGraphErrors *graph_sigma = new TGraphErrors();
          graph_sigma->SetName(graphname_sigma.c_str());
          graph_sigma->SetTitle(";p_{T} [GeV];#sigma");
          for (int iPt = 0; iPt < nPtBins; iPt++) {
            // h_AN_xf_shuffle_eta_xf_3_sqrt
            std::stringstream hname;
            hname << "h_AN_pt_shuffle_"
                  << beams[iB] << "_"
                  << particles[iP] << "_"
                  << directions[iDir] << "_pt_"
                  << iPt << "_"
                  << configs[iC] << "_sqrt";
            TH1F *h_AN = (TH1F*)inputfile->Get(hname.str().c_str());
            if (!h_AN) {
              std::cerr << "Error. Could not open histogram " << hname.str() << " in file " << inputfilename << "." << std::endl;
              continue;
            }
            TF1 *fit_amplitude = new TF1("fit_amplitude", "gaus", -4, 4);
            h_AN->Fit(fit_amplitude);

            graph_mean->SetPoint(iPt, pTMeans[iP][iPt], fit_amplitude->GetParameter(1));
            graph_mean->SetPointError(iPt, 0, fit_amplitude->GetParError(1));
            graph_sigma->SetPoint(iPt, pTMeans[iP][iPt], fit_amplitude->GetParameter(2));
            graph_sigma->SetPointError(iPt, 0, fit_amplitude->GetParError(2));

            // Plot histogram
            {
              std::stringstream canvas_name;
              canvas_name << "canvas_" << hname.str();
              TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
              canvas->cd();
              fit_amplitude->SetLineColor(kRed);
              h_AN->Draw("HIST");
              fit_amplitude->Draw("SAME");

              std::stringstream stream;
              TLatex latex;
              latex.SetNDC();
              latex.SetTextColor(kBlue);
              std::cout << "text size: " << latex.GetTextSize() << std::endl;
              latex.SetTextSize(0.04);
              stream << "#mu = " << std::setprecision(2) << fit_amplitude->GetParameter(1) << " #pm " << fit_amplitude->GetParError(1);
              latex.DrawLatex(0.43, 0.25, stream.str().c_str());
              stream.str(""); stream << std::setprecision(2) << "#sigma = " << fit_amplitude->GetParameter(2) << " #pm " << fit_amplitude->GetParError(2);
              latex.DrawLatex(0.43,0.20, stream.str().c_str());
              latex.SetTextColor(kBlack);
              latex.SetTextSize(0.05);
              stream.str(""); stream << std::fixed << std::setprecision(0) << pTBins[iPt] << " < p_{T} < " << pTBins[iPt+1] << " GeV";
              latex.DrawLatex(0.2, 0.80, stream.str().c_str());

              canvas->Draw();
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
              delete canvas;
            }
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_mean->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_mean" << std::endl;
            graph_mean->Print();
            graph_mean->SetLineColor(color);
            graph_mean->SetMarkerColor(color);
            graph_mean->SetMarkerStyle(kFullCircle);
            graph_mean->Draw("AP");

            graph_mean->SetMinimum(-0.05);
            graph_mean->SetMaximum(0.05);

            if (iP == 0) graph_mean->GetXaxis()->SetRangeUser(1.0, 10.0);
            else graph_mean->GetXaxis()->SetRangeUser(2.0, 20.0);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 0, max_x, 0);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_sigma->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_sigma" << std::endl;
            graph_sigma->Print();
            graph_sigma->SetLineColor(color);
            graph_sigma->SetMarkerColor(color);
            graph_sigma->SetMarkerStyle(kFullCircle);
            graph_sigma->Draw("AP");

            graph_sigma->SetMinimum(0.95);
            graph_sigma->SetMaximum(1.05);

            if (iP == 0) graph_sigma->GetXaxis()->SetRangeUser(1.0, 10.0);
            else graph_sigma->GetXaxis()->SetRangeUser(2.0, 20.0);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 1, max_x, 1);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
        }
        {
          std::string local_plots_folder = plots_folder + "/" + beams[iB] + "/" + configs[iC] + "/eta/";
          gSystem->Exec(("mkdir -p " + local_plots_folder).c_str());
          std::stringstream graphname_template;
          graphname_template << "graph_shuffle_eta_"
                             << beams[iB] << "_"
                             << particles[iP] << "_"
                             << configs[iC] << "_";
          std::string graphname_mean = graphname_template.str() + "mean";
          std::string graphname_sigma = graphname_template.str() + "sigma";
          TGraphErrors *graph_mean = new TGraphErrors();
          graph_mean->SetName(graphname_mean.c_str());
          graph_mean->SetTitle(";#eta;#mu");
          TGraphErrors *graph_sigma = new TGraphErrors();
          graph_sigma->SetName(graphname_sigma.c_str());
          graph_sigma->SetTitle(";#eta;#sigma");
          for (int iEta = 0; iEta < nEtaBins; iEta++) {
            // h_AN_xf_shuffle_eta_xf_3_sqrt
            std::stringstream hname;
            hname << "h_AN_eta_shuffle_"
                  << beams[iB] << "_"
                  << particles[iP] << "_eta_"
                  << iEta << "_"
                  << configs[iC] << "_sqrt";
            TH1F *h_AN = (TH1F*)inputfile->Get(hname.str().c_str());
            if (!h_AN) {
              std::cerr << "Error. Could not open histogram " << hname.str() << " in file " << inputfilename << "." << std::endl;
              continue;
            }
            TF1 *fit_amplitude = new TF1("fit_amplitude", "gaus", -4, 4);
            h_AN->Fit(fit_amplitude);

            graph_mean->SetPoint(iEta, etaMeans[iP][iEta], fit_amplitude->GetParameter(1));
            graph_mean->SetPointError(iEta, 0, fit_amplitude->GetParError(1));
            graph_sigma->SetPoint(iEta, etaMeans[iP][iEta], fit_amplitude->GetParameter(2));
            graph_sigma->SetPointError(iEta, 0, fit_amplitude->GetParError(2));

            // Plot histogram
            {
              std::stringstream canvas_name;
              canvas_name << "canvas_" << hname.str();
              TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
              canvas->cd();
              fit_amplitude->SetLineColor(kRed);
              h_AN->Draw("HIST");
              fit_amplitude->Draw("SAME");

              std::stringstream stream;
              TLatex latex;
              latex.SetNDC();
              latex.SetTextColor(kBlue);
              std::cout << "text size: " << latex.GetTextSize() << std::endl;
              latex.SetTextSize(0.04);
              stream << "#mu = " << std::setprecision(2) << fit_amplitude->GetParameter(1) << " #pm " << fit_amplitude->GetParError(1);
              latex.DrawLatex(0.43, 0.25, stream.str().c_str());
              stream.str(""); stream << std::setprecision(3) << "#sigma = " << fit_amplitude->GetParameter(2) << " #pm " << fit_amplitude->GetParError(2);
              latex.DrawLatex(0.43,0.20, stream.str().c_str());
              latex.SetTextColor(kBlack);
              latex.SetTextSize(0.05);
              stream.str(""); stream << std::fixed << std::setprecision(2) << etaBins[iEta] << " < #eta < " << etaBins[iEta+1];
              latex.DrawLatex(0.2, 0.80, stream.str().c_str());

              canvas->Draw();
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
              delete canvas;
            }
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_mean->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_mean" << std::endl;
            graph_mean->Print();
            graph_mean->SetLineColor(color);
            graph_mean->SetMarkerColor(color);
            graph_mean->SetMarkerStyle(kFullCircle);
            graph_mean->Draw("AP");

            graph_mean->SetMinimum(-0.05);
            graph_mean->SetMaximum(0.05);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 0, max_x, 0);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_sigma->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_sigma" << std::endl;
            graph_sigma->Print();
            graph_sigma->SetLineColor(color);
            graph_sigma->SetMarkerColor(color);
            graph_sigma->SetMarkerStyle(kFullCircle);
            graph_sigma->Draw("AP");

            graph_sigma->SetMinimum(0.95);
            graph_sigma->SetMaximum(1.05);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 1, max_x, 1);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
        }
        {
          std::string local_plots_folder = plots_folder + "/" + beams[iB] + "/" + configs[iC] + "/xF/";
          gSystem->Exec(("mkdir -p " + local_plots_folder).c_str());
          std::stringstream graphname_template;
          graphname_template << "graph_shuffle_xf_"
                             << beams[iB] << "_"
                             << particles[iP] << "_"
                             << configs[iC] << "_";
          std::string graphname_mean = graphname_template.str() + "mean";
          std::string graphname_sigma = graphname_template.str() + "sigma";
          TGraphErrors *graph_mean = new TGraphErrors();
          graph_mean->SetName(graphname_mean.c_str());
          graph_mean->SetTitle(";x_{F};#mu");
          TGraphErrors *graph_sigma = new TGraphErrors();
          graph_sigma->SetName(graphname_sigma.c_str());
          graph_sigma->SetTitle(";x_{F};#sigma");
          for (int iXf = 0; iXf < nXfBins; iXf++) {
            // h_AN_xf_shuffle_xf_xf_3_sqrt
            std::stringstream hname;
            hname << "h_AN_xf_shuffle_"
                  << beams[iB] << "_"
                  << particles[iP] << "_xf_"
                  << iXf << "_"
                  << configs[iC] << "_sqrt";
            TH1F *h_AN = (TH1F*)inputfile->Get(hname.str().c_str());
            if (!h_AN) {
              std::cerr << "Error. Could not open histogram " << hname.str() << " in file " << inputfilename << "." << std::endl;
              continue;
            }
            TF1 *fit_amplitude = new TF1("fit_amplitude", "gaus", -4, 4);
            h_AN->Fit(fit_amplitude);

            graph_mean->SetPoint(iXf, xfMeans[iP][iXf], fit_amplitude->GetParameter(1));
            graph_mean->SetPointError(iXf, 0, fit_amplitude->GetParError(1));
            graph_sigma->SetPoint(iXf, xfMeans[iP][iXf], fit_amplitude->GetParameter(2));
            graph_sigma->SetPointError(iXf, 0, fit_amplitude->GetParError(2));

            // Plot histogram
            {
              std::stringstream canvas_name;
              canvas_name << "canvas_" << hname.str();
              TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
              canvas->cd();
              fit_amplitude->SetLineColor(kRed);
              h_AN->Draw("HIST");
              fit_amplitude->Draw("SAME");

              std::stringstream stream;
              TLatex latex;
              latex.SetNDC();
              latex.SetTextColor(kBlue);
              std::cout << "text size: " << latex.GetTextSize() << std::endl;
              latex.SetTextSize(0.04);
              stream << "#mu = " << std::setprecision(2) << fit_amplitude->GetParameter(1) << " #pm " << fit_amplitude->GetParError(1);
              latex.DrawLatex(0.43, 0.25, stream.str().c_str());
              stream.str(""); stream << std::setprecision(3) << "#sigma = " << fit_amplitude->GetParameter(2) << " #pm " << fit_amplitude->GetParError(2);
              latex.DrawLatex(0.43,0.20, stream.str().c_str());
              latex.SetTextColor(kBlack);
              latex.SetTextSize(0.05);
              stream.str(""); stream << std::fixed << std::setprecision(3) << xfBins[iXf] << " < x_{F} < " << xfBins[iXf+1];
              latex.DrawLatex(0.2, 0.80, stream.str().c_str());

              canvas->Draw();
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());
              canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
              delete canvas;
            }
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_mean->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_mean" << std::endl;
            graph_mean->Print();
            graph_mean->SetLineColor(color);
            graph_mean->SetMarkerColor(color);
            graph_mean->SetMarkerStyle(kFullCircle);
            graph_mean->Draw("AP");

            graph_mean->SetMinimum(-0.05);
            graph_mean->SetMaximum(0.05);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 0, max_x, 0);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
          {
            std::stringstream canvas_name;
            canvas_name << "canvas_" << graph_sigma->GetName();
            TCanvas *canvas = new TCanvas(canvas_name.str().c_str());
            canvas->cd();
            int color = (iP == 0 ? kRed : kBlue);
            std::cout << "Print graph_sigma" << std::endl;
            graph_sigma->Print();
            graph_sigma->SetLineColor(color);
            graph_sigma->SetMarkerColor(color);
            graph_sigma->SetMarkerStyle(kFullCircle);
            graph_sigma->Draw("AP");

            graph_sigma->SetMinimum(0.95);
            graph_sigma->SetMaximum(1.05);

            gPad->Update();
            gPad->Modified();
            double min_x = gPad->GetUxmin();
            double max_x = gPad->GetUxmax();

            TLine *tline = new TLine();
            tline->SetLineWidth(2);
            tline->SetLineColor(kBlack);
            tline->SetLineStyle(kDashed);
            tline->DrawLine(min_x, 1, max_x, 1);

            canvas->Draw();
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".png").c_str());
            canvas->SaveAs((local_plots_folder + "/" + canvas_name.str() + ".pdf").c_str());

            delete canvas;
          }
        }
      }
    }
  }
  gSystem->Exit(0);
}
