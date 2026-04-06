void compare_gpr_fit()
{
  TFile *file_par = TFile::Open("full_average_xf_high_vtx.root");
  TFile *file_gpr = TFile::Open("full_average_xf_high_vtx_gpr.root");

  std::string canvas_name = "canvas_par_vs_gpr";
  TCanvas *canvas = new TCanvas(canvas_name.c_str());
  canvas->cd();
  std::string graph_name = "graph_AN_peak_blue_pi0_sqrt_efficiency_30";
  TGraphErrors *h_par = (TGraphErrors*)file_par->Get(graph_name.c_str());
  TGraphErrors *h_gpr = (TGraphErrors*)file_gpr->Get(graph_name.c_str());

  h_par->SetMarkerColor(kRed);
  h_par->SetLineColor(kRed);
  h_gpr->SetMarkerColor(kViolet);
  h_gpr->SetLineColor(kViolet);

  h_par->SetMinimum(-0.006);
  h_par->SetMaximum(0.006);
  h_par->Draw("AP E1");
  h_gpr->Draw("P E1 SAME");

  canvas->Draw();
  canvas->SaveAs((canvas_name + ".pdf").c_str());
  canvas->SaveAs((canvas_name + ".png").c_str());

}
