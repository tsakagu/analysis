void get_an_stats()
{
  TFile *inputfile_pt_nodir_phenix = TFile::Open("full_average_pt_nodir_phenix.root");
  TFile *inputfile_pt_nodir_phenix_gpr = TFile::Open("full_average_pt_nodir_phenix_gpr.root");
  TFile *inputfile_pt_nodir = TFile::Open("full_average_pt_nodir.root");
  TFile *inputfile_pt_nodir_gpr = TFile::Open("full_average_pt_nodir_gpr.root");
  TFile *inputfile_pt = TFile::Open("full_average_pt.root");
  TFile *inputfile_pt_gpr = TFile::Open("full_average_pt_gpr.root");
  TFile *inputfile_eta = TFile::Open("full_average_eta.root");
  TFile *inputfile_eta_gpr = TFile::Open("full_average_eta_gpr.root");
  TFile *inputfile_xf = TFile::Open("full_average_xf.root");
  TFile *inputfile_xf_gpr = TFile::Open("full_average_xf_gpr.root");
  TFile *inputfile_xf_low_vtx = TFile::Open("full_average_xf_low_vtx.root");
  TFile *inputfile_xf_low_vtx_gpr = TFile::Open("full_average_xf_low_vtx_gpr.root");
  TFile *inputfile_xf_high_vtx = TFile::Open("full_average_xf_high_vtx.root");
  TFile *inputfile_xf_high_vtx_gpr = TFile::Open("full_average_xf_high_vtx_gpr.root");
  TFile *inputfile_pt_low_xf = TFile::Open("full_average_pt_low_xf.root");
  TFile *inputfile_pt_low_xf_gpr = TFile::Open("full_average_pt_low_xf_gpr.root");
  TFile *inputfile_pt_high_xf = TFile::Open("full_average_pt_high_xf.root");
  TFile *inputfile_pt_high_xf_gpr = TFile::Open("full_average_pt_high_xf_gpr.root");

  std::string outputfolder="asym_values/";
  gSystem->Exec(("mkdir -p " + outputfolder).c_str());
  
  // Show stat uncertainties
  
  const int nParticles = 2;
  const std::string particles[nParticles] = {"pi0", "eta"};
  const int nDirections = 2;
  const std::string directions[nDirections] = {"forward", "backward"};
  const int nMethods = 2;
  const std::string methods[nMethods] = {"sqrt", "rellum"};
  for (int iP = 0; iP < nParticles; iP++) {
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_pT_low_xf" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_pt_low_xf->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_pT_high_xf" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_pt_high_xf->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_pT" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_pt_nodir->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_pT_phenix" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_pt_nodir_phenix->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    for (int iDir = 0; iDir < nDirections; iDir++) {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_" << directions[iDir] << "_pT" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_" << directions[iDir] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_pt->Get(graphname.str().c_str());
      if (!graph_asym) {
        std::cout << "Error. Cannot find graph " << graphname.str() << std::endl;
      }
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_eta" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_eta->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_xf" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_xf->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_xf_low_vtx" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_xf_low_vtx->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_stat_" << particles[iP] << "_xf_high_vtx" << ".txt";
      std::ofstream out_stat(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym = (TGraphErrors*)inputfile_xf_high_vtx->Get(graphname.str().c_str());
      int N = graph_asym->GetN();
      for (int i = 0; i < N; i++) {
        out_stat << std::fixed << std::setprecision(5) << graph_asym->GetPointX(i) << " & " << graph_asym->GetPointY(i) << " & " << graph_asym->GetErrorY(i) << "\n";
      }
      out_stat.close();
    }
  }

  for (int iP = 0; iP < nParticles; iP++) {
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_pT_low_xf.txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_pt_low_xf->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_pt_low_xf->Get(graphname_sqrt.str().c_str());
      
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_pT_high_xf.txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_pt_high_xf->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_pt_high_xf->Get(graphname_sqrt.str().c_str());
      
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_pT.txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_pt_nodir->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_pt_nodir->Get(graphname_sqrt.str().c_str());
      
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_pT_phenix.txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_pt_nodir_phenix->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_pt_nodir_phenix->Get(graphname_sqrt.str().c_str());
      
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    for (int iDir = 0; iDir < nDirections; iDir++) {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_" << directions[iDir] << "_pT.txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_" << directions[iDir] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_" << directions[iDir] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_pt->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_pt->Get(graphname_sqrt.str().c_str());
      
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_eta" << ".txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_eta->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_eta->Get(graphname_sqrt.str().c_str());
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_xf" << ".txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_xf->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_xf->Get(graphname_sqrt.str().c_str());
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_xf_low_vtx" << ".txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_xf_low_vtx->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_xf_low_vtx->Get(graphname_sqrt.str().c_str());
      int N = graph_asym_rellum->GetN();
      std::cout << "N = " << N << std::endl;
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_calc_" << particles[iP] << "_xf_high_vtx" << ".txt";
      std::ofstream out_sys_calc(outname.str().c_str());
      std::stringstream graphname_rellum;
      graphname_rellum << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      std::stringstream graphname_sqrt;
      graphname_sqrt << "graph_AN_peak_" << particles[iP] << "_sqrt_final";
      TGraphErrors *graph_asym_rellum = (TGraphErrors*)inputfile_xf_high_vtx->Get(graphname_rellum.str().c_str());
      TGraphErrors *graph_asym_sqrt = (TGraphErrors*)inputfile_xf_high_vtx->Get(graphname_sqrt.str().c_str());
      int N = graph_asym_rellum->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_rellum->GetPointY(i) - graph_asym_sqrt->GetPointY(i));
        out_sys_calc << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_calc.close();
    }
  }

  for (int iP = 0; iP < nParticles; iP++) {
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_pT_low_xf.txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";;
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_pt_low_xf->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_pt_low_xf_gpr->Get(graphname.str().c_str());

      std::cout << "files = " << inputfile_pt << ", " << inputfile_pt_gpr << std::endl;
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_pT_high_xf.txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";;
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_pt_high_xf->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_pt_high_xf_gpr->Get(graphname.str().c_str());

      std::cout << "files = " << inputfile_pt << ", " << inputfile_pt_gpr << std::endl;
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_pT_phenix.txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";;
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_pt_nodir_phenix->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_pt_nodir_phenix_gpr->Get(graphname.str().c_str());

      std::cout << "files = " << inputfile_pt << ", " << inputfile_pt_gpr << std::endl;
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_pT.txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";;
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_pt_nodir->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_pt_nodir_gpr->Get(graphname.str().c_str());

      std::cout << "files = " << inputfile_pt << ", " << inputfile_pt_gpr << std::endl;
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    for (int iDir = 0; iDir < nDirections; iDir++) {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_" << directions[iDir] << "_pT.txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_" << directions[iDir] << "_rellum_final";;
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_pt->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_pt_gpr->Get(graphname.str().c_str());

      std::cout << "files = " << inputfile_pt << ", " << inputfile_pt_gpr << std::endl;
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_eta" << ".txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_eta->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_eta_gpr->Get(graphname.str().c_str());
      
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_xf" << ".txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_xf->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_xf_gpr->Get(graphname.str().c_str());
      
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_xf_low_vtx" << ".txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_xf_low_vtx->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_xf_low_vtx_gpr->Get(graphname.str().c_str());
      int N = graph_asym_parametric->GetN();
      std::cout << "N = " << N << std::endl;
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
    {
      std::stringstream outname; outname << outputfolder << "/" << "asym_sys_fit_" << particles[iP] << "_xf_high_vtx" << ".txt";
      std::ofstream out_sys_fit(outname.str().c_str());
      std::stringstream graphname;
      graphname << "graph_AN_peak_" << particles[iP] << "_rellum_final";
      TGraphErrors *graph_asym_parametric = (TGraphErrors*)inputfile_xf_high_vtx->Get(graphname.str().c_str());
      TGraphErrors *graph_asym_gpr = (TGraphErrors*)inputfile_xf_high_vtx_gpr->Get(graphname.str().c_str());
      
      int N = graph_asym_parametric->GetN();
      for (int i = 0; i < N; i++) {
        float difference = std::abs(graph_asym_parametric->GetPointY(i) - graph_asym_gpr->GetPointY(i));
        out_sys_fit << std::fixed << std::setprecision(5) << difference << "\n";
      }
      out_sys_fit.close();
    }
  }
  gSystem->Exit(0);
}

