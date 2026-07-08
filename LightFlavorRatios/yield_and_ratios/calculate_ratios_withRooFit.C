#include "TFile.h"
#include "TH1F.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooDataHist.h"
#include "RooPlot.h"

#include <phool/PHRandomSeed.h>
#include <gsl/gsl_rng.h>

#include "../util/binning.h"
#include "../util/DifferentialContainer.h"
#include "../corrections/EfficiencyCorrection.h"
#include "../corrections/GeoAcceptanceCorrection.h"
#include "../corrections/LambdaFeedDownCorrection.h"

void get_yield(TH1F* h_yield, int i, RooDataSet& ds, RooAddPdf model, std::string bkg_name, RooRealVar& mass, RooRealVar& nSignal, RooRealVar& nBackground)
{
  model.fitTo(ds);

  double nsignal = nSignal.getVal();
  double nsignal_err = nSignal.getError();

  double nbkg = nBackground.getVal();
  double nbkg_err = nBackground.getError();

  std::cout << "nsignal val " << nsignal << std::endl;
  std::cout << "nsignal err " << nsignal_err << std::endl;
  std::cout << "nbkg " << nbkg << std::endl;
  std::cout << "nbkg err " << nbkg_err << std::endl;

  h_yield->SetBinContent(i,nsignal);
  h_yield->SetBinError(i,nsignal_err);

  std::string name;
  if(i>=0) name = std::string(h_yield->GetName())+"_"+std::to_string(i);
  else name = std::string(h_yield->GetName())+"_";
  std::string title;
  if(i>=0) title = std::string(h_yield->GetTitle())+" bin "+std::to_string(i);
  else title = std::string(h_yield->GetTitle());

  RooPlot* plot = mass.frame(RooFit::Title(title.c_str()));
  plot->SetName(name.c_str());
  ds.plotOn(plot);
  model.plotOn(plot,RooFit::DrawOption("FL"),RooFit::FillColor(kAzure+1),RooFit::MoveToBack());
  model.plotOn(plot,RooFit::Components(bkg_name.c_str()),RooFit::DrawOption("FL"),RooFit::LineStyle(kDashed),RooFit::FillColor(kGray));
  plot->Write();
}

void get_yield_allbins(TH1F* h_yield, HistogramInfo& hinfo, RooDataSet& ds, RooAddPdf model, std::string bkg_name, RooRealVar& mass, RooRealVar& nSignal, RooRealVar& nBackground)
{
  for(int i=1; i<=h_yield->GetNbinsX(); i++)
  {
    std::cout << "bin " << i << " of " << h_yield->GetNbinsX() << std::endl;
    std::string selection = hinfo.get_bin_selection(std::string(ds.GetName())+"_"+hinfo.name,i);
    std::cout << "selection: " << selection << std::endl;

    RooDataSet ds_selected = static_cast<RooDataSet&>(*(ds.reduce({mass},selection.c_str())));

    get_yield(h_yield,i,ds_selected,model,bkg_name,mass,nSignal,nBackground);
  }
}

void calculate_ratios_withRooFit()
{
  TFile* Ks_file = TFile::Open("/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/Kshort_3runs.root");
  TFile* lambda_file = TFile::Open("/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/Lambda_3runs.root");

  TTree* Ks_tree = (TTree*)Ks_file->Get("DecayTree");
  TTree* lambda_tree = (TTree*)lambda_file->Get("DecayTree");
/*
  TH1F* Ks_mass = (TH1F*)Ks_file->Get("Kshort_mass");
  TH1F* lambda_mass = (TH1F*)lambda_file->Get("lambda_mass");

  std::vector<DifferentialContainer> Ks_distributions =
  {
    DifferentialContainer(Ks_file,"K_S0",BinInfo::final_pt_bins),
    DifferentialContainer(Ks_file,"K_S0",BinInfo::final_rapidity_bins),
    DifferentialContainer(Ks_file,"K_S0",BinInfo::final_phi_bins),
    DifferentialContainer(Ks_file,"K_S0",BinInfo::final_ntrack_bins)
  };

  std::vector<DifferentialContainer> lambda_distributions =
  {
    DifferentialContainer(lambda_file,"Lambda0",BinInfo::final_pt_bins),
    DifferentialContainer(lambda_file,"Lambda0",BinInfo::final_rapidity_bins),
    DifferentialContainer(lambda_file,"Lambda0",BinInfo::final_phi_bins),
    DifferentialContainer(lambda_file,"Lambda0",BinInfo::final_ntrack_bins)
  };
*/
  const float maxval = 10.;
  bool blind = true;

  std::vector<HistogramInfo> yield_hinfo =
  {
    BinInfo::final_pt_bins,
    BinInfo::final_eta_bins,
    BinInfo::final_rapidity_bins,
    BinInfo::final_phi_bins,
  };

  std::vector<TH1F*> Ks_yields;
  std::vector<TH1F*> lambda_yields;
  RooArgList Ks_args;
  RooArgList lambda_args;

  RooRealVar m_ks("K_S0_mass","K_S0_mass",0.45,0.55);
  RooRealVar m_lambda("Lambda0_mass","Lambda0_mass",1.1,1.14);

  std::vector<RooRealVar> Ks_vars;
  std::vector<RooRealVar> lambda_vars;

  Ks_args.add(m_ks);
  lambda_args.add(m_lambda);

  for(HistogramInfo& hinfo : yield_hinfo) 
  {
    Ks_yields.push_back(makeHistogram("Ks_yield","K_{S}^{0} yield",hinfo));
    lambda_yields.push_back(makeHistogram("Lambda_yield","#Lambda yield",hinfo));
    std::string Ks_branchname = "K_S0_"+hinfo.name;
    std::string lambda_branchname = "Lambda0_"+hinfo.name;
    std::cout << Ks_branchname << " " << lambda_branchname << std::endl;
    RooRealVar Ks_var(Ks_branchname.c_str(),Ks_branchname.c_str(),hinfo.bins.front(),hinfo.bins.back());
    RooRealVar lambda_var(lambda_branchname.c_str(),lambda_branchname.c_str(),hinfo.bins.front(),hinfo.bins.back());
    Ks_vars.push_back(Ks_var);
    lambda_vars.push_back(lambda_var);
  }

  for(int i=0; i<Ks_vars.size(); i++)
  {
    Ks_args.add(Ks_vars[i]);
    lambda_args.add(lambda_vars[i]);
  }

  Ks_args.Print();
  lambda_args.Print();

  RooDataSet Ks_ds("K_S0","K_S0",Ks_args,RooFit::Import(*Ks_tree));
  RooDataSet lambda_ds("Lambda0","Lambda0",lambda_args,RooFit::Import(*lambda_tree));

  // Kshort signal and background models
  RooRealVar ks_mean("ks_mean","mean",0.49,0.48,0.5);
  RooRealVar ks_width("ks_width","width",0.01,0.0001,0.2);
  RooGaussian ks_signal("ks_signal","signal",m_ks,ks_mean,ks_width);

  RooRealVar ks_lin_coef("ks_lin_coef","lin_coef",0.,0.,1.);
  RooRealVar ks_quad_coef("ks_quad_coef","quad_coef",0.,0.,1.);
  RooRealVar ks_cubic_coef("ks_cubic_coef","cubic_coef",0.,0.,1.);
  RooRealVar ks_quartic_coef("ks_quartic_coef","quartic_coef",0.,0.,1.);
  RooRealVar ks_quintic_coef("ks_quintic_coef","quintic_coef",0.,0.,1.);
  RooChebychev ks_bkg("ks_bkg","background",m_ks,RooArgList(ks_lin_coef,ks_quad_coef,ks_cubic_coef));

  //RooRealVar ks_lin_coef("ks_a0","a0",0.,-maxval,maxval);
  //RooRealVar ks_quad_coef("ks_a1","a1",0.,-maxval,maxval);
  //RooRealVar ks_cub_coef("ks_a2","a2",0.,-maxval,maxval);
  //RooPolynomial ks_bkg("ks_bkg","background",m_ks,RooArgList(ks_lin_coef,ks_quad_coef,ks_cub_coef));

  RooRealVar ks_nsig("ks_nsig","nsignal",100.,0.,1e9);
  RooRealVar ks_nbkg("ks_nbkg","nbkg",100.,0.,1e9);
  //RooAddPdf ks_massfit("ks_massfit","combined signal and background",RooArgList(ks_bkg,ks_signal),RooArgList(ks_nsig,ks_nbkg));

  //RooRealVar ks_bkgfrac("ks_bkgfrac","background fraction",0.01,0.,1.);
  RooAddPdf ks_model("ks_model","combined signal and background",RooArgList(ks_bkg,ks_signal),RooArgList(ks_nsig,ks_nbkg));

  // Lambda signal and background models
  RooRealVar lambda_mean("lambda_mean","mean",1.113,1.11,1.12);
  RooRealVar lambda_width("lambda_width","width",0.01,0.0001,0.2);
  RooGaussian lambda_signal("lambda_signal","signal",m_lambda,lambda_mean,lambda_width);

  RooRealVar lambda_lin_coef("lambda_lin_coef","lin_coef",0.,0.,1.);
  RooRealVar lambda_quad_coef("lambda_quad_coef","quad_coef",0.,0.,1.);
  RooRealVar lambda_cubic_coef("lambda_cubic_coef","cubic_coef",0.,0.,1.);
  RooRealVar lambda_quartic_coef("lambda_quartic_coef","quartic_coef",0.,0.,1.);
  RooRealVar lambda_quintic_coef("lambda_quintic_coef","quintic_coef",0.,0.,1.);
  RooChebychev lambda_bkg("lambda_bkg","background",m_lambda,RooArgList(lambda_lin_coef,lambda_quad_coef,lambda_cubic_coef));

  //RooRealVar lambda_lin_coef("lambda_a0","a0",1.,-maxval,maxval);
  //RooRealVar lambda_quad_coef("lambda_a1","a1",0.1,-maxval,maxval);
  //RooRealVar lambda_cub_coef("lambda_a2","a2",0.01,-maxval,maxval);
  //RooPolynomial lambda_bkg("lambda_bkg","background",m_lambda,RooArgList(lambda_lin_coef,lambda_quad_coef,lambda_cub_coef));

  //RooRealVar lambda_bkg_frac("lambda_bkg_frac","background fraction",0.5,0.,1.);
  //RooAddPdf lambda_massfit("lambda_massfit","combined signal and background",RooArgList(lambda_bkg,lambda_signal),lambda_bkg_frac);

  RooRealVar lambda_nsig("lambda_nsig","nsignal",100,0.,1e9);
  RooRealVar lambda_nbkg("lambda_nbkg","nbkg",100,0.,1e9);
  //RooRealVar lambda_bkgfrac("lambda_bkgfrac","background fraction",0.1,0.,1.);
  RooAddPdf lambda_model("lambda_model","combined signal and background",RooArgList(lambda_bkg,lambda_signal),RooArgList(lambda_nsig,lambda_nbkg));
  /*
  // Import histograms as RooDataHists
  RooDataHist ks_hist("ks_hist","ks_hist",m_ks,RooFit::Import(*Ks_mass));
  RooDataHist lambda_hist("lambda_hist","lambda_hist",m_lambda,RooFit::Import(*lambda_mass));
  */

  TFile* fout = new TFile("fits.root","RECREATE");

  TH1F* all_Ks = makeHistogram(BinInfo::mass_bins.at("K_S0"));
  TH1F* all_Lambda = makeHistogram(BinInfo::mass_bins.at("Lambda0"));
  get_yield(all_Ks,-1,Ks_ds,ks_model,"ks_bkg",m_ks,ks_nsig,ks_nbkg);
  get_yield(all_Lambda,-1,lambda_ds,lambda_model,"lambda_bkg",m_lambda,lambda_nsig,lambda_nbkg);

  all_Ks->SetTitle("#Lambda^{0} mass;mass [GeV/c^{2}];Candidates");
  all_Lambda->SetTitle("K_{S}^{0} mass;mass [GeV/c^{2}];Candidates");

  for(size_t i=0; i<yield_hinfo.size(); i++)
  {
    std::cout << "======= Differential " << Ks_yields[i]->GetName() << " =======" << std::endl;
    get_yield_allbins(Ks_yields[i],yield_hinfo[i],Ks_ds,ks_model,"ks_bkg",m_ks,ks_nsig,ks_nbkg);
    std::cout << "======= Differential " << lambda_yields[i]->GetName() << " =======" << std::endl;
    get_yield_allbins(lambda_yields[i],yield_hinfo[i],lambda_ds,lambda_model,"lambda_bkg",m_lambda,lambda_nsig,lambda_nbkg);
  }

  // Lambda yield feed-down correction from cascade decays

  std::vector<TH1F*> lambda_yields_corrected;
  for(TH1F* h : lambda_yields)
  {
    std::string newname = std::string(h->GetName())+"_corrected";
    std::string newtitle = std::string(h->GetTitle())+", corrected";
    TH1F* new_h = (TH1F*)h->Clone(newname.c_str());
    new_h->SetTitle(newtitle.c_str());
    lambda_yields_corrected.push_back(new_h);
  }

  std::cout << "CLONED!" << std::endl;

  std::string fd_filename = "/sphenix/tg/tg01/hf/hjheng/HF-analysis/simulation/Pythia_ppMinBias/cascade_feeddown/Cascade_feeddown_fraction.root";
  std::vector<LambdaFeedDownCorrection> fd_corr =
  {
    LambdaFeedDownCorrection(fd_filename,"h_feeddown_frac_xi_all"),
    LambdaFeedDownCorrection(fd_filename,"h_feeddown_frac_xi_eta_all"),
    LambdaFeedDownCorrection(fd_filename,"h_feeddown_frac_xi_rapidity_all"),
    LambdaFeedDownCorrection(fd_filename,"h_feeddown_frac_xi_phi_all")
  };

  std::cout << "GOT CORRECTIONS!" << std::endl;

  for(int idiff=0; idiff<yield_hinfo.size(); idiff++)
  {
    for(int ibin=1; ibin<=lambda_yields[idiff]->GetNbinsX(); ibin++)
    {
      float yield = lambda_yields[idiff]->GetBinContent(ibin);
      float lowedge = lambda_yields[idiff]->GetBinLowEdge(ibin);
      float upedge = lowedge + lambda_yields[idiff]->GetBinWidth(ibin);
      std::cout << yield << ", " << lowedge << ", " << upedge << std::endl;

      std::pair<float,float> corr_and_err = fd_corr[idiff].get_correction_and_error(lowedge,upedge);
      std::cout << corr_and_err.first << ", " << corr_and_err.second << std::endl;
      float corr_yield = yield*(1.-corr_and_err.first);
      float corr_err = corr_yield * sqrt(pow(lambda_yields[idiff]->GetBinError(ibin)/yield,2.)+pow(corr_and_err.second/corr_and_err.first,2.));

      lambda_yields_corrected[idiff]->SetBinContent(ibin,corr_yield);
      lambda_yields_corrected[idiff]->SetBinError(ibin,corr_err);
    }
  }

  TH1F* lambdaKsratio = new TH1F("integrated_lambdaKs_ratio","Integrated #Lambda/K_{S}^{0} Ratio",1,0.,1.);

  std::vector<TH1F*> lambdaKs_diffratios;
  std::vector<TH1F*> lambdaKs_diffratios_corrected;
  for(HistogramInfo& hinfo : yield_hinfo) lambdaKs_diffratios.push_back(makeHistogram("lambdaKsratio","#Lambda/K_{S}^{0} Ratio",hinfo));
  for(HistogramInfo& hinfo : yield_hinfo) lambdaKs_diffratios_corrected.push_back(makeHistogram("lambdaKsratio_corrected","#Lambda/K_{S}^{0} Ratio, corrected",hinfo));

  for(size_t idiff=0; idiff<yield_hinfo.size(); idiff++)
  {
    for(size_t ibin=1; ibin<=lambdaKs_diffratios[idiff]->GetNbinsX(); ibin++)
    {
      float ratio = lambda_yields[idiff]->GetBinContent(ibin)/(2.*Ks_yields[idiff]->GetBinContent(ibin));
      float ratio_err = ratio * sqrt(pow(lambda_yields[idiff]->GetBinError(ibin)/lambda_yields[idiff]->GetBinContent(ibin),2.)+pow(Ks_yields[idiff]->GetBinError(ibin)/Ks_yields[idiff]->GetBinContent(ibin),2.));
      lambdaKs_diffratios[idiff]->SetBinContent(ibin,ratio);
      lambdaKs_diffratios[idiff]->SetBinError(ibin,ratio_err);
      lambdaKs_diffratios_corrected[idiff]->SetBinContent(ibin,ratio);
      lambdaKs_diffratios_corrected[idiff]->SetBinError(ibin,ratio_err);
    }
  }

  EfficiencyCorrection eff_corr; // pt only for now

  std::vector<GeoAcceptanceCorrection> acc_corr = {
    GeoAcceptanceCorrection("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio.root","pT"),
    GeoAcceptanceCorrection("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_eta.root","Lambda0_inGeo_#eta"),
    GeoAcceptanceCorrection("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_rap.root","Lambda0_inGeo_y"),
    GeoAcceptanceCorrection("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_phi.root","Lambda0_inGeo_#phi")
  };

  // apply efficiency correction to pT-dependent ratio
  for(size_t idiff = 0; idiff<yield_hinfo.size(); idiff++)
  {
    for(size_t i=1; i<=lambdaKs_diffratios_corrected[idiff]->GetNbinsX(); i++)
    {
      float ratio = lambdaKs_diffratios_corrected[idiff]->GetBinContent(i);
      float ratio_err = lambdaKs_diffratios_corrected[idiff]->GetBinError(i);
      float xval = lambdaKs_diffratios_corrected[idiff]->GetBinCenter(i);
      float xlow = lambdaKs_diffratios_corrected[idiff]->GetBinLowEdge(i);
      float xhigh = xlow + lambdaKs_diffratios_corrected[idiff]->GetBinWidth(i);

      float eff = 1.;
      float eff_err = 0.;

      if(idiff==0) // eff correction only for vs. pT for now
      {
        eff = eff_corr.get_eff(xval);
        eff_err = eff_corr.get_eff_error(xval);
      }
      std::pair<float,float> acc_and_err = acc_corr[idiff].get_acceptance_and_error(xlow,xhigh);
      float acc = acc_and_err.first;
      float acc_err = acc_and_err.second;

      float corrected_ratio = ratio / (eff*acc);
      float corrected_err = ratio / (eff*acc) * sqrt(pow(ratio_err/ratio,2.)+pow(eff_err/eff,2.)+pow(acc_err/acc,2.));

      std::cout << "ratio " << ratio << " err " << ratio_err << " xval " << xval << " eff " << eff << " err " << eff_err << " acc " << acc << " err " << acc_err << std::endl;
      std::cout << "new val " << corrected_ratio << " err " << corrected_err << std::endl;
      lambdaKs_diffratios_corrected[idiff]->SetBinContent(i,corrected_ratio);
      lambdaKs_diffratios_corrected[idiff]->SetBinError(i,corrected_err);
    }
  }

  if (blind)
  {
    const uint seed = PHRandomSeed();
    std::unique_ptr<gsl_rng> m_rng;
    m_rng.reset(gsl_rng_alloc(gsl_rng_mt19937));
    gsl_rng_set(m_rng.get(), seed);
    float blind_par = fabs(gsl_rng_uniform_pos(m_rng.get()) - 0.5); //Take a value anywhere between -0.5 and 0.5

    for(size_t i=0; i<yield_hinfo.size(); i++)
    {
      for (int j=1;  j<=lambdaKs_diffratios[i]->GetNbinsX(); j++)
      {
        float content = lambdaKs_diffratios[i]->GetBinContent(j);
        float content_corrected = lambdaKs_diffratios_corrected[i]->GetBinContent(j);
        lambdaKs_diffratios[i]->SetBinContent(j, content + blind_par);
        lambdaKs_diffratios_corrected[i]->SetBinContent(j,content_corrected + blind_par);
      }
    }
  }

  fout->cd();

  all_Ks->Write();
  all_Lambda->Write();

  for(size_t i=0; i<yield_hinfo.size(); i++)
  {
    //phi_yields[i]->Write();
    Ks_yields[i]->Write();
    lambda_yields[i]->Write();
    lambda_yields_corrected[i]->Write();

    //phiKs_diffratios[i]->Write();
    lambdaKs_diffratios[i]->Write();
    lambdaKs_diffratios_corrected[i]->Write();
  }

/*
  // Do the fit
  ks_model.fitTo(ks_hist);
  lambda_model.fitTo(lambda_hist);

  TCanvas* ck = new TCanvas("ck","ck",600,600);
  RooPlot* ks_plot = m_ks.frame(RooFit::Title("Kshort mass"));
  ks_hist.plotOn(ks_plot);
  ks_model.plotOn(ks_plot);
  ks_plot->Draw();

  TCanvas* cl = new TCanvas("cl","cl",600,600);
  RooPlot* lambda_plot = m_lambda.frame(RooFit::Title("Lambda mass"));
  lambda_hist.plotOn(lambda_plot);
  lambda_model.plotOn(lambda_plot);
  lambda_plot->Draw();
*/
}
