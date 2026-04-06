#ifndef __SHUFFLE_BUNCHES_H__
#define __SHUFFLE_BUNCHES_H__

#include "AsymmetryCalc/Constants.h"
#include "AsymmetryCalc/WrapperArrays.h"
#include "AsymmetryCalc/WrapperAsymArrays.h"
#include "AsymmetryCalc/WrapperCorrAsymArrays.h"
#include "AsymmetryCalc/WrapperFitAsymArrays.h"


#include <map>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>

namespace AsymmetryCalc
{
  class ShuffleBunches {
    
  public:

    ShuffleBunches(const int iMin,
                   const int iMax,
                   const std::string& inputfilename_fill_run,
                   const std::string& inputfilename_spin_pattern,
                   const std::string& inputfile_template,
                   const std::string& input_csv_pt_pi0_name,
                   const std::string& input_csv_eta_pi0_name,
                   const std::string& input_csv_xf_pi0_name,
                   const std::string& input_csv_pt_eta_name,
                   const std::string& input_csv_eta_eta_name,
                   const std::string& input_csv_xf_eta_name,
                   const std::string& input_csv_pt_ratios_pi0_name,
                   const std::string& input_csv_eta_ratios_pi0_name,
                   const std::string& input_csv_xf_ratios_pi0_name,
                   const std::string& input_csv_pt_ratios_eta_name,
                   const std::string& input_csv_eta_ratios_eta_name,
                   const std::string& input_csv_xf_ratios_eta_name);
                   
    void get_average_bin_pt_pi0(const std::string &input_csv_name);
    void get_average_bin_eta_pi0(const std::string &input_csv_name);
    void get_average_bin_xf_pi0(const std::string &input_csv_name);

    void get_average_bin_pt_eta(const std::string &input_csv_name);
    void get_average_bin_eta_eta(const std::string &input_csv_name);
    void get_average_bin_xf_eta(const std::string &input_csv_name);

    void get_pt_ratios_pi0(const std::string &input_csv_name);
    void get_eta_ratios_pi0(const std::string &input_csv_name);
    void get_xf_ratios_pi0(const std::string &input_csv_name);

    void get_pt_ratios_eta(const std::string &input_csv_name);
    void get_eta_ratios_eta(const std::string &input_csv_name);
    void get_xf_ratios_eta(const std::string &input_csv_name);

    ~ShuffleBunches();

    void reserve_fill_yields();
    void reserve_fill_raw_asyms();
    void reserve_average_raw_asyms();
    void reserve_corr_asyms();
    void reserve_fit_asyms();

    void clear_fill_yields();
    void clear_fill_raw_asyms();
    void clear_average_raw_asyms();
    void clear_corr_asyms();
    void clear_fit_asyms();
    
    void reset();

    //! Shuffle bunches and get raw asymmetries fill by fill
    void run_1();

    //! Average raw asymmetries and get a single value of corrected asymmetry per iteration
    void run_2();
    
    void compute_fill(const int ifill);
    
    void read_binary_array(const std::string& inputfilename);

    void check_array();
    
    void shuffle_bunch_indices(const int seednumber);

    void sum_yields(const int iIter);

    void set_shuffle(bool val) { do_shuffle = val; }

    void set_store_fill_histos(bool val, const std::string& outputfilename_template);

    void set_store_iter_histos(bool val, const std::string& outputfilename_template);

    void get_yield_histograms();

    void get_fill_raw_asym_histograms();

    void book_outfile_fill_histograms(const std::string& outputfilename);
    
    void save_outfile_fill_histograms();

    void get_iter_raw_asym_histograms();

    void get_iter_corr_asym_histograms();

    void get_iter_fit_asym_histograms();

    void book_outfile_iter_histograms(const std::string& outputfilename);
    
    void save_outfile_iter_histograms();

    void compute_raw_asyms(const int iIter);

    void average_asyms(const int iter);

    void fit_raw_asyms(const int iter);

    void fit_asym(double& asym_amplitude,
                  double& asym_amplitude_err,
                  const double (&asym_azimuth_values)[ASYM_CONSTANTS::nPhiBins],
                  const double (&asym_azimuth_uncertainties)[ASYM_CONSTANTS::nPhiBins],
                  bool fit_geom = true);      

    void compute_corrected_asyms(const int iter);

    void compute_corrected_asym(double &corrected_asym,
                                double &corrected_asym_err,
                                const double &raw_asym,
                                const double &bkg_asym,
                                const double &raw_asym_err,
                                const double &bkg_asym_err,
                                const double &R,
                                const double &R_err);

    void compute_rellum_asym(const uint32_t (&array_yield)[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins],
                             const double pol,
                             double (&array_asym)[ASYM_CONSTANTS::nPhiBins],
                             double (&array_asym_err)[ASYM_CONSTANTS::nPhiBins]);

    void compute_sqrt_asym(const uint32_t (&array_yield)[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins],
                           const double pol,
                           double (&array_asym)[ASYM_CONSTANTS::nPhiBins],
                           double (&array_asym_err)[ASYM_CONSTANTS::nPhiBins]);

    const std::vector<pt_corr_array>& get_corrected_mean_pt() { return corrected_mean_pt; }
    const std::vector<pt_corr_array>& get_corrected_unc_pt() { return corrected_unc_pt; }
    const std::vector<pt_nodir_corr_array>& get_corrected_mean_pt_nodir() { return corrected_mean_pt_nodir; }
    const std::vector<pt_nodir_corr_array>& get_corrected_unc_pt_nodir() { return corrected_unc_pt_nodir; }
    const std::vector<eta_corr_array>& get_corrected_mean_eta() { return corrected_mean_eta; }
    const std::vector<eta_corr_array>& get_corrected_unc_eta() { return corrected_unc_eta; }
    const std::vector<xf_corr_array>& get_corrected_mean_xf() { return corrected_mean_xf; }
    const std::vector<xf_corr_array>& get_corrected_unc_xf() { return corrected_unc_xf; }

  private:

    int iterMin = 0;
    int iterMax = 0;
    int nIterations = 0;
    int nFills = 0;
    
    std::string infile_template = "";

    std::vector<std::pair<int, std::vector<int>>> fill_to_runs;

    // Read the Spin Pattern
    TFile *inputfile_spin = nullptr;
    TTree *spin_patterns = nullptr;
    int fill_spin;
    double blue_polarization;
    double yellow_polarization;
    double polarization[2];
    int8_t yellow_spin_pattern[ASYM_CONSTANTS::nBunches];
    int8_t blue_spin_pattern[ASYM_CONSTANTS::nBunches];
    int8_t spin_pattern[2][ASYM_CONSTANTS::nBunches];
    int8_t shuffled_spin_pattern[2][ASYM_CONSTANTS::nBunches];
    double relative_luminosity = 1;
    
    bool do_shuffle = false;
    int shuffled_indices[ASYM_CONSTANTS::nBunches] = {0};
    bool spin_flip = false;
    int global_irun = 0;

    uint32_t array_yield_pt_bunches[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins][120];
    uint32_t array_yield_eta_bunches[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins][120];
    uint32_t array_yield_xf_bunches[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins][120];

    uint32_t array_yield_pt[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins] = {0};
    uint32_t array_yield_eta[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins] = {0};
    uint32_t array_yield_xf[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins] = {0};

    std::vector<pt_yield_array> vec_yield_pt;
    std::vector<eta_yield_array> vec_yield_eta;
    std::vector<xf_yield_array> vec_yield_xf;

    std::vector<pt_asym_array> fill_pt_asyms;
    std::vector<pt_nodir_asym_array> fill_pt_nodir_asyms;
    std::vector<eta_asym_array> fill_eta_asyms;
    std::vector<xf_asym_array> fill_xf_asyms;

    std::vector<pt_asym_array> fill_pt_asym_errs;
    std::vector<pt_nodir_asym_array> fill_pt_nodir_asym_errs;
    std::vector<eta_asym_array> fill_eta_asym_errs;
    std::vector<xf_asym_array> fill_xf_asym_errs;

    std::vector<pt_asym_array> mean_pt;
    std::vector<pt_nodir_asym_array> mean_pt_nodir;
    std::vector<eta_asym_array> mean_eta;
    std::vector<xf_asym_array> mean_xf;

    std::vector<pt_asym_array> unc_pt;
    std::vector<pt_nodir_asym_array> unc_pt_nodir;
    std::vector<eta_asym_array> unc_eta;
    std::vector<xf_asym_array> unc_xf;

    std::vector<pt_fit_array> fit_mean_pt;
    std::vector<pt_nodir_fit_array> fit_mean_pt_nodir;
    std::vector<eta_fit_array> fit_mean_eta;
    std::vector<xf_fit_array> fit_mean_xf;

    std::vector<pt_fit_array> fit_unc_pt;
    std::vector<pt_nodir_fit_array> fit_unc_pt_nodir;
    std::vector<eta_fit_array> fit_unc_eta;
    std::vector<xf_fit_array> fit_unc_xf;

    std::vector<pt_corr_array> corrected_mean_pt;
    std::vector<pt_nodir_corr_array> corrected_mean_pt_nodir;
    std::vector<eta_corr_array> corrected_mean_eta;
    std::vector<xf_corr_array> corrected_mean_xf;

    std::vector<pt_corr_array> corrected_unc_pt;
    std::vector<pt_nodir_corr_array> corrected_unc_pt_nodir;
    std::vector<eta_corr_array> corrected_unc_eta;
    std::vector<xf_corr_array> corrected_unc_xf;

    // Mean kinematic values
    double pTMeans[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins] = {0};
    double etaMeans[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins] = {0};
    double xfMeans[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins] = {0};

    // background ratios
    double bkg_ratio_pt[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins] = {0};
    double bkg_ratio_err_pt[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins] = {0};
    double bkg_ratio_eta[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins] = {0};
    double bkg_ratio_err_eta[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins] = {0};
    double bkg_ratio_xf[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins] = {0};
    double bkg_ratio_err_xf[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins] = {0};
    
    // Store fill-dependent histograms (NOT for bunch shuffling)
    bool store_fill_histos = false;
    std::string outfilename_fill_template = "";
    TFile *outfile_fill_histograms = nullptr;

    bool store_iter_histos = false;
    std::string outfilename_iter_template = "";
    TFile *outfile_iter_histograms = nullptr;

    // Verification histograms
    TH1F *h_yield_pt[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nSpins] = {nullptr}; // forward vs backward
    TH1F *h_yield_eta[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nSpins] = {nullptr};
    TH1F *h_yield_xf[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nSpins] = {nullptr};
    
  };
};

#endif
