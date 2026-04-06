#!/bin/bash

production_date="01312026"
selection_label="$1"
inputfolder="${PWD}"
inputfolder_mass="../../invariant_mass_fit/parametric_fit/macros"
inputfolder_polarization="/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/polarization_analysis"
inputfolder_analysis="/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/analysis_files/analysis_ana509"
folder_csv_inputs="/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/folder_mass/csv_inputs"
folder_csv_ratios="/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/folder_mass/csv_ratios"
outputfolder_asymmetry="/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/analysis_asymmetry/analysis_ana509/"
plots_folder_template="${inputfolder_mass}"

# Extract average pT/eta/xF in each bin
root -b "${inputfolder}/read_average_bin.C(\"${production_date}\", \"${selection_label}\", \"${inputfolder_mass}\", \"${inputfolder_analysis}\", \"${folder_csv_inputs}\")";

# # Sum forward and backward invariant mass spectra
# root -b "${inputfolder_mass}/macro_sum_forward_backward.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_analysis}\")";

# # Fit the invariant mass spectrum in each bin using the parametric method
# root -b "${inputfolder_mass}/fit_argusmodified_2.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_mass}\",\"${inputfolder_analysis}\",\"${folder_csv_inputs}\", \"${folder_csv_ratios}\",\"${plots_folder_template}\")";
# root -b "${inputfolder_mass}/fit_argusmodified_2_eta.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_mass}\",\"${inputfolder_analysis}\",\"${folder_csv_inputs}\", \"${folder_csv_ratios}\",\"${plots_folder_template}\")";
# root -b "${inputfolder_mass}/fit_argusmodified_2_xf.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_mass}\",\"${inputfolder_analysis}\",\"${folder_csv_inputs}\", \"${folder_csv_ratios}\",\"${plots_folder_template}\")";

# Fit forward and backward pt-dependent asymmetries
root -b "${inputfolder}/fit_asym_direction_all_fills_ana509.C(\"${production_date}\", \"${selection_label}\", true, \"${inputfolder_mass}\", \"${inputfolder_polarization}\",\"${inputfolder_analysis}\", \"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_average_asym_allfills_ana509.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_polarization}\",\"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_sinusoid_asym_direction_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", false)";
root -b "${inputfolder}/fit_sinusoid_asym_direction_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", true)";

# Fit non-directional pt-dependent asymmetries
root -b "${inputfolder}/fit_asym_pt_all_fills_ana509.C(\"${production_date}\", \"${selection_label}\", true, \"${inputfolder_mass}\", \"${inputfolder_polarization}\",\"${inputfolder_analysis}\", \"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_average_asym_pt_allfills_ana509.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_polarization}\",\"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_sinusoid_asym_pt_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", false)";
root -b "${inputfolder}/fit_sinusoid_asym_pt_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", true)";

# Fit eta-dependent asymmetries
root -b "${inputfolder}/fit_asym_eta_all_fills_ana509.C(\"${production_date}\", \"${selection_label}\", true, \"${inputfolder_mass}\", \"${inputfolder_polarization}\",\"${inputfolder_analysis}\", \"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_average_asym_eta_allfills_ana509.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_polarization}\",\"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_sinusoid_asym_eta_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", false)";
root -b "${inputfolder}/fit_sinusoid_asym_eta_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", true)";

# Fit xf-dependent asymmetries
root -b "${inputfolder}/fit_asym_xf_all_fills_ana509.C(\"${production_date}\", \"${selection_label}\", true, \"${inputfolder_mass}\", \"${inputfolder_polarization}\",\"${inputfolder_analysis}\", \"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_average_asym_xf_allfills_ana509.C(\"${production_date}\", \"${selection_label}\",\"${inputfolder_polarization}\",\"${outputfolder_asymmetry}\")";
root -b "${inputfolder}/fit_sinusoid_asym_xf_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", false)";
root -b "${inputfolder}/fit_sinusoid_asym_xf_ana509.C(\"${production_date}\", \"${selection_label}\", \"${outputfolder_asymmetry}\", \"${folder_csv_inputs}\", \"${folder_csv_ratios}\", true)";
