#!/bin/bash

declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad" "efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")

inputfolder_mass="/sphenix/user/virgilemahaut/analysis/AnNeutralMeson/folder_mass"
inputfolder_parametric="../parametric_fit/"
inputfolder_parametric_pt="${inputfolder_parametric}/plots_mass_argusmodified_2/analysis_ana509/ana509_01312026_sigma_3/"
inputfolder_parametric_eta="${inputfolder_parametric}/plots_mass_argusmodified_2_eta/analysis_ana509/ana509_01312026_sigma_3/"
inputfolder_parametric_xf="${inputfolder_parametric}/plots_mass_argusmodified_2_xf/analysis_ana509/ana509_01312026_sigma_3/"
inputfolder_gpr="$PWD"
inputfolder_gpr_fit="${inputfolder_gpr}/result_figures/analysis_ana509/01312026/"
outputfolder="${inputfolder_gpr}/result_comparison/analysis_ana509/01312026/"
mkdir -p ${outputfolder}

cut_label=""
# cut_label="phenix_"
# cut_label="low_vtx_"
# cut_label="high_vtx_"
# cut_label="low_xf_"
# cut_label="high_xf_"

for label in "${arr_labels[@]}"; do
    echo "outputfolder_base = ${outputfolder}"
    root -b "compare_pt_fit_ratios.C(\"${cut_label}${label}\",\"${inputfolder_mass}\",\"${inputfolder_parametric_pt}\",\"${inputfolder_gpr_fit}\",\"${outputfolder}\")"
    root -b "compare_eta_fit_ratios.C(\"${cut_label}${label}\",\"${inputfolder_mass}\",\"${inputfolder_parametric_eta}\",\"${inputfolder_gpr_fit}\",\"${outputfolder}\")"
    root -b "compare_xf_fit_ratios.C(\"${cut_label}${label}\",\"${inputfolder_mass}\",\"${inputfolder_parametric_xf}\",\"${inputfolder_gpr_fit}\",\"${outputfolder}\")"
done

outputfolder_csv="csv_ratios"

for label in "${arr_labels[@]}"; do
    root -b "write_csv_ratios.C(\"${cut_label}${label}\",\"01312026\",\"${inputfolder_gpr}\",\"${outputfolder_csv}\")"
done
