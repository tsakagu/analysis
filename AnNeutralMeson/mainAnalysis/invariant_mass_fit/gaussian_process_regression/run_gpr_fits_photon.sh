#!/bin/bash

production_date="01312026"
selection_label="$1"
inputfolder="/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/analysis_files/analysis_ana509/analysis_ana509"
outputfolder="result_figures/analysis_ana509/${production_date}/${selection_label}"
mkdir -p ${outputfolder}

inputfilename="analysis_complete_ana509_${production_date}_${selection_label}.root"
inputfilename_sym="analysis_complete_ana509_${production_date}_${selection_label}_symmetric.root"

declare -a arr_pt_indices=("2" "3" "4" "5" "6" "7" "8");
for pt_index in "${arr_pt_indices[@]}"; do
    echo "Fit pT index ${pt_index} (${production_date}, ${selection_label})"
    histoname="h_pair_mass_pt_${pt_index}"
    outname="canvas_pi0_eta_mass_pt_${pt_index}"
    outname_overlay="canvas_pi0_eta_mass_pt_${pt_index}_pass1"
    outname_root="gpr_fit_mass_pt_${pt_index}.root"
    graphname_template="graph_gpr_mass_fit_pt_${pt_index}";
    python invariant_mass_fit_gpr.py --input ${inputfolder}/${inputfilename} --hist ${histoname} --output ${outputfolder}/${outname} --output-root ${outputfolder}/${outname_root} --graph-template ${graphname_template} --k-pi0-mask 4.0 --k-eta-mask 4.0 --min-pi0-lengthscale 0.08 --min-eta-lengthscale 0.29 --verbose
done

declare -a arr_eta_indices=("0" "1" "2" "3")
for eta_index in "${arr_eta_indices[@]}"; do
    echo "Fit eta index ${eta_index} (${production_date}, ${selection_label})"
    histoname="h_pair_mass_eta_${eta_index}_sym"
    outname="canvas_pi0_eta_mass_eta_${eta_index}"
    outname_overlay="canvas_pi0_eta_mass_eta_${eta_index}_pass1"
    outname_root="gpr_fit_mass_eta_${eta_index}.root"
    graphname_template="graph_gpr_mass_fit_eta_${eta_index}";
    python invariant_mass_fit_gpr.py --input ${inputfolder}/${inputfilename_sym} --hist ${histoname} --output ${outputfolder}/${outname} --output-root ${outputfolder}/${outname_root} --graph-template ${graphname_template} --k-pi0-mask 4.0 --k-eta-mask 4.0 --min-pi0-lengthscale 0.08 --min-eta-lengthscale 0.29 --verbose
done

declare -a arr_xf_indices=("0" "1" "2" "3")
for xf_index in "${arr_xf_indices[@]}"; do
    echo "Fit xf index ${xf_index} (${production_date}, ${selection_label})"
    histoname="h_pair_mass_xf_${xf_index}_sym"
    outname="canvas_pi0_xf_mass_xf_${xf_index}"
    outname_overlay="canvas_pi0_xf_mass_xf_${xf_index}_pass1"
    outname_root="gpr_fit_mass_xf_${xf_index}.root"
    graphname_template="graph_gpr_mass_fit_xf_${xf_index}";
    python invariant_mass_fit_gpr.py --input ${inputfolder}/${inputfilename_sym} --hist ${histoname} --output ${outputfolder}/${outname} --output-root ${outputfolder}/${outname_root} --graph-template ${graphname_template} --k-pi0-mask 4.0 --k-eta-mask 4.0 --min-pi0-lengthscale 0.08 --min-eta-lengthscale 0.29 --verbose
done
