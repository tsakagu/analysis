#!/bin/bash

declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad")
for label in "${arr_labels[@]}"; do
    ./run_gpr_fits_mbd.sh "${label}";
done

declare -a arr_labels=("efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")
for label in "${arr_labels[@]}"; do
    ./run_gpr_fits_photon.sh "${label}";
done

# declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_mbd.sh "low_vtx_${label}";
# done

# declare -a arr_labels=("efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_photon.sh "low_vtx_${label}";
# done

# declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_mbd.sh "high_vtx_${label}";
# done

# declare -a arr_labels=("efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_photon.sh "high_vtx_${label}";
# done

# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_mbd.sh "low_xf_${label}";
# done

# declare -a arr_labels=("efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_photon.sh "low_xf_${label}";
# done

# declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_mbd.sh "high_xf_${label}";
# done

# declare -a arr_labels=("efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_photon.sh "high_xf_${label}";
# done

# declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_mbd.sh "phenix_${label}";
# done

# declare -a arr_labels=("efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")
# for label in "${arr_labels[@]}"; do
#     ./run_gpr_fits_photon.sh "phenix_${label}";
# done
