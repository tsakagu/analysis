#!/bin/bash

declare -a arr_suffixes=("pi0_pT" "pi0_forward_pT" "pi0_backward_pT" "pi0_eta" "pi0_xf" "eta_pT" "eta_forward_pT" "eta_backward_pT" "eta_eta" "eta_xf")
second_suffix="_phenix"

for first_suffix in "${arr_suffixes[@]}"; do
    suffix=${first_suffix}${second_suffix}
    paste asym_sys_calc_${suffix}.txt asym_sys_fit_${suffix}.txt asym_sys_shuffle_${suffix}.txt | awk '{printf "%.5f\n", sqrt($1*$1 + $2*$2 + $3*$3)}' > asym_sys_tot_${suffix}.txt
done
