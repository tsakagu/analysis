#!/bin/bash

declare -a arr_suffixes=("pi0_pT" "pi0_forward_pT" "pi0_backward_pT" "pi0_eta" "pi0_xf" "eta_pT" "eta_forward_pT" "eta_backward_pT" "eta_eta" "eta_xf")
second_suffix="_phenix"

for first_suffix in "${arr_suffixes[@]}"; do
    suffix=${first_suffix}${second_suffix}
    paste -d'&' asym_stat_${suffix}.txt norm_asym_sys_shuffle_${suffix}.txt | awk -F"&" '{printf "%.5f\n", ($3*$4)}' > asym_sys_shuffle_${suffix}.txt
done
