#!/bin/bash

declare -a arr_labels=("efficiency_30_pt_g3_0mrad")
declare -a arr_labels=("MBD_pt_l3_0mrad" "MBD_pt_l3_15mrad" "MBD_pt_l4_0mrad" "MBD_pt_l4_15mrad" "efficiency_30_pt_g3_0mrad" "efficiency_30_pt_g3_15mrad" "efficiency_30_pt_g4_0mrad" "efficiency_30_pt_g4_15mrad")

cut_label=""
# cut_label="low_vtx_"
# cut_label="high_vtx_"
# cut_label="low_xf_"
# cut_label="high_xf_"

folder_execution="${PWD}"

pushd ${folder_execution} > /dev/null
for label in "${arr_labels[@]}"; do
    ./execution.sh "${cut_label}${label}"
done
popd > /dev/null
