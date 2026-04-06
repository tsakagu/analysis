# PPG07 Analysis: $\pi^0/\eta$ Transverse Single-Spin Asymmetries in $p^\uparrow + p$ collisions

This is the analysis code of the PPG07 analysis ([Conf Note](https://www.sphenix.bnl.gov/PublicResults/sPH-CONF-COLDQCD-2025-02), [IAN](https://sphenix-invenio.sdcc.bnl.gov/records/psd2w-6kd75)).

The code is divided into several sections:

  - [Run Selection](runSelection): Scripts for extracting the good run list from the SQL database.
  - [Trigger Bit Selection](triggerBitSelection): Pre-select events for which the photon and MBD scaled trigger bits are fired
  - [Trigger Emulator](triggerEmulator): In photon-triggered events, determine which trigger window (8x8 towers non-overlapping regions) have fired the trigger. Adapted from D. Lis original [trigger emulator](https://github.com/sPHENIX-Collaboration/coresoftware/blob/master/offline/packages/trigger/CaloTriggerEmulator.cc)
  - [Polarization Analyis](polarizationAnalysis): Extraction of fill-wise Polarization and Relative Luminosity (+ average) 
  - [Main Analysis](mainAnalysis):
    - [Pre-Analysis](https://github.com/virgile-mahaut/analysis/blob/anneutralmeson_main/AnNeutralMeson/mainAnalysis/macros/Fun4All_AnNeutralMeson.C): Extract relevant cluster information from the pre-selected triggered events and store them.
    - [Micro-Analysis](https://github.com/virgile-mahaut/analysis/tree/anneutralmeson_main/AnNeutralMeson/mainAnalysis/macros/micro_analysis): Extract candidate photon pair/meson candidates + QA plots, after passing all cuts.
    - [Nano-Analysis](github.com/virgile-mahaut/analysis/tree/anneutralmeson_main/AnNeutralMeson/mainAnalysis/macros/nano_analysis): Extract yields from selected meson candidates (nano-analysis)
    - [Invariant mass spectrum fit](https://github.com/virgile-mahaut/analysis/tree/anneutralmeson_main/AnNeutralMeson/mainAnalysis/invariant_mass_fit) (parametric fit or gaussian process regression). Gaussian Process Regression requires python installation with matplotlib and [Gpytorch](https://gpytorch.ai/)
    - [Macros for computing/fitting/averaging the single-spin asymmetries](https://github.com/virgile-mahaut/analysis/tree/anneutralmeson_main/AnNeutralMeson/mainAnalysis/asymmetry_workflow/analysis_code)
    - [Bunch Shuffling](https://github.com/virgile-mahaut/analysis/tree/anneutralmeson_main/AnNeutralMeson/mainAnalysis/bunch_shuffling): Re-compute (fully-corrected and fully-averaged) asymmetries many times, for different permutations of the spin pattern (+ [systematic macro](https://github.com/virgile-mahaut/analysis/tree/anneutralmeson_main/AnNeutralMeson/mainAnalysis/macros/bunch_shuffling_systematics)
  