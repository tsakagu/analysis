# Light-Flavor Ratio Analysis

Repository for determination of ratios of light-flavor resonances, originally constructed for PPG 16 (Lambda/2Kshort ratio).

## Input

The current codebase uses two KFParticle output files, one for Kshort candidates and one for Lambda candidates.
The locations of these files are, for now, `/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/Kshort_3runs.root` and 
`/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/Lambda_3runs.root`.

## Usage

### Running standard analysis

In `yield_and_ratios`:

`root -l calculate_ratios_withRooFit.C`

This runs on the aforementioned input files and generates the file `fits.root` containing the mass peak fits, yield histograms, and ratio histograms.

To generate plots:

In `yield_and_ratios`:

`root -l -b plot_results.C`

This function extracts the histograms from `fits.root`. 
Inside this macros there is a variable `finalize`. 
If set to `true`, it saves the plots to `/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/plots`.
If set to `false`, it saves the plots to `plots/` in your local directory.
(Recommended best practice, if you modify the plotting macro, is to turn `finalize` to `false` until you're certain you like the results.)


PDF versions of the plots are saved to the `pdf/` subfolder, and PNG versions are saved to `png/`.

### Changing binning scheme

Modify the file `util/binning.h`. The static structures at the bottom of this file are the bins used in the analysis. 
Options for linear, logarithmic, and custom binning are provided.

### Adding or modifying corrections

The `corrections/` folder holds the structures and settings through which various corrections are applied. 
Currently, this only contains the relative tracking efficiency correction, which uses Tony's 7-bin histogram, only applied to pT differential results.

## Caveats/To-Dos

- Macros repository is almost certainly outdated and needs to be synchronized with latest light-flavor production
- Unbinned fit will not work for larger datasets, needs to be supplanted by a binned-fit option.
- The plots of the individual mass peaks may have very bad-looking jumps if the background function becomes negative. 
As far as I can tell, this is because RooFit doesn't like when that happens.
Working on finding a fix -- for now, just fiddle with the background function a bit until this doesn't happen.
- Background model tends to misbehave at the edges of the mass window.
