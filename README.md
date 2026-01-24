# README

To run, make sure you are in the base directory (`etaMuMu/`) and type: `lb-run DaVinci/v45r8 ipython scripts/etaMuMu.py`. There are three boolean values in `etaMuMu.py` that allow you to set the data type (MC vs real), data source (local vs analysis production), and decay type (eta -> mu+ mu- (gamma)).

A `.root` file will be generated in the `ntuples/` directory, which is broken into the `tag`, `prt`, and `mc`. `tag` and `prt` are reconstruction-level events; `mc` contains generator-level events. There are a series of indexes which match reco-level events to gen-level events for the purposes of background analysis.

`plots/` contains scripts to plot values int the ntuple, e.g. `mc_pid`. These plots should be stored in `figs/`.

`src/` contains the scripts used to perform this analysis (fill the ntuple and other accessory operations).
