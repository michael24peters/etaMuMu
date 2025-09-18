# README

To run, make sure you are in the base directory (`etaMuMu/`) and type: `lb-run DaVinci/v45r8 ipython -- scripts/etaMuMu.py -m -g`, where `-m` is the MC data option, `-g` is the eta -> mu+ mu- gamma decay mode option. These are the only options that work right now.

A `.root` file will be generated in the base directory, which is broken into the `tag`, `prt`, `mctag`, `mcprt`, `,gentag`, and `genprt`. `tag` and `prt` are reconstruction-level events; `mc` contains MC-matched events (reconstructed -> generator-level mapping, e.g. "after applying my selection criteria, I think this is an eta which I will put in `tag`, now I'm going to check if it was *actually* an eta particle by finding its closest generator-level match and comparing, the result of which I will put in `mctag`); and `gen` contains all generator-level events which were eta -> mu+ mu- (gamma).

`plots/` contains scripts to plot the leaves within the ntuple, e.g. `mcprt_pid`. These plots should be stored in `figs/`.

`scripts/` contains the scripts used to perform this analysis (fill the ntuple and other accessory operations).

