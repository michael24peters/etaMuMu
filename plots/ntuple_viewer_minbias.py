import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import os
import argparse

# parse command line arguments
parser = argparse.ArgumentParser(description="Plot histogram from ROOT ntuple leaf")
parser.add_argument("leaf", help="Name of leaf to plot")
args = parser.parse_args()

# specify ntuple leaf to view
leaf = args.leaf

infile = "../../private/MC_2018_MinBias/magup/00322941_00000001_1.etamumugamma.root"
outdir = "../figs/mc_minbias/"

# open ROOT file
f = uproot.open(infile)

# get the tree
tree = f["tree"]

# load the leaf into a NumPy array
vals = ak.flatten(tree[leaf].array())

# determine number of bins
#bins = np.arange(-350-0.5, 350+1.5, 1)
bins = np.arange(min(vals)-0.5, max(vals)+1.5, 1)

# plot a histogram
if leaf in ["mctag_pid", "mcprt_pid"]: plt.figure(figsize=(16,8))
# plt.hist(vals, bins=bins, edgecolor='black', linewidth=1)
plt.hist(vals, bins=bins)
plt.xlabel('Values')
plt.ylabel('Events')
plt.title(leaf + " plot")

ADD_LABELS = False
if ADD_LABELS:
    # label filled bins
    counts, edges = np.histogram(vals, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    for c, x in zip(counts, centers):
        if c > 0: plt.text(x, c, f"{int(x)}: {c}", ha="center", va="bottom", fontsize=8, rotation=90)

# save plot to file
plt.savefig(os.path.join(outdir, f"{leaf}_hist.png"))
plt.close()
