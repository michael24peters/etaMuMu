import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np

# Specify ntuple leaf to view
leaf = "mcprt_pid"

# open your ROOT file
f = uproot.open("eta2MuMuGamma_mc.root")

# get the tree
tree = f["tree"]        # replace with actual tree name

# load the leaf into a NumPy array
pid = ak.flatten(tree[leaf].array())

# Determine number of bins
# bins = np.arange(min(pid)-0.5, max(pid)+1.5, 1)
bins = np.arange(-13-0.5, 22+1.5, 1)

# plot a histogram
plt.hist(pid, range=[-14,23], edgecolor='black', linewidth=1)
plt.xlabel('Values')
plt.ylabel('Events')
plt.title(leaf + " plot")
plt.savefig("figs/" + leaf + "_focused.png")
plt.close()
