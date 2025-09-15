# =============================================================================

# Count up number of candidate particles
def count(path, counter):
  try:
      count = len(tes[path])
      counter.append(count)
  except: counter.append(0)

# =============================================================================

# Count up number of MC matched candidates
def count_mc_matched(path, counter):
  try:
    prts = tes[path]
    match_count = 0
    for prt in prts:
      mcp = genTool.relatedMCP(prt)
      if mcp and abs(mcp.particleID().pid()) == 221:
        match_count += 1
        print("test")
    counter.append(match_count)
  except: counter.append(0)

# =============================================================================

# Plot counters
def plot_counters(filename, counters):
  # Plot counters on histogram
  import matplotlib.pyplot as plt
  import numpy as np

  # Convert arrays to numpy arrays
  # dcount = np.array(counters[0])
  ccount = np.array(counters[1])
  mcount = np.array(counters[2])
  fcount = np.array(counters[3])

  # Keep entries where not all arrays are 0
  # i.e. remove dead bins
  # mask = ~((dcount == 0) & (ccount == 0) & (mcount == 0) & (fcount == 0))
  mask = ~((ccount == 0) & (mcount == 0) & (fcount == 0))
  # dcount = dcount[mask]
  ccount = ccount[mask]
  mcount = mcount[mask]
  fcount = fcount[mask]

  # X locations
  x = np.arange(len(fcount))

  # Width of each bar
  width = 1.0

  # Create the plot
  fig, ax = plt.subplots(figsize=(max(12, .01*len(fcount)), 6))

  # Plot each cut's counts
  # bar1 = ax.bar(x, dcount, width=width, label='DaughtersCut', color='blue')
  bar2 = ax.bar(x, ccount, width=width-.2, label='CombinationCut', color='firebrick')
  bar3 = ax.bar(x, mcount, width=width-.4, label='MotherCut', color='seagreen')
  bar4 = ax.bar(x, fcount, width=width-.6, label='Final Selection', color='purple')

  # Labels, title, legend
  ax.set_ylabel('Candidates')
  ax.set_title('Candidate Survival at Each Cut Stage')
  ax.legend(loc='upper right')

  plt.tight_layout()
  plt.savefig(filename + ".png", dpi=300)
  plt.show()

# =============================================================================