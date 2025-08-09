

dcount = [100, 200, 100, 50, 90]
ccount = [80, 160, 50, 20, 70]
mcount = [20, 40, 10, 11, 15]
fcount = [18, 37, 8, 9, 14]

# Plot counters on histogram
import matplotlib.pyplot as plt
import numpy as np

# Convert arrays to numpy arrays
dcount = np.array(dcount)
ccount = np.array(ccount)
mcount = np.array(mcount)
fcount = np.array(fcount)

# X locations
x = np.arange(len(dcount))

# Width of each bar
width = 0.2

# Create the plot
fig, ax = plt.subplots(figsize=(max(12, .01*len(dcount)), 6))

# Plot each cut's counts
bar1 = ax.bar(x, dcount, label='DaughtersCuts', color='steelblue')
bar2 = ax.bar(x, ccount, label='CombinationCut', color='firebrick')
bar3 = ax.bar(x, mcount, label='MotherCut', color='seagreen')
bar4 = ax.bar(x, fcount, label='Final Selection', color='purple')

# Labels, title, legend
ax.set_ylabel('Number of Candidates')
ax.set_title('Candidate Survival at Each Cut Stage')
ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig("cutflow_matplotlib.png", dpi=300)
plt.show()