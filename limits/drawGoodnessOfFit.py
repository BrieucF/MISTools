import argparse
import numpy as np
import matplotlib.pyplot as plt

import ROOT as R


parser = argparse.ArgumentParser(description='Draw goodness of fit results')

parser.add_argument('toys', action='store', type=str, help='Goodness of fit outputs on toys')
parser.add_argument('data', action='store', type=str, help='Goodness of fit output on data')

args = parser.parse_args()

def tree_to_dataset(tree):
    entries = tree.GetEntries()

    dataset = []
    for i in range(0, entries):
        tree.GetEntry(i)
        dataset.append(tree.limit)

    return np.asarray(dataset)

# Load toys results
f = R.TFile.Open(args.toys)
toys = tree_to_dataset(f.Get("limit"))

f = R.TFile.Open(args.data)
data = tree_to_dataset(f.Get("limit"))

# Create a figure instance
fig = plt.figure(1, figsize=(7, 7), dpi=300)

# Create an axes instance
ax = fig.add_subplot(111)
ax.margins(0.05, 0.05)
ax.minorticks_on()

# Bin toys data
# toys, bin_edges = np.histogram(toys, bins=10)

# Draw toys (markers)
toys, bin_edges, patches = ax.hist(toys, bins=15)
# ax.plot((bin_edges[1:] + bin_edges[:-1]) / 2, toys, 'o', label='Toys')

ax.errorbar((bin_edges[1:] + bin_edges[:-1]) / 2, toys, yerr=np.sqrt(toys), markersize=0, linewidth=0, elinewidth=1, capsize=3)

# ax.arrow(data[0], 1, 0, -1, head_width=0.05, length_includes_head=True)
ax.annotate('data', xytext=(data[0], max(toys) * 0.2), xy=(data[0], 0), ha='center', arrowprops=dict(arrowstyle='->', mutation_scale=35, connectionstyle='arc3'))

ax.set_xlabel("Best fit test statistic")

fig.tight_layout()
fig.savefig('goodnessoffit.pdf')

