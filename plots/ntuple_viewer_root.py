import ROOT
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description="Plot histogram from ROOT ntuple leaf")
parser.add_argument("leaf", help="Name of leaf to plot")
args = parser.parse_args()

# Specify ntuple leaf to view
leaf = args.leaf

# Open ROOT file and get tree
file = ROOT.TFile.Open("eta2MuMuGamma_mc.root")
#file.ls() # print TFile structure
tree = file.Get("tree")

# Define histogram
hist = ROOT.TH1F("hits", f"{leaf} plot", 38, -14.5, 23.5)
hist.GetXaxis().SetTitle("Values")
hist.GetYaxis().SetTitle("Events")

# Fill histogram directly from tree
tree.Draw(f"{leaf} >> hits")

# Draw plot and save to file
c = ROOT.TCanvas("c", "", 800, 600)
hist.Draw()
c.Update()
c.Print(f"figs/{leaf}_focused.png")

