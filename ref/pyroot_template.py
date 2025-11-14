import ROOT
tfile = ROOT.TFile("output.root")
tfile.ls()
ttree = tfile.Get("data")
hist = ROOT.TH1D("hits", "Candidate (eta) Mass", 75, 400, 700)
name = "tag_m"
ttree.Draw("tag_m >> hits")
hist.GetXaxis().SetTitle("Mass")
hist.GetYaxis().SetTitle("Count")
hist.SetMaximum(hist.GetMaximum() * 1.2)  # Set y-axis max to 120% of the maximum bin value
print(hist.GetEntries())
hist.Draw()
ROOT.gPad.Print(name + ".png")

