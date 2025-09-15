import ROOT
tfile = ROOT.TFile("output.root")
tfile.ls()
ttree = tfile.Get("data")
hist = ROOT.TH1D("hits", "Candidate (eta) Mass", 56, 447.5, 670.5)
name = "tag_m"
ttree.Draw("tag_m >> hits")
hist.GetXaxis().SetTitle("Mass")
hist.GetYaxis().SetTitle("Count")
print(hist.GetEntries())
hist.Draw()
ROOT.gPad.Print(name + ".png")
