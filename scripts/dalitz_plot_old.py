import ROOT
tfile = ROOT.TFile("eta2MuMuGamma_mc.root")
tfile.ls()
ttree = tfile.Get("eta2mumugamma")

# Create 2d histogram.
dalitz = ROOT.TH2F(
  "dalitz", "eta->mumugamma; m^2_{gamma,mu+} [GeV]; m^2_{mu+,mu-} [GeV]",
  100, 0, 120,
  100, 0, 120
)

m12s = []
m23s = []

# Populate plot.
for i in range(ttree.GetEntries()):
  ttree.GetEntry(i)
  p1 = p2 = p3 = None

  for j, pid in enumerate(ttree.prt_pid):
    # Calculate invariant mass (in GeV).
    px, py, pz, e = (ttree.prt_px[j] / 1000., ttree.prt_py[j] / 1000., 
                      ttree.prt_pz[j] / 1000., ttree.prt_e[j] / 1000.)
    if pid == 22: p1 = ROOT.TLorentzVector(px, py, pz, e)
    elif pid == -13: p2 = ROOT.TLorentzVector(px, py, pz, e)
    elif pid == 13: p3 = ROOT.TLorentzVector(px, py, pz, e)

  # Fill data entries into plot.
  if None not in (p1, p2, p3):
    # Debug mass ranges.
    m12s.append((p1 + p2).M2())
    m23s.append((p2 + p3).M2())

    m12 = (p1 + p2).M2()
    m23 = (p2 + p3).M2()
    dalitz.Fill(m12, m23)

# Debug mass ranges.
print(f"m12 range: [{min(m12s):.4f}, {max(m12s):.4f}]")
print(f"m23 range: [{min(m23s):.4f}, {max(m23s):.4f}]")

# Generate plot.
ROOT.gStyle.SetOptStat(0) # 0 = none, 1 = entries, 1111 = entries + mean + RMS
c = ROOT.TCanvas()
dalitz.Draw("COLZ")
c.SaveAs("dalitz_all.png")