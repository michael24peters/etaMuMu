#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu analysis                                #
# Author: Michael Peters                                                      #
###############################################################################

# Data.
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['00169948_00000003_7.AllStreams.dst',
                            # '00169948_00000138_7.AllStreams.dst',
                            # '00169948_00000186_7.AllStreams.dst',
                            # '00169948_00000319_7.AllStreams.dst',
                            # '00169948_00000411_7.AllStreams.dst'
                            ],
                            clear = True)

# DaVinci configuration.
from Configurables import DaVinci
DaVinci().DataType = '2018'
# N = cross-section * luminosity
DaVinci().Lumi = False # Processing of luminosity data. Set to False for MC.
DaVinci().Simulation = True # True = MC simulation data.
# Found using lb-dirac dirac-bookkeeping-production-information 00169948.
DaVinci().DDDBtag = 'dddb-20210528-8'
DaVinci().CondDBtag = 'sim-20201113-8-vc-md100-Sim10'

# Output file.
outfile = 'etaMuMuGamma.root'

# Tag configuration.
TrkCats = [('ve', 1), ('tt', 2), ('it', 3), ('ot', 4), ('mu', 7)]
TrgLocs = ['Hlt2ExoticaPrmptDiMuonSSTurbo', 
           'Hlt2ExoticaPrmptDiMuonTurbo',
           'Hlt2ExoticaDiMuonNoIPTurbo',
           'Hlt2ExoticaDisplDiMuon']

# Reconstruction.
from Configurables import CombineParticles, FilterDesktop
from StandardParticles import StdLooseMuons as muons
# TODO: only loose photon cut applied
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

# Apply daughter cuts
sel_dtrs = Selection(
  "SelDaughters",
  Algorithm = FilterDesktop('DaughterCuts', Code='ALL'),
  RequiredSelections=[muons, photons])

# Apply combination cuts
comb = CombineParticles(
  'combEtaMuMuGamma',
  DecayDescriptor= "eta -> mu+ mu- gamma",
  DaughtersCuts = {
    # IPCHI2 3.2 ~= IP < 18um
    "mu+" : "(PT > 500*MeV) & (P > 3*GeV)", 
    "mu-" : "(PT > 500*MeV) & (P > 3*GeV)", 
    "gamma" : "(PT > 300*MeV) & (P > 1.5*GeV)" 
  },
  CombinationCut = (
    "(ADAMASS('eta') < 150*MeV) & "
    "(AMAXDOCA('') < 0.4*mm) & "
    "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) & "
    "(AMAXCHILD('mu+' == ABSID, TRCHI2DOF) < 3) & "
    "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.4) & "
    "(AMINCHILD('mu+' == ABSID, PROBNNmu) > 0.4)"
  ),
  MotherCut = "ALL")
sel_comb = Selection(
  'SelComb',
  Algorithm = comb,
  RequiredSelections = [sel_dtrs])

# Apply mother cuts
cuts_mother = FilterDesktop('MotherCuts',
  # TODO: & (PT > 2000*MeV)
  Code = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")
sel_mother = Selection(
  "SelMother",
  Algorithm = cuts_mother,
  RequiredSelections=[sel_comb])

# Final selection sequence
seq = SelectionSequence('seqMuMuGamma', TopSelection = sel_mother)

# DaVinci algorithm sequence.
DaVinci().appendToMainSequence([seq])

# TisTos configuration.
from Configurables import ToolSvc, TriggerTisTos
for stage in ('Hlt1', 'Hlt2', 'Strip/Phys'):
    ToolSvc().addTool(TriggerTisTos, stage + "TriggerTisTos")
    tool = getattr(ToolSvc(), stage + "TriggerTisTos")
    tool.HltDecReportsLocation = '/Event/' + stage + '/DecReports'
    tool.HltSelReportsLocation = '/Event/' + stage + '/SelReports'

# GaudiPython configuration.
import GaudiPython, ROOT
from Configurables import ApplicationMgr 
gaudi = GaudiPython.AppMgr()
tes = gaudi.evtsvc() # Transient Event Storage

# Tools.
rndTool = ROOT.TRandom3(0)
genTool = gaudi.toolsvc().create(
    'DaVinciSmartAssociator',
    interface = 'IParticle2MCWeightedAssociator')
rftTool = gaudi.toolsvc().create(
    'PVOfflineTool',
    interface = 'IPVOfflineTool')
pvrTool = gaudi.toolsvc().create(
    'GenericParticle2PVRelator<_p2PVWithIPChi2, '
    'OfflineDistanceCalculatorName>/P2PVWithIPChi2',
    interface = 'IRelatedPVFinder')
dstTool = gaudi.toolsvc().create(
    'LoKi::TrgDistanceCalculator',
    interface = 'IDistanceCalculator')
trkTool = gaudi.toolsvc().create(
    'TrackMasterExtrapolator',
    interface = 'ITrackExtrapolator')
l0Tool = gaudi.toolsvc().create(
    'L0TriggerTisTos',
    interface = 'ITriggerTisTos')
hlt1Tool = gaudi.toolsvc().create(
    'TriggerTisTos/Hlt1TriggerTisTos',
    interface = 'ITriggerTisTos')
hlt2Tool = gaudi.toolsvc().create(
    'TriggerTisTos/Hlt2TriggerTisTos',
    interface = 'ITriggerTisTos')
physTool = gaudi.toolsvc().create(
    'TriggerTisTos/Strip/PhysTriggerTisTos',
    interface = 'ITriggerTisTos')
docaTool = GaudiPython.gbl.LoKi.Particles.DOCA(0, 0, dstTool)

# Daughter count overwhelms the graph and makes it unreadable. 
# Just note that about ~65 candidates pass through each event.
dcount, ccount, mcount, fcount = [], [], [], []
mc_ccount, mc_mcount, mc_fcount = [], [], []

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

# Initialize the tuple.
from Ntuple import Ntuple
ntuple = Ntuple('output.root', tes, genTool, rftTool, pvrTool, None, dstTool,
                None, trkTool, l0Tool, hlt1Tool, hlt2Tool)

# Run.
import sys, ROOT
try: evtmax = int(sys.argv[1])
except: evtmax = float('inf')
evtnum = 0

while evtnum < evtmax:
  gaudi.run(1) # Advance Gaudi by one event
  if not bool(tes['/Event']): break # Exit if no data found in TES
  evtnum += 1
  ntuple.clear()

  # Begin counters
  count('Phys/SelDaughters/Particles', dcount)
  count('Phys/SelComb/Particles', ccount)
  count('Phys/SelMother/Particles', mcount)
  count(seq.outputLocation(), fcount)

  count_mc_matched('Phys/SelComb/Particles', mc_ccount)
  count_mc_matched('Phys/SelMother/Particles', mc_mcount)
  count_mc_matched(seq.outputLocation(), mc_fcount)

  # Fill event info.
  daq = tes['DAQ/ODIN'] 
  try:
    ntuple.ntuple['run_n'][0] = daq.runNumber()
    ntuple.ntuple['evt_n'][0] = daq.eventNumber()
    ntuple.ntuple['evt_tck'][0] = daq.triggerConfigurationKey()
  except: continue
  # Save number of primary vertices
  try: ntuple.ntuple['pvr_n'][0] = len(tes['Rec/Vertex/Primary'])
  except: continue
  try: ntuple.ntuple['evt_spd'][0] = GaudiPython.gbl.LoKi.L0.DataValue(
    'Spd(Mult)')(tes['Trig/L0/L0DUReport'])
  except: pass

  # Create tools.
  if not ntuple.detTool: ntuple.detTool = gaudi.detSvc()[
    '/dd/Structure/LHCb/BeforeMagnetRegion/Velo']

  # Fill candidates.
  fill = False

  pvrs = tes['Rec/Vertex/Primary']
  trks = tes['Rec/Track/Best']
  prts = tes[seq.outputLocation()]
  sigs = []
  try: len(prts); run = True
  except: run = False
  if run:
    for prt in prts:
      if (prt.measuredMass() < 1500 or 'CL' not in Type
        or rndTool.Rndm() < 0.1):
        sigs += [prt]; ntuple.fillPrt(prt, pvrs, trks); fill = True

  # MC-only muons for truth matching.
  # prts = tes['Phys/StdAllLooseMuons/Particles'] # Save all loose selected muons
  prts = tes['Phys/StdAllLooseMuons/Particles'] # Save loose selected muons
  try: len(prts); run = True
  except: run = False
  if run:
    for prt in prts: ntuple.fillPrt(prt); fill = True
  
  # MC-only photons for truth matching
  prts = tes['Phys/StdLoosePhotons/Particles']
  try: len(prts); run = True
  except: run = False
  if run:
    for prt in prts: ntuple.fillPrt(prt); fill = True

  # Fill MC.
  mcps = tes['MC/Particles']
  try: len(mcps); run = True
  except: run = False
  if run:
    # Loop over MC truth particles
    for mcp in mcps:
      pid = mcp.particleID()
      # Save only muons and etas
      if abs(pid.pid()) in [13, 221]: # 13 = mu, 221 = eta
        ntuple.fillMcp(mcp)
        fill = True

  # Fill ntuple if there is information for this event
  if fill: ntuple.fill()

# Close and write the output.
ntuple.close()

plot_counters("candidate_survival", [[], ccount, mcount, fcount])
plot_counters("mc_match", [[], mc_ccount, mc_mcount, mc_fcount])
