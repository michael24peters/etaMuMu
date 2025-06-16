#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu analysis                                #
# Author: Michael Peters                                                      #
###############################################################################

# Data.
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['00169948_00000003_7.AllStreams.dst',
                            '00169948_00000138_7.AllStreams.dst',
                            '00169948_00000186_7.AllStreams.dst',
                            '00169948_00000319_7.AllStreams.dst',
                            '00169948_00000411_7.AllStreams.dst'
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
outfile = 'etaMuMuGamma_v3.root'

### ?? RELATED TO TURBO WHICH I DID NOT INCLUDE HERE ??
### ?? SHOULD I HAVE BOTH THE TRIGGER AND MC TYPES ??
# Tag configuration.
TrkCats = [('ve', 1), ('tt', 2), ('it', 3), ('ot', 4), ('mu', 7)]
TrgLocs = ['Hlt2ExoticaPrmptDiMuonSSTurbo', 
           'Hlt2ExoticaPrmptDiMuonTurbo',
           'Hlt2ExoticaDiMuonNoIPTurbo',
           'Hlt2ExoticaDisplDiMuon']

# Reconstruction.
from Configurables import CombineParticles
from StandardParticles import StdLooseMuons as muons
# TODO: only loose photon cut applied
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

# Create the di-muon + gammas.
comb = CombineParticles( # Combine daughters to reconstruct composite particle.
  'combMuMuGamma',
  DecayDescriptor = "eta -> mu+ mu- gamma",
  CombinationCut = (
    # Requires absolute difference between invariant mass of combined 
    # daughters and PDG eta mass (~547 MeV) must be within 100 MeV.
    ## TODO: Test 30*MeV
    "(ADAMASS('eta') < 100*MeV) & "
    # Requires distance of closest approach within 0.2 to ensure daughters 
    # close enough to likely come from common decay vertex.
    "(AMAXDOCA('') < 0.2*mm) & "
    "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) & "
    "(AMAXCHILD('mu+' == ABSID, TRCHI2DOF) < 3) & "
    # Requires muon PID probability greater than 0.5.
    "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.5) & "
    "(AMINCHILD('mu+' == ABSID, PROBNNmu) > 0.5)"
  ),
  # Higher end thresholds chosen to favor boosted eta (since eta is itself 
  # boosted from parent decay, e.g. B -> eta X). The cost of higher thresholds 
  # is signal loss.
  ## TODO: plot P and PT distributions to determine the best values.
  DaughtersCuts = {
    # Muon thresholds P > 3-5 GeV, PT > 300-600 MeV
    ## TODO: Add TRGHOSTPROB < 0.4? Muon ghost probability, i.e. likely fake 
    ## tracks.
    ## !! Do I need ISMUON here? or is it redundant?
    "mu+" : "(PT > 500*MeV) & (P > 10*GeV) & (ISMUON)", 
    "mu-" : "(PT > 500*MeV) & (P > 10*GeV) & (ISMUON)", 
    # Photon reconstruction (ECAL) thresholds P > 2-5 GeV, PT > 300-500 MeV
    ## TODO: Add CL > 0.2? ECAL photon confidence level
    "gamma" : "(PT > 650*MeV) & (P > 5*GeV)" 
  },
  # VFASPF == Vertex Fit Algorithm from the Standard Particle Finder, i.e. 
  # vertex fitting information.
  # VCHI2PDOF == Vertex chi^2 per degree of freedom for the fitted vertex, i.e.
  # how well track fitting matches actual detector hits. Lower equals better
  # agreement. Typically 5 DOF (x, y, z, px/p, py/p, curvature).
  ## TODO: consider tightening to 5 (and as low as 3).
  MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)"
)
sel = Selection(
  'SelMuMuGamma',
  Algorithm = comb,
  RequiredSelections = [muons, photons])
seq = SelectionSequence('seqMuMuGamma', TopSelection = sel)

# DaVinci algorithm sequence.
DaVinci().appendToMainSequence([seq])

# Create the jets.
### ?? DO I NEED TO DO THIS ??
# Create the jets.
from Configurables import HltParticleFlow, HltJetBuilder
pf = HltParticleFlow('pf')
pf.Inputs += [
['ProtoParticle',  'best',  'Rec/ProtoP/Charged'],
['ProtoParticle',  'gamma', 'Rec/ProtoP/Neutrals']]
pf.Inputs += [['Particle', 'daughters', seq.outputLocation()]]
pf.Output = 'Phys/PF/Particles'
pf.ProBestNames = ['mu+']
pf.ProBestKeys  = [701]
pf.ProBestMins  = [0.5]
pf.EcalBest = True
pf.SprRecover = False
pf.TrkLnErrMax = 10
pf.TrkUpErrMax = 10
pf.TrkDnErrMax = 10
jb = HltJetBuilder('jb')
jb.JetEcPath = ''
jb.Inputs = [pf.Output]
jb.Output = 'Phys/JB/Particles'
jb.JetPtMin = 0
DaVinci().appendToMainSequence([pf, jb])

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

  # Fill event info.
  ### ?? WHAT IS THE PURPOSE OF THIS LINE ??
  daq = tes['DAQ/ODIN'] 
  try:
    ntuple.ntuple['run_n'][0] = daq.runNumber()
    ntuple.ntuple['evt_n'][0] = daq.eventNumber()
    ntuple.ntuple['evt_tck'][0] = daq.triggerConfigurationKey()
  except: continue
  # Save number of primary vertices
  try: ntuple.ntuple['pvr_n'][0] = len(tes['Rec/Vertex/Primary'])
  except: continue
  ### ?? WHAT IS THE PURPOSE OF THIS LINE ??
  try: ntuple.ntuple['evt_spd'][0] = GaudiPython.gbl.LoKi.L0.DataValue(
    'Spd(Mult)')(tes['Trig/L0/L0DUReport'])
  except: pass

  # Create tools.
  ### ?? WHAT IS THE PURPOSE OF THIS LINE ??
  if not ntuple.detTool: ntuple.detTool = gaudi.detSvc()[
    '/dd/Structure/LHCb/BeforeMagnetRegion/Velo']

  # Fill candidates.
  fill = False # Track if successfully filled

  ### ?? REMOVED SAME-SIGN DI-MUONS SELECTIONS (seqSsmus, DiMuon.py line 88) ??
  # prts = tes[seq.outputLocation()]
  # try: len(prts); run = True
  # except: run = False
  # if run:
  #   for prt in prts: ntuple.fillPrt(prt, pvrs, trks); fill = True

  pvrs = tes['Rec/Vertex/Primary']
  trks = tes['Rec/Track/Best']
  prts = tes[seq.outputLocation()]
  sigs = []
  try: len(prts); run = True
  except: run = False
  if run:
    for prt in prts:
      ### ?? I DONT UNDERSTAND THIS IF STATEMENT ??
      ### Attempt to answer: filters out high mass decays not relevant
      ### CL = control line, unsure of purpose
      ### Take 10% of events at random, unsure of purpose
      if (prt.measuredMass() < 1500 or 'CL' not in Type
        or rndTool.Rndm() < 0.1):
        sigs += [prt]; ntuple.fillPrt(prt, pvrs, trks); fill = True
  ### ?? WHAT IS THE PURPOSE OF THIS LINE ??
  # pf = Particle Flow
  prts = tes[pf.Output]
  try: len(prts); run = len(sigs) > 0 # Need at least one signal candidate.
  except: run = False
  if run:
    for prt in prts:
      ### ?? If this fails, what happens? How does continue work?
      if abs(prt.charge()) != 1: continue # Filter out neutral particles
      # Try to associate PF particle with PV
      try: pvr = pvrTool.relatedPV(prt, 'Rec/Vertex/Primary')
      except: pvr = None
      if not pvr: continue
      ### ?? WAS ROOT.Double(-1) ??
      from ctypes import c_double
      ip, ipChi2 = c_double(-1.0), c_double(-1.0)
      dstTool.distance(prt, pvr, ip, ipChi2)
      ### ?? CHECK MY UNDERSTANDING OF THIS SECTION ??
      # Require large impact parameter: track of shortest distance between 
      # particle's trajectory and PV. In this case this is the measure of how 
      # inconsistent the track is with originating from the PV. If it comes
      # from the PV, it shoul dhave a small IP and ipChi2; if from a SV 
      # (decay), it should have a larger IP and ipChi2. ipChi2 = 9 is at 3
      # sigma level.
      # Pt() > 200 is a soft cut to exclude very low energy clutter.
      # Need to convert c_double to python float via .value
      if ipChi2.value > 9 and prt.momentum().Pt() > 200:
        # Check if PF candidate is near a daughter of signal candidate
        minDoca = 100 # Placeholder; 100mm is a very large value
        # Loop through PF candidates
        for sig in sigs:
          # Loop through PF candidate daughters
          for dtr in sig.daughters():
            doca = docaTool.doca(prt, dtr)
            if doca < minDoca: minDoca = doca
          # If DOCA is small to daughter of known signal, muon might be near
          # the decay and is thus worth filling. In other words, it's checking
          # whether there are good quality, displaced muons near signal decay.
          if minDoca < 0.2: ntuple.fillPrt(prt); fill = True

  # MC-only muons for truth matching.
  prts = tes['Phys/StdAllLooseMuons/Particles'] # Save all loose selected muons
  try: len(prts); run = True
  except: run = False
  if run:
    for prt in prts: ntuple.fillPrt(prt); fill = True

  # Fill MC.
  ### ?? WHAT IS THE PURPOSE OF THESE LINES ??
  mcps = tes['MC/Particles']
  try: len(mcps); run = True
  except: run = False
  if run:
    # Loop over MC truth particles
    for mcp in mcps:
      pid = mcp.particleID()
      ### ?? DO I NEED TO FILTER ALL OF THESE PARTICLES ??
      # if (abs(pid.pid()) in [13, 23, 32, 36] or pid.hasCharm() or 
      #   pid.hasBottom()): ntuple.fillMcp(mcp); fill = True
      # Save only muons and etas
      if abs(pid.pid()) in [13, 221]: # 13 = mu, 221 = eta
        ntuple.fillMcp(mcp)
        fill = True

  # Fill ntuple if there is information for this event
  if fill: ntuple.fill()

# Close and write the output.
ntuple.close()