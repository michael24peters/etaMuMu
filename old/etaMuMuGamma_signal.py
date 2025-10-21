#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu analysis                                #
# Author: Michael Peters                                                      #
###############################################################################

IS_MC = True
IS_MUMUGAMMA = True

# DaVinci configuration.
from Configurables import DaVinci
from GaudiConf import IOHelper

import argparse
# Command line arguments
parser = argparse.ArgumentParser(description="LHCb decay mode and data selection")
parser.add_argument("-m", "--mc",
                    action="store_true",
                    help="Run on MC / simulation")

parser.add_argument("-g", "--mumugamma",
                    action="store_true",
                    help="Analyze η→μμγ (default: η→μμ)")

parser.add_argument("--evtmax", type=int, default=-1,
                    help="Max events to process (-1 = all)")

args = parser.parse_args()

# Set flags
IS_MC = args.mc                  # True = simulation, False = real data
IS_MUMUGAMMA = args.mumugamma    # True = η→μμγ, False = η→μμ

# MC or real data.
if IS_MC:
  DaVinci().DataType = '2018'
  DaVinci().Lumi = False # Processing of luminosity data.
  DaVinci().Simulation = True # MC simulation data.
  # Found using "lb-dirac dirac-bookkeeping-production-information 00169948".
  DaVinci().DDDBtag = 'dddb-20210528-8'
  DaVinci().CondDBtag = 'sim-20201113-8-vc-md100-Sim10'
  IOHelper('ROOT').inputFiles([
                               'data/00169948_00000003_7.AllStreams.dst',
                               'data/00169948_00000138_7.AllStreams.dst',
                              ],
                              clear = True)
else:
  DaVinci().DataType = '2018' # !! TODO
  DaVinci().Lumi = True
  DaVinci().Simulation = False
  # TODO: IOHelper input pointing to data DSTs from bookkeeping

# Reconstruction.
from Configurables import CombineParticles, FilterDesktop
from StandardParticles import StdLooseMuons as muons
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

# Decay mode config
if IS_MUMUGAMMA == True:
  # Output file
  outfile = 'eta2MuMuGamma' + ('_mc' if IS_MC else '_data') + '.root'
  # Decay descriptor
  decay_descriptor = "eta -> mu+ mu- gamma"
  # Daughter cuts
  daughter_cuts = {
    "mu+" : "(PT > 500*MeV) & (P > 3*GeV)",
    "mu-" : "(PT > 500*MeV) & (P > 3*GeV)",
    "gamma" : "(PT > 500*MeV) & (CL > 0.2)" # ?? PT > 300*MeV & P>1.5*GeV & PROBNNgamma > 0>
  }
  required_selections = [muons, photons]
# Repeat for eta -> mu+ mu-
else:
  outfile = 'eta2MuMu' + ('_mc' if IS_MC else '_data') + '.root'
  decay_descriptor = "eta -> mu+ mu-"
  daughter_cuts = {
    "mu+" : "(PT > 500*MeV) & (P > 3*GeV)",
    "mu-" : "(PT > 500*MeV) & (P > 3*GeV)",
  }
  required_selections = [muons]

# Combination cuts
combination_cuts = (
  "(ADAMASS('eta') < 150*MeV) & "
  "(AMAXDOCA('') < 0.4*mm) & "
  "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) & "
  "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.4)"
)
print(f"outfile name: {outfile}")

# Apply cuts
comb = CombineParticles(
  'combEtaMuMuGamma',
  DecayDescriptor = decay_descriptor,
  DaughtersCuts = daughter_cuts,
  CombinationCut = combination_cuts,
  MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")

# Selection
sel_comb = Selection(
  'selEtaMuMuGamma',
  Algorithm = comb,
  RequiredSelections = required_selections)

# Final selection sequence
seq = SelectionSequence('seqMuMuGamma', TopSelection = sel_comb)

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

# Initialize the tuple.
from scripts.Ntuple import Ntuple
ntuple = Ntuple(outfile, IS_MC, IS_MUMUGAMMA, tes, genTool, rftTool, pvrTool, None, dstTool,
                None, trkTool, l0Tool, hlt1Tool, hlt2Tool)

# Run.
import sys, ROOT
evtmax = args.evtmax if args.evtmax > 0 else float("inf")
evtnum = 0

while evtnum < evtmax:
  gaudi.run(1) # Advance Gaudi by one event
  if not bool(tes['/Event']): break # Exit if no data found in TES
  evtnum += 1
  ntuple.clear()

  # Fill event info.
  daq = tes['DAQ/ODIN'] 
  try:
    ntuple.ntuple['run_n'][0] = daq.runNumber()
    ntuple.ntuple['evt_n'][0] = daq.eventNumber()
    ntuple.ntuple['evt_tck'][0] = daq.triggerConfigurationKey()
  except: continue
  # Save number of primary vertices
  try: ntuple.ntuple['pvr_n'][0] = len(tes['Rec/Vertex/Primary'])
  except: pass
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
      sigs += [prt]; ntuple.fillPrt(prt, pvrs, trks); fill = True

  # Fill MC.
  if IS_MC:
    mcps = tes['MC/Particles']
    try: len(mcps); run = True
    except: run = False
    if run:
      # Loop over MC truth particles
      for mcp in mcps:
        # Look at every eta
        if abs(mcp.particleID().pid()) == 221:
          ntuple.fillGen(mcp)
          fill = True

  # Fill ntuple if there is information for this event
  if fill: ntuple.fill() # Debug

# Close and write the output.
ntuple.close()
