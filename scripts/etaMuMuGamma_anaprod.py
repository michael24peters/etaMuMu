#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu analysis                                #
# Author: Michael Peters                                                      #
###############################################################################

# DaVinci configuration.
from Configurables import DaVinci
from GaudiConf import IOHelper

from Configurables import CombineParticles, FilterDesktop
from StandardParticles import StdLooseMuons as muons
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

# Analysis Production configuration
from ProdConf import ProdConf

def parseArgs() -> bool:
  """ 
  Parser method required for GaudiPython to work with AnalysisProductions

  Argument parser method required for GaudiPython to work with
  AnalysisProductions because the option files get unsorted.
  """

  parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    description = __doc__,
  )
  parser.add_argument(
    "options_paths",
    nargs = "+",
    type = Path,
    help = "Additional options files to load before starting Gaudi Python",
  )
  args = parser.parse_args()

  opt = None
  for options_path in args.options_paths:
    print("Adding options file:", options_path)
    loader = importlib.machinery.SourceFileLoader(
      "option", str(options_path.resolve())
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    try: opt = mod.backwards
    except AttributeError: pass

  if opt is None:
    raise parser.error("One options file must set the BDTTagger mode")
  else:
    print(f"Running with BDTTagger.Backwards == '{opt}'.")
  
  if   opt == 'True' : return True
  elif opt == 'False': return False
  else: raise ValueError("Invalid input, must be 'True' or 'False'")

DaVinci().Lumi = False # Processing of luminosity data.

# Output file
backwards = parseArgs()
outfile = f"{ProdConf().OutputFilePrefix}.{ProdConf().OutputFileTypes[0]}"

# Decay descriptor
decay_descriptor = "eta -> mu+ mu- gamma"
# Daughter cuts
daughter_cuts = {
  "mu+" : "(PT > 500*MeV) & (P > 3*GeV)",
  "mu-" : "(PT > 500*MeV) & (P > 3*GeV)",
  # ?? PT > 300*MeV & P>1.5*GeV & PROBNNgamma > 0>
  "gamma" : "(PT > 500*MeV) & (CL > 0.2)" 
}
required_selections = [muons, photons]

# Combination cuts
combination_cuts = (
  "(ADAMASS('eta') < 150*MeV) & " # change based on side bands
  "(AMAXDOCA('') < 0.4*mm) & " # doca btwn children
  "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) & " # track
  "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.4)" # muon weights
)
print(f"outfile name: {outfile}")

# Apply cuts
comb = CombineParticles(
  'combEtaMuMuGamma',
  DecayDescriptor = decay_descriptor,
  DaughtersCuts = daughter_cuts,
  CombinationCut = combination_cuts,
  # vertex fit can fail, must pass
  # vertex, hard to model this cut eff from data
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
ntuple = Ntuple(outfile, tes, genTool, rftTool, pvrTool, None, dstTool,
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
        # Look at every eta, fill eta -> mu+ mu- gamma decay instances
        if abs(mcp.particleID().pid()) == 221:
          ntuple.fillGen(mcp)
          fill = True

  # Fill ntuple if there is information for this event
  if fill: ntuple.fill() # Debug

# Close and write the output.
ntuple.close()
