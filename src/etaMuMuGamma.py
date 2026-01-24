#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu (gamma) analysis                        #
# Author: Michael Peters                                                      #
###############################################################################

# DaVinci configuration.
from Configurables import DaVinci
from GaudiConf import IOHelper

# Analyisis Production configuration.
try: from ProdConf import ProdConf
except ImportError: pass
from pathlib import Path
import argparse
import sys
import importlib.machinery
import importlib.util

# =============================================================================

def parseArgs() -> bool:
    """
    Parser method required for GaudiPython to work with AnalysisProductions

    Argument parser method required for GaudiPython to work with
    AnalysisProductions because the option files get unsorted.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "options_paths",
        nargs="+",
        type=Path,
        help="Additional options files to load before starting Gaudi Python",
    )
    args = parser.parse_args()

    opt = None
    for options_path in args.options_paths:
        print("Adding options file:", options_path)
        loader = importlib.machinery.SourceFileLoader(
            "option", str(options_path.resolve()))
        spec = importlib.util.spec_from_loader(loader.name, loader)
        mod = importlib.util.module_from_spec(spec)
        loader.exec_module(mod)
        try: opt = mod.backwards
        except AttributeError: pass

    if opt is None:
        opt = 'False'
        print("Warning: No options file set BDTTagger.Backwards." +
              "Defaulting to 'False'.")
    else: print(f"Running with BDTTagger.Backwards == '{opt}'.")

    if opt == 'True': return True
    elif opt == 'False': return False
    else: raise ValueError("Invalid input, must be 'True' or 'False'")


# =============================================================================

# Set flags
IS_MC = True  # True = MC, False = real data
IS_SAMPLE = False # True = sample data, False = analysis production
IS_MUMUGAMMA = True  # True = η→μμγ, False = η→μμ

# MC or real data.
if IS_MC and IS_SAMPLE:
    DaVinci().DataType = '2018'
    DaVinci().Lumi = False  # Processing of luminosity data.
    DaVinci().Simulation = True  # MC simulation data.
    # Found using "lb-dirac dirac-bookkeeping-production-information 00169948".
    # DaVinci().DDDBtag = 'dddb-20210528-8' # for 00169948
    # DaVinci().CondDBtag = 'sim-20201113-8-vc-md100-Sim10' # for 00169948
    DaVinci().DDDBtag = 'dddb-20170721-3'  # for 00090844
    DaVinci().CondDBtag = 'sim-20190128-vc-md100'  # for 00090844
    IOHelper('ROOT').inputFiles([
        'data/minbias/00090844_00000001_7.AllStreams.dst',  # minbias
        # 'data/minbias/00090844_00000048_7.AllStreams.dst',
        # 'data/minbias/00090844_00000055_7.AllStreams.dst',
        # 'data/minbias/00090844_00000075_7.AllStreams.dst',
        # 'data/minbias/00090844_00000079_7.AllStreams.dst',
        # 'data/minbias/00090844_00000108_7.AllStreams.dst',
        # 'data/minbias/00090844_00000186_7.AllStreams.dst',
        # 'data/minbias/00090844_00000193_7.AllStreams.dst',
        # 'data/minbias/00090844_00000207_7.AllStreams.dst',
        # 'data/minbias/00090844_00000227_7.AllStreams.dst',
        # 'data/minbias/00090844_00000054_7.AllStreams.dst',
        # 'data/minbias/00090844_00000176_7.AllStreams.dst',
        # 'data/norm/00169948_00000003_7.AllStreams.dst',  # eta->mumugamma
        # 'data/norm/00169948_00000138_7.AllStreams.dst'  # 39112231, sim10b, magdown
    ],
        clear=True)
    
    # For logging date and time
    from datetime import datetime

    # Get current date and time, append .root file extension
    extension = "_" + str(datetime.now().strftime("%Y%m%d")) + ".root"
elif IS_MC and not IS_SAMPLE:
    DaVinci().Lumi = False  # Processing of luminosity data.
    # Output file
    backwards = parseArgs()
    outfile = f"{ProdConf().OutputFilePrefix}.{ProdConf().OutputFileTypes[0]}"

# Reconstruction.
from Configurables import CombineParticles
from StandardParticles import StdLooseMuons as muons
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

# Decay mode config
daughter_cuts = {
    "mu+": "(PT > 500*MeV) & (P > 3*GeV)",
    "mu-": "(PT > 500*MeV) & (P > 3*GeV)"
}
required_selections = [muons]
if IS_MUMUGAMMA:
    # Append gamma cuts and selection
    daughter_cuts["gamma"] = "(PT > 500*MeV) & (CL > 0.2)"
    required_selections.append(photons)
    if IS_SAMPLE: outfile = 'ntuples/eta2MuMuGamma' + ('_mc' if IS_MC else '') + extension
    decay_descriptor = "eta -> mu+ mu- gamma"
else:
    if IS_SAMPLE: outfile = 'ntuples/eta2MuMu' + ('_mc' if IS_MC else '') + extension
    decay_descriptor = "eta -> mu+ mu-"

print(f"Writing output to {outfile}")  # debug

# Combination cuts
combination_cuts = (
    "(ADAMASS('eta') < 150*MeV) & "  # change based on side bands
    "(AMAXDOCA('') < 0.4*mm) & "  # doca btwn children
    # possibly change TRCHI2DOF to 2.5
    "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) & "  # track
    "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.95)"  # muon weights
)

# Apply cuts
comb = CombineParticles(
    'combEtaMuMuGamma',
    DecayDescriptor=decay_descriptor,
    DaughtersCuts=daughter_cuts,
    CombinationCut=combination_cuts,
    # vertex fit can fail, must pass
    # vertex, hard to model this cut eff from data
    MotherCut="(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")

# Selection
sel_comb = Selection(
    'selEtaMuMuGamma',
    Algorithm=comb,
    RequiredSelections=required_selections)

# Final selection sequence
seq = SelectionSequence('seqMuMuGamma', TopSelection=sel_comb)

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
import GaudiPython
import ROOT
gaudi = GaudiPython.AppMgr()
tes = gaudi.evtsvc()  # Transient Event Storage

# Tools.
rndTool = ROOT.TRandom3(0)
genTool = gaudi.toolsvc().create(
    'DaVinciSmartAssociator',
    interface='IParticle2MCWeightedAssociator')
rftTool = gaudi.toolsvc().create(
    'PVOfflineTool',
    interface='IPVOfflineTool')
pvrTool = gaudi.toolsvc().create(
    'GenericParticle2PVRelator<_p2PVWithIPChi2, '
    'OfflineDistanceCalculatorName>/P2PVWithIPChi2',
    interface='IRelatedPVFinder')
dstTool = gaudi.toolsvc().create(
    'LoKi::TrgDistanceCalculator',
    interface='IDistanceCalculator')
trkTool = gaudi.toolsvc().create(
    'TrackMasterExtrapolator',
    interface='ITrackExtrapolator')
l0Tool = gaudi.toolsvc().create(
    'L0TriggerTisTos',
    interface='ITriggerTisTos')
hlt1Tool = gaudi.toolsvc().create(
    'TriggerTisTos/Hlt1TriggerTisTos',
    interface='ITriggerTisTos')
hlt2Tool = gaudi.toolsvc().create(
    'TriggerTisTos/Hlt2TriggerTisTos',
    interface='ITriggerTisTos')
physTool = gaudi.toolsvc().create(
    'TriggerTisTos/Strip/PhysTriggerTisTos',
    interface='ITriggerTisTos')
docaTool = GaudiPython.gbl.LoKi.Particles.DOCA(0, 0, dstTool)

# Initialize the tuple.
# local sample
try: from scripts.Ntuple import Ntuple
# analysis production
except: from Ntuple import Ntuple
ntuple = Ntuple(outfile, IS_MC, tes, genTool, rftTool, pvrTool,
                None, dstTool, None, trkTool, l0Tool, hlt1Tool, hlt2Tool)

# Run.
# local sample
try: evtmax = args.evtmax if args.evtmax > 0 else float("inf")
except: 
    # analysis production
    try: evtmax = int(sys.argv[1])
    except: evtmax = float("inf")
evtnum = 0

while evtnum < evtmax:
    gaudi.run(1)  # Advance Gaudi by one event
    if not bool(tes['/Event']): break  # Exit if no data found in TES
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

    # Fill MC.
    if IS_MC:
        mcps = tes['MC/Particles']
        try: len(mcps); run = True
        except: run = False
        if run:
            for mcp in mcps:
                if abs(mcp.particleID().pid()) == 221:
                    ntuple.fillMcp(mcp)
                    fill = True

    pvrs = tes['Rec/Vertex/Primary']
    trks = tes['Rec/Track/Best']
    prts = tes[seq.outputLocation()]
    sigs = []
    try: len(prts); run = True
    except: run = False
    if run:
        for prt in prts:
            sigs += [prt]
            ntuple.fillPrt(prt, pvrs, trks)
            fill = True

    # Fill ntuple if there is information for this event
    if fill: ntuple.fill()  # Debug

# Close and write the output.
ntuple.close()
