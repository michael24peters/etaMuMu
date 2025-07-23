#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu analysis                                #
# Author: Michael Peters                                                      #
###############################################################################

## TODO: ind constituents, muon pid (cuts) and photon, pt photon, 
## TODO: linking: take recon. obj, what hits prod rec items, can it be linked back to gen particle, i.e. did muon come from eta?
#### truth matching
## TODO: # cand per event, possible getting photons from elsewhere (not eta), what is faking signal, what % from eta or other things (what are other things)
## TODO: talk to christian abt the dec file

## use daughters to see number of photons reconstructed; use pid to see number of photons in MC; compare difference

# Data.
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['data/00169948_00000003_7.AllStreams.dst',
                            # 'data/00169948_00000138_7.AllStreams.dst',
                            # 'data/00169948_00000186_7.AllStreams.dst',
                            # 'data/00169948_00000319_7.AllStreams.dst',
                            # 'data/00169948_00000411_7.AllStreams.dst'
                            ],
                           clear = True)

# Data type configuration.
from GaudiKernel import SystemOfUnits as Units
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = True
DaVinci().DataType = '2018'
# DDDBtag and CONDBtag found using: https://lhcb.github.io/starterkit-lessons/first-analysis-steps/minimal-dv-job.html
DaVinci().DDDBtag = 'dddb-20210528-8'
DaVinci().CondDBtag = 'sim-20201113-8-vc-md100-Sim10'

# Ouput file
outfile = 'etaMuMuGamma.root'

# Onia reconstruction
from Configurables import CombineParticles
from StandardParticles import StdLooseMuons as muons
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

comb = CombineParticles('etaComb') # acting on reconstructed data (tracks, calo, rich, etc.)
comb.DecayDescriptor = '[eta -> mu+ mu- gamma]cc'
comb.DaughtersCuts = {'mu+' : 'ISMUON', 'mu-' : 'ISMUON'}
comb.CombinationCut = "ADAMASS('eta') < 100*MeV"
comb.MotherCut = '(VFASPF(VCHI2PDOF) < 20) & (PT > 2.0*GeV)'
selection = Selection('SelEta', Algorithm = comb, 
                    RequiredSelections = [muons, photons])
etaSeq = SelectionSequence('SeqEta', TopSelection = selection)

# DaVinci algorithm sequence
DaVinci().appendToMainSequence([etaSeq])

# GaudiPython configuration.
import GaudiPython
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc() # transient event storage

# Initialize the tuple
from Ntuple import Ntuple
ntuple = Ntuple(outfile, tes, gaudi.toolsvc(), gaudi.detSvc())
mom = ['px', 'py', 'pz', 'e']
pos = ['x', 'y', 'z']
cov = ['dx', 'dy', 'dz', 'chi2', 'ndof']
ntuple.init('gen', ['pid'] + mom) # gen == generated (mc) particles ; (prefix, array of information you want to store)
ntuple.init('gen_pvr', pos + cov) # pvr == primary vertex
ntuple.init('rec_pvr', pos + cov) # rec == reconstructed
ntuple.init('can', ['m'] + mom) # can == eta (candidates ; m==mass, t==lifetime
ntuple.init()

# Run
import sys
try: evtmax = int(sys.argv[1]) # number of events to run over (from terminal)
except: evtmax = float('inf')
evtnum = 0
while evtnum < evtmax:
  gaudi.run(1)
  if not bool(tes['/Event']): break
  evtnum += 1
  ntuple.clear()
  fill = False

  # mcpart.endVertices()[0].products() # vector of MCParticle daughters

  # Fill generated info
  try:
    gens = tes['MC/Particles']
    for gen in gens: # Each gen is of type MCParticle
      vrs = {}
      # Fill eta, eta', mu+, mu-
      if abs(gen.particleID().pid()) in [221, 331, 13]:
        try:
          # Fill pid
          vrs['pid'] = gen.particleID().pid()
          # Fill momentum
          ntuple.fillMom(gen.momentum(), vrs)
          ntuple.fill('gen', vrs = vrs)
        except Exception as e:
          print(f"Error while processing generated pid and mom: {e}")
          pass
        # Fill the PVs for eta, eta'
        try:
          if gen.particleID().pid() == 221 or gen.particleID().pid() == 331:
            ntuple.fillPvr(gen.primaryVertex(), vrs)
        except Exception as e:
          print(f"Error while processing generated PV info: {e}")
          pass
        fill = True
  except Exception as e:
    print(f"Error while processing generated info: {e}")
    pass
  
  # Fill reconstructed info
  try:
    etas = tes[etaSeq.outputLocation()]
    # Loop over the eta candidates
    for cand in etas: # Each cand is of type Particle (or DataObject)
      vrs = {}
      # Fill m
      vrs['m'] = cand.measuredMass()
      # Fill mom
      ntuple.fillMom(cand.momentum(), vrs)
      ntuple.fill('can', vrs = vrs)
      # Fill the PVs
      ntuple.fillPvr(cand, vrs, pre='rec_pvr')
  except: pass

  # Fill the ntuple
  if fill: ntuple.fill()
# Close the ntuple
ntuple.close()
