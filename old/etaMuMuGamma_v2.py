#!/bin/env python
###############################################################################
# GaudiPython script for eta -> mu mu analysis                                #
# Author: Michael Peters                                                      #
###############################################################################

'''
look for patterns in false reconstruction
eta -> pi pi gamma
break into categories

how many events run over
how many events after reconstr eta -> mu mu gamma
how many '' after truth matching
apply cuts

repeat on minbias mc
look for background

'''

# Data.
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['00169948_00000003_7.AllStreams.dst',
                            '00169948_00000138_7.AllStreams.dst',
                            '00169948_00000186_7.AllStreams.dst',
                            '00169948_00000319_7.AllStreams.dst',
                            '00169948_00000411_7.AllStreams.dst'
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

# Ouput file.
outfile = 'etaMuMuGamma.root'

# Onia reconstruction.
from Configurables import CombineParticles
from StandardParticles import StdLooseMuons as muons
from StandardParticles import StdLooseAllPhotons as photons
from PhysSelPython.Wrappers import Selection, SelectionSequence

# Below acts on reconstructed data (tracks, calo, rich, etc.)
comb = CombineParticles('etaComb') 
comb.DecayDescriptor = '[eta -> mu+ mu- gamma]cc'
comb.DaughtersCuts = {'mu+' : 'ISMUON', 'mu-' : 'ISMUON',
  'mu+' : '(PT > 500*MeV) & (P > 10000*MeV)', 'mu-' : '(PT > 500*MeV) & (P > 10000*MeV)',
  'gamma' : '(PT > 650*MeV) & (P > 5000*MeV)'
}
# TODO: 2 < n < 4.5
comb.CombinationCut = "ADAMASS('eta') < 100*MeV"
comb.MotherCut = '(VFASPF(VCHI2PDOF) < 20) & (PT > 2.0*GeV)'
selection = Selection('SelEta', Algorithm = comb, 
                    RequiredSelections = [muons, photons])
etaSeq = SelectionSequence('SeqEta', TopSelection = selection)

# DaVinci algorithm sequence.
DaVinci().appendToMainSequence([etaSeq])

# GaudiPython configuration.
import GaudiPython
gaudi = GaudiPython.AppMgr()
tes = gaudi.evtsvc() # transient event storage

# TODO: move truth matching to Ntuple.py
# Create instance of the DaVinci truth matching tool
mctool = gaudi.toolsvc().create(
    'DaVinciSmartAssociator',
    interface = 'IParticle2MCWeightedAssociator')

# Initialize the tuple
from Ntuple import Ntuple
ntuple = Ntuple(outfile, tes, gaudi.toolsvc(), gaudi.detSvc())
mom = ['px', 'py', 'pz', 'e']
pos = ['x', 'y', 'z']
cov = ['dx', 'dy', 'dz', 'chi2', 'ndof']
ntuple.init('gen', ['pid'] + mom) # generated (mc) particles
ntuple.init('gen_pvr', pos + cov) # generated (mc) eta primary vertices
ntuple.init('gen_dtr', ['pid'] + mom) # generated (mc) daughters
ntuple.init('can', ['m'] + mom) # reconstructed eta candidates
ntuple.init('rec_pvr', pos + cov) # reconstructed eta primary vertices
ntuple.init('can_dtr', ['m'] + mom) # reconstructed daughter candidates

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

  # TODO: more specfic try-except (for now just passes if there is a failure)
  try:
    for eta in tes[etaSeq.outputLocation()]: #this would be inside event loop, looping through (reconstructed) candidates in TES location
      breakpoint() # DEBUG; breakpoint 1
      vrs = {}
      # Fill mass, momentum, and primary vertices (RECONSTRUCTED - eta)
      vrs['m'] = eta.measuredMass()
      ntuple.fillMom(eta.momentum(), vrs)
      ntuple.fill('can', vrs = vrs)
      ntuple.fillPvr(eta, vrs, pre='rec_pvr')
      fill = True

      # Loop through the daughter particles associated with the parent
      for dp in eta.daughtersVector(): # type <vector<Particle>>
        breakpoint() # DEBUG; breakpoint 2
        vrs = {}
        # Fill mass and momentum (RECONSTRUCTED - daughters)
        vrs['m'] = dp.measuredMass()
        ntuple.fillMom(dp.momentum(), vrs)
        # TODO: new ntuple category for daughters?
        ntuple.fill('can_dtr', vrs = vrs)

        # DavinciSmartAssociator (preferred method of truth matching)
        # Get (generated) MCParticle associated with daughter and their weights
        rels = mctool.relatedMCPs(dp)
        w = 0 # Association weights
        relp = None # type <MCParticle>
        # Loop through all relations to find best fit
        for rel in rels: 
          # Replace w if new weight is better, setting relp to the MCParticle
          # TODO: make tuple and sort instead of this
          if rel.weight() > w: 
            w = rel.weight()
            relp = rel.to()
        # Check if related particle was found before trying to fill ntuple
        if relp:
          vrs = {}
          # Fill pid, momentum (GENERATED - daughters)
          # TODO filter by pid
          # if abs(relp.particleID().pid()) in [13, 22]
          vrs['pid'] = relp.particleID().pid()
          ntuple.fillMom(relp.momentum(), vrs)
          ntuple.fill('gen_dtr', vrs = vrs)
      
      # Repeat truth matching for etas
      rels = mctool.relatedMCPs(eta)
      w = 0
      relp = None
      for rel in rels:
        if rel.weight() > w: 
          w = rel.weight()
          relp = rel.to()
      # Fill pid, momentum (GENERATED - eta)
      if relp:
        # TODO: remove this cut
        if relp.particleID().pid() in [221, 331]:
          vrs = {}
          vrs['pid'] = relp.particleID().pid()
          ntuple.fillMom(relp.momentum(), vrs)
          ntuple.fill('gen', vrs = vrs)
          ntuple.fillPvr(relp.primaryVertex(), vrs, pre = 'gen_pvr')
  except: pass

  # Fill the ntuple
  if fill: ntuple.fill()
# Close the ntuple
ntuple.close()

'''
from LoKiMC.decorators import MCID, MCABSID, MCPT, MCM, MCETA, MCP, MCPHI
from LoKiPhys.decorators import M, PX, PY, PZ, E, ID, ABSID, PROBNNe, PROBNNghost, IPCHI2, P, PT, CL, ETA, TRGHOSTPROB, PHI
###Sometimes DaVinciSmartAssociator fails (no match found at all, Particle and MCParticle are not the same type, etc.), so the backup I used is delta r matching
for pi in <TES LOC>: #same as above this would be inside event loop
    for dp in pi.daughtersVector():
        mindr = sys.float_info.max #variable that stores the minimum delta r from the loop (starts very large since we want the minimum delta r)
        relp = None
        for mcp in tes['MC/Particles']: #loop through the MCParticle TES location (probably something like 'MC/Particles')
            if MCID(mcp) == ID(dp): #requiring MCParticle and Particle to have same PID
                dphi = MCPHI(mcp) - PHI(dp)
                deta = MCETA(mcp) - ETA(dp)
                deltar = sqrt(dphi**2 + deta**2) #delta r for this MCParticle and Particle
                if deltar < mindr: mindr = deltar; relp = mcp #checking if this is smaller than the current minimum delta r and updating mindr and relp if that is the case
        #this is where you add info to your ntuple for this daughter's linked MCParticle
'''
