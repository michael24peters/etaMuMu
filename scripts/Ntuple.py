###############################################################################
# Ntuple class for eta -> mu+ mu- gamma analysis                              #
# Author: Michael Peters                                                      #
###############################################################################

from collections import OrderedDict
import ROOT, array, GaudiPython
from GaudiPython.Bindings import gbl
STD  = gbl.std
LHCB = gbl.LHCb

# Tag configuration.
TrkCats = [('ve', 1), ('tt', 2), ('it', 3), ('ot', 4), ('mu', 7)]
# Sprucing lines for Run 2.
TrgLocs = ['Hlt2ExoticaPrmptDiMuonSSTurbo', 
           'Hlt2ExoticaPrmptDiMuonTurbo',
           'Hlt2ExoticaDiMuonNoIPTurbo',
           'Hlt2ExoticaDisplDiMuon']

# =============================================================================

class Ntuple:
  """
  Class to store an ntuple.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, name, IS_MC, IS_MUMUGAMMA, tes, genTool, rftTool, pvrTool, velTool, dstTool, 
                detTool, trkTool, l0Tool, hlt1Tool, hlt2Tool):
    from collections import OrderedDict
    import ROOT, array
    self.tes      = tes
    self.genTool  = genTool
    self.rftTool  = rftTool
    self.pvrTool  = pvrTool
    self.velTool  = velTool
    self.dstTool  = dstTool
    self.detTool  = detTool
    self.trkTool  = trkTool
    self.l0Tool   = l0Tool
    self.hlt1Tool = hlt1Tool
    self.hlt2Tool = hlt2Tool
    self.saved    = {}
    self.ntuple   = OrderedDict()
    self.tfile    = ROOT.TFile(name, 'RECREATE')
    self.ttree    = ROOT.TTree('tree', 'data')
    vrsVrt = ['x', 'y', 'z', 'dx', 'dy', 'dz']
    vrsMom = ['pid', 'px', 'py', 'pz', 'e', 'p', 'pT']
    vrsMom = ['idx_gen'] + vrsMom # MC
    vrsPrt = ['iso', 'mu', 'loose_mu', 'tight_mu', 'pnn_mu', 'pnn_pi',
              'pnn_k', 'pnn_p', 'pnn_ghost', 'prb_ghost', 'ip', 'ip_chi2',
              'x0', 'y0', 'z0', 't0', 'p0', 'id0', 'z1', 'id1', 'id2',
              'id3', 'xm2', 'ym2', 'zm2']
    # TODO not used: tag_doca, tag_chi2, tag_ve_iso0, tag_ve_iso1, tag_ln_iso0, tag_ln_iso1
    # prt_iso
    vrsTag = ['m', 'dtf_m', 'dtf_dm', 'ip', 'ip_chi2', 'fd', 'fd_chi2', 
              'doca', 'dtf_chi2', 'chi2', 've_ns', 'tt_ns', 'it_ns', 
              'ot_ns', 'mu_ns']
    vrsTrg = ['l0_tos0', 'l0_tos1', 'l0_tis', 'hlt1_tos0', 'hlt1_tos1',
              'hlt1_tos2', 'hlt1_tos3', 'hlt1_tis'] + [
              'hlt2_tos%i' % i for i in range(0, len(TrgLocs))] + ['hlt2_tis']
    self.vrsInit('pvr', vrsVrt)
    self.vrsInit('tag', ['idx_pvr', 'idx_prt0', 'idx_prt1', 'idx_prt2'] +
                  vrsMom + vrsVrt + vrsTag + vrsTrg)
    self.vrsInit('tag', ['ve_iso0', 've_iso1', 'ln_iso0', 'ln_iso1'])
    self.vrsInit('prt', ['idx_pvr'] + vrsMom + vrsPrt)
    
    # MC data.
    if IS_MC: 
      vrsMcp = ['pid', 'q', 'px', 'py', 'pz', 'e', 'p', 'pT', 'x', 'y', 'z']
      self.vrsInit('mcpvr', ['x', 'y', 'z'])
      self.vrsInit('mctag', ['idx_pvr', 'idx_prt0', 'idx_prt1', 'idx_prt2', 
                              'pid_mom'] + vrsMcp)
      self.vrsInit('mcprt', ['idx_pvr', 'eta'] + vrsMcp)
      
      # Generator-level variables
      vrsGen =  ['pid', 'q', 'px', 'py', 'pz', 'e', 'p', 'pT', 'x', 'y', 'z']
      self.vrsInit('genpvr', ['x', 'y', 'z'])
      self.vrsInit('gentag', ['idx_pvr', 'idx_prt0', 'idx_prt1', 'idx_prt2',
                                'pid_mom'] + vrsGen)
      self.vrsInit('genprt', ['idx_pvr', 'eta'] + vrsGen)

      # Background variables
#      vrsBg =  ['pid', 'q', 'px', 'py', 'pz', 'e', 'p', 'pT', 'x', 'y', 'z']
#      self.vrsInit('bgpvr', ['x', 'y', 'z'])
#      self.vrsInit('bgtag', ['idx_pvr', 'idx_prt0', 'idx_prt1', 'idx_prt2',
#                              'pid_mom'] + vrsBg)
#      self.vrsInit('bgprt', ['idx_pvr'] + vrsBg) 

    self.ntuple['pvr_n'] = array.array('d', [-1])
    self.ntuple['run_n'] = array.array('d', [-1])
    self.ntuple['evt_n'] = array.array('d', [-1])
    self.ntuple['evt_tck'] = array.array('d', [-1])
    self.ntuple['evt_spd'] = array.array('d', [-1])
    
    for key, val in self.ntuple.items():
      if type(val) is array.array: self.ttree.Branch(key, val, key + '/D')
      else: self.ttree.Branch(key, val)    

  # ---------------------------------------------------------------------------

  def vrsInit(self, pre, vrs):
    import ROOT
    for v in vrs: self.ntuple['%s_%s' % (pre, v)] = ROOT.vector('double')()

  # ---------------------------------------------------------------------------

  """
  Return key for particle based on its physical properties, in priority order of 
  momentum, track, calorimeter.
  """
  def key(self, obj):
    """
    Generate the key for an object.
    """

    key = None
    try: 
      key = (obj.momentum().Px(), obj.momentum().Py(), 
             obj.momentum().Pz())
      try:
        trk = obj.proto().track()
        key = (trk.momentum().X(), trk.momentum().Y(), 
               trk.momentum().Z())
      except:
        try:
          pos = obj.proto().calo()[0].position()
          key = (pos.x(), pos.y(), pos.z(), pos.e())
        except: pass
    except:
        try: key = (obj.position().X(), obj.position().Y(),
                    obj.position().Z())
        except: pass
    return key
  
  # ---------------------------------------------------------------------------

  def close(self):
    """
    Close the ntuple.
    """

    self.tfile.Write('data', ROOT.TObject.kOverwrite)
    self.tfile.Close()

  # ---------------------------------------------------------------------------

  def clear(self):
    """
    Clear the ntuple.
    """

    self.saved.clear()
    for key, val in self.ntuple.items():
      if isinstance(val, array.array): val[0] = -1
      else: val.clear()

  # ---------------------------------------------------------------------------

  def fill(self, key = None, val = None, idx = None, vrs = None):
    """
    Fill the ntuple for either an event or an object.
    """

    if key == None or val == None: self.tfile.Cd(''); self.ttree.Fill()
    elif key in self.ntuple:
      if idx == None: self.ntuple[key].push_back(val)
      elif idx < len(self.ntuple[key]): self.ntuple[key][idx] = val

  # ---------------------------------------------------------------------------

  def fillVrt(self, pre, vrt, cov = None, pos = None):
    """
    Fill an object's vertex information
    """

    from math import sqrt
    if cov == None: cov = vrt.covMatrix()
    if pos == None: pos = vrt.position()
    self.fill('%s_x'  % pre, pos.X())
    self.fill('%s_y'  % pre, pos.Y())
    self.fill('%s_z'  % pre, pos.Z())
    self.fill('%s_dx' % pre, sqrt(abs(cov[0][0])))
    self.fill('%s_dy' % pre, sqrt(abs(cov[1][1])))
    self.fill('%s_dz' % pre, sqrt(abs(cov[2][2])))

  # ---------------------------------------------------------------------------

  def fillMom(self, pre, mom):
    """
    Fill an object's momentum
    """

    self.fill('%s_px' % pre, mom.Px())
    self.fill('%s_py' % pre, mom.Py())
    self.fill('%s_pz' % pre, mom.Pz())
    self.fill('%s_e'  % pre, mom.E())
    ## ?? decay tree fitter momentum (or mass?): attempt to fit momentum based
    ## on mother and daughter particle measurement and energy conservation. ??
    try:
        self.fill('%s_dtf_m'  % pre, mom.m().value())
        self.fill('%s_dtf_dm' % pre, mom.m().error())
    except:
        self.fill('%s_dtf_m'  % pre, -1)
        self.fill('%s_dtf_dm' % pre, -1)

  # ---------------------------------------------------------------------------

  """
  prt  : Particle object
  pvrs : primary vertices
  prts : vector<Particle>, vector of Particle objects
  """
  def fillPrt(self, prt, pvrs = None, prts = None):
    pid = prt.particleID().pid()
    pro = prt.proto() # ProtoParticle
    # Try to get raw Particle data
    try: prt = prt.data()
    except: pass
    # Try to get primary vertex associated with particle
    try: pvr = self.pvrTool.relatedPV(prt, 'Rec/Vertex/Primary')
    except: pvr = None
    # Recursive base case; check if this is a composite particle that decays.
    vrt = prt.endVertex()
    mom = prt.momentum()
    key = self.key(prt) # Generate unique key for particle.
    pre = 'tag' if vrt else 'prt'

    # Daughters.
    if pre == 'tag': # tag => candidate => has daughters
      trks = []
      idxs = []
      hits = [0]*13 # Default (to prevent segfault)
      for dtr in prt.daughters():
        (dtrPre, dtrIdx) = self.fillPrt(dtr) # Recursively loop thru daughters
        idxs += [dtrIdx]
        # Extract daughter track if able
        if dtr.proto() and dtr.proto().track():
          trk, mu = dtr.proto().track(), None
          # Store each hit in each sub-detector in a dictionary
          if trk: hits = self.hits(trk)
        # 7 -> calorimeter, so if the daughter has no hits in the calorimeter
        # but does have a muon PID, add the muon track. (Muons should skip the
        # calorimeter.)
	    # muonPID() contains info from muon system
        if hits[7] == 0 and dtr.proto().muonPID():
          mu = dtr.proto().muonPID().muonTrack()
          trks += [(trk, mu)]
      ## !! Index daughters by linking candiate to prtIdx-th daughter.
      ## Example: 'tag_idx_prt0 = 42' means the first daughter of the candidate
      ## is entry 42 in the tree. !!
      for prtIdx, dtr in enumerate(idxs):
        self.fill('%s_idx_prt%i' % (pre, prtIdx), dtr)

    # Save current index of TTree array (for building decay tree?)
    idx = self.ntuple['%s_px' % pre].size()
    self.saved[key] = idx

    if pre == 'tag':
      # Find shared hits.
      n = self.share(trks)
      # Maximum Distance of Closest Approach (DOCA) among daughters
      from LoKiArrayFunctors.decorators import AMAXDOCA
      # Use PV as constraint
      if pvr:
        dtf = GaudiPython.gbl.DecayTreeFitter.Fitter(prt, pvr)
        dtf.fit()
      else: dtf = None
      if dtf and dtf.status() == 0:
        par = dtf.fitParams(prt) # Vertex information
        mom = par.momentum() # Fitted momentum
        # Save fitted values, including uncertainty
        self.fill('%s_dtf_chi2' % pre, dtf.chiSquare())
        self.fillVrt(pre, prt, par.posCovMatrix(), par.position())
      # Use original vertex info from Particle if no dtf.
      else:
        self.fill('%s_dtf_chi2' % pre, -1)
        self.fillVrt(pre, vrt)

    # Momentum and mass.
    self.fill('%s_m' % pre, prt.measuredMass())
    self.fillMom(pre, mom)
    self.fill('%s_p' % pre, prt.p())
    self.fill('%s_pT' % pre, prt.pt())

    # Trigger.
    # TOS == Trigger On Signal, TIS == Trigger Independent of Signal,
    # TOB == Trigger On Both.
    # setOfflineInput(prt) sets candidate to analyze for the trigger tool.
    self.l0Tool.setOfflineInput(prt)
    self.hlt1Tool.setOfflineInput(prt)
    self.hlt2Tool.setOfflineInput(prt)
    # setTriggerInput(trigger) sets which decision to evaluate.
    # Check if candidate fired L0DiMuon trigger.
    self.l0Tool.setTriggerInput('L0DiMuonDecision')
    trg = self.l0Tool.tisTosTobTrigger()
    self.fill('%s_l0_tos0' % pre, trg.tos())
    # Repeat for single muon.
    self.l0Tool.setTriggerInput('L0MuonDecision')
    trg = self.l0Tool.tisTosTobTrigger()
    self.fill('%s_l0_tos1' % pre, trg.tos())
    self.l0Tool.setTriggerInput('L0HadronDecision')
    trg = self.l0Tool.tisTosTobTrigger()
    self.fill('%s_l0_tis' % pre, trg.tis())
    # HLT1 software trigger
    self.hlt1Tool.setTriggerInput('Hlt1DiMuonNoIPDecision')
    trg = self.hlt1Tool.tisTosTobTrigger()
    self.fill('%s_hlt1_tos0' % pre, trg.tos())
    self.hlt1Tool.setTriggerInput('Hlt1MultiDiMuonNoIPDecision')
    trg = self.hlt1Tool.tisTosTobTrigger()
    self.fill('%s_hlt1_tos1' % pre, trg.tos())
    self.hlt1Tool.setTriggerInput('Hlt1DiMuonLowMassDecision')
    trg = self.hlt1Tool.tisTosTobTrigger()
    self.fill('%s_hlt1_tos2' % pre, trg.tos())
    self.hlt1Tool.setTriggerInput('Hlt1DiMuonHighMassDecision')
    trg = self.hlt1Tool.tisTosTobTrigger()
    self.fill('%s_hlt1_tos3' % pre, trg.tos())
    self.hlt1Tool.setTriggerInput('Hlt1.*TrackMVA.*')
    trg = self.hlt1Tool.tisTosTobTrigger()
    self.fill('%s_hlt1_tis' % pre, trg.tis())
    # HLT2 software trigger
    self.hlt2Tool.setTriggerInput('Hlt2Topo.*')
    trg = self.hlt2Tool.tisTosTobTrigger()
    self.fill('%s_hlt2_tis' % pre, trg.tis())
    for locIdx, loc in enumerate(TrgLocs):
      self.hlt2Tool.setTriggerInput(loc + 'Decision')
      trg = self.hlt2Tool.tisTosTobTrigger()
      self.fill('%s_hlt2_tos%i' % (pre, locIdx), trg.tos())

    # Particle ID.
    self.fill('%s_pid' % pre, pid)
    # If there is a proto track, try to fill what information is available.
    if pro:
      # Fill muon information
      try:
        # Muon hypothesis
        self.fill('%s_mu' % pre, pro.muonPID().IsMuon())
        # More inclusive muon hypothesis
        self.fill('%s_loose_mu' % pre, pro.muonPID().IsMuonLoose())
        # More exclusive muon hypothesis
        self.fill('%s_tight_mu' % pre, pro.muonPID().IsMuonTight())
      except: pass
      # Likelihood particle is muon, pion, kaon, proton, or ghost particle.
      self.fill('%s_pnn_mu'    % pre, pro.info(701, -100))
      self.fill('%s_pnn_pi'    % pre, pro.info(702, -100))
      self.fill('%s_pnn_k'     % pre, pro.info(703, -100))
      self.fill('%s_pnn_p'     % pre, pro.info(704, -100))
      self.fill('%s_pnn_ghost' % pre, pro.info(705, -100))
      self.fill('%s_prb_ghost' % pre, pro.track().ghostProbability()
                if pro.track() else -100) # Info unavailable

      # Track
      trk = pro.track()
      ids = []
      # Fill VELO hits.
      if trk:
        for i in trk.lhcbIDs():
          if i.isVelo():
	    ## ?? What does detTool do? sensor()?
            d = self.detTool.sensor(i.veloID())
            ids += [(d.z(), d, i)]
          ids.sort()
      # Up to 4 VELO hits allowed per track. If less found, fill rest with -1.
      for hit in range(0, 4):
        if hit >= len(ids):
          self.fill('%s_x%i'  % (pre, hit), -1)
          self.fill('%s_y%i'  % (pre, hit), -1)
          self.fill('%s_z%i'  % (pre, hit), -1)
          self.fill('%s_id%i' % (pre, hit), -1)
          self.fill('%s_t%i'  % (pre, hit),  -1)
          self.fill('%s_p%i'  % (pre, hit),  -1)
          continue
        ## ?? ADDITIONAL EXPLANATION NEEDED ??
        # Calculate track parameters based on z position
        z, d, i = ids[hit]
        s = i.veloID()
        v = GaudiPython.gbl.LHCb.StateVector()
        self.trkTool.propagate(trk, z, v, prt.particleID())
        # Fill track coordinate vertex info
        self.fill('%s_x%i'  % (pre, hit), v.x())
        self.fill('%s_y%i'  % (pre, hit), v.y())
        self.fill('%s_z%i'  % (pre, hit), v.z())
        # Phi type strip == -1, r type strip == 1.
        self.fill('%s_id%i' % (pre, hit), i.lhcbID()*
                  (-1 if s.isPhiType() else 1))
        # Fill spatial coordinate and strip width
        if (s.isPhiType()):
          self.fill('%s_t%i' % (pre, hit), d.globalPhi(s.strip(), 0))
          self.fill('%s_p%i' % (pre, hit), d.phiPitch(s.strip()))
        if (s.isRType()):
          self.fill('%s_t%i' % (pre, hit), d.globalR(s.strip(), 0))
          self.fill('%s_p%i' % (pre, hit), d.rPitch(s.strip()))

      ## ?? What does this do? ??
      # Muon extrapolation.
      z = 15270.0 # Muon system region
      # Propagated state vector
      v = GaudiPython.gbl.LHCb.StateVector()
      # Save extrapolated position
      ### ?? ADDED THIS CHECK TO PREVENT SEGFAULT ERROR ??
      if trk:
        self.trkTool.propagate(trk, z, v, prt.particleID())
        self.fill('%s_xm2' % pre, v.x())
        self.fill('%s_ym2' % pre, v.y())
        self.fill('%s_zm2' % pre, v.z())

    # Find all linked MC particle matches.
    try:
      # Relate reconstructed particle to generator-level particle.
      gen = None; wgt = 0; rels = self.genTool.relatedMCPs(prt)
      # Select match with heighest weight
      for rel in rels: gen = rel.to() if rel.weight() > wgt else gen
      # Look at daughters of closest match gen-level particle. If they are
      # muons or photons, save them.
      if gen: (genPre, genIdx) = self.fillMcp(gen, pre)
      else: genIdx = -1
      self.fill('%s_idx_gen' % pre, genIdx)
    except: pass

    # IP.
    from ctypes import c_double
    ip, ipChi2 = c_double(-1.0), c_double(-1.0)
    self.dstTool.distance(prt, pvr, ip, ipChi2)
    # Compute Impact Paramter (IP)
    # Measures shortest distance from particle to PV. High IP means particle
    # likely isn't prompt (PV) eta and instead displaced (SV) other particle
    # (e.g. B or D meson) producing decay that fakes an eta.
    if pvr: self.dstTool.distance(prt, pvr, ip, ipChi2)
    self.fill('%s_ip' % pre, ip.value) # need to convert c_double to python float
    self.fill('%s_ip_chi2' % pre, ipChi2.value) # same here
    fd, fdChi2 = c_double(-1.0), c_double(-1.0)
    # Compute Flight Distance (FD)
    # Measures how far reconstructed candidate traveled from PV to SV. Should be
    # approximately zero since prompt eta decays are basically at production
    # point. High FD means desired decay products likely coming from LLP (e.g.
    # B or D meson) and eta candidate is random combo of displaced tracks that
    # happen to look right. Or it could mean the tracks were mis-constructed and
    # really should have been identified as an eta candidate.
    if pvr and vrt: self.dstTool.distance(vrt, pvr, fd, fdChi2)
    self.fill('%s_fd' % pre, fd.value) # need to convert c_double to python float
    self.fill('%s_fd_chi2' % pre, fdChi2.value) # same here

    # Isolation.
    # !! Removed for now.

    # Primary vertex.
    if pvr:
        key = self.key(pvr)
        if not key in self.saved:
            self.saved[key] = self.ntuple['pvr_x'].size()
            self.fillVrt('pvr', pvr)
        # Save index of primary vertex
        self.fill('%s_idx_pvr' % pre, self.saved[key])
    else: self.fill('%s_idx_pvr' % pre, -1)

    # Return prefix used for labeling and index where the particle's data
    # is stored.
    return (pre, idx)

  # ---------------------------------------------------------------------------

  """
  hits() and share() combined find the tracks that have common hits. If they
  have many shared hits, they may be an identical track and thus are thus a
  fake particle. This is useful for 'clone rejection'.
  """
  def hits(self, trk1, trk2 = None, mu1 = None, mu2 = None, dct = None):
      if dct == None: dct = {n:0 for n in range(0, 13)}
      if not trk1: return dct
      id1s = [id1 for id1 in trk1.lhcbIDs()]
      id2s = [id2 for id2 in trk2.lhcbIDs()] if trk2 else []
      if mu1: id1s += [id1 for id1 in mu1.lhcbIDs()]
      if mu2: id2s += [id2 for id2 in mu2.lhcbIDs()]
      for id1 in id1s:
          cat1 = id1.detectorType()
          chn1 = id1.channelID()
          if not trk2: dct[cat1] += 1; dct[0] += 1
          for id2 in id2s:
              cat2 = id2.detectorType()
              chn2 = id2.channelID()
              if chn1 == chn2:
                  dct[0] += 1
                  if cat1 == cat2: dct[cat1] += 1
      return dct

  # ---------------------------------------------------------------------------

  def share(self, trks):
      ns = {n:0 for n in range(0, 13)}
      for trk1 in trks:
          for trk2 in trks:
              if trk1[0] == trk2[0]: continue
              hits = self.hits(trk1[0], trk2[0], trk1[1], trk2[1])
              for key, val in hits.items():
                  if val > ns[key]: ns[key] = val
      for cat in TrkCats:
          n = ns[cat[1]]
          self.fill('tag_%s_ns' % (cat[0]), n)
      return ns[TrkCats[0][1]]

  # ---------------------------------------------------------------------------

  def pseudorapidity(self, prt):
    from math import log
    p = prt.p() # |p|
    pz = prt.momentum().Pz()
    return 0.5 * log((p + pz) / (p - pz))

  # ---------------------------------------------------------------------------


  # Fill all MC-matched particles (reconstruction -> generator-level mapping).
  def fillMcp(self, prt, rec):
    pid = prt.particleID().pid()
    mom = prt.momentum()
    pos = None
    key = self.key(prt)

    pre = 'mctag' if rec == 'tag' else 'mcprt'
    # Prevent filling same particle more than once.
    if key in self.saved: return (pre, self.saved[key])

    # Save current index of TTree array, get next row index of TTree.
    idx = self.ntuple['%s_px' % pre].size()
    self.saved[key] = idx

    # Momentum.
    self.fillMom(pre, mom)
    self.fill('%s_p' % pre, prt.p())
    self.fill('%s_pT' % pre, prt.pt())

    # PID.
    self.fill('%s_q' % pre, float(prt.particleID().threeCharge()) / 3.0)
    self.fill('%s_pid' % pre, pid)
    # Fill PID of parent if possible
    try: self.fill('%s_pid_mom' % pre, prt.originVertex().mother().
                    particleID().pid())
    # Failure could be from primary particles or incorrect MC information
    except: self.fill('%s_pid_mom' % pre, 0)

    # Vertex.
    if not pos: pos = prt.originVertex().position()
    self.fill('%s_x' % pre, pos.X())
    self.fill('%s_y' % pre, pos.Y())
    self.fill('%s_z' % pre, pos.Z())

    # Primary vertex.
    pvr = prt.primaryVertex()
    if pvr:
      key = self.key(pvr)
      if not key in self.saved:
        self.saved[key] = self.ntuple['mcpvr_x'].size()
        self.fill('mcpvr_x', pvr.position().X())
        self.fill('mcpvr_y', pvr.position().Y())
        self.fill('mcpvr_z', pvr.position().Z())
      self.fill('%s_idx_pvr' % pre, self.saved[key])
    else: self.fill('%s_idx_pvr' % pre, -1)

#    # TODO: This will not catch falsely identified etas
#    # e.g. was actually a D meson but was identified as eta
#    if pre == 'mctag' and pid != 221: fillBkgd(prt, rec)
#    # TODO: This will not catch falsely identified daughters
#    # e.g. two rec muons actually one muon, falsely identified photon, etc.
#    elif pre == 'mcprt' and pid not in [-13, 13, 22]: fillBkgd(prt, rec)

    return (pre, idx)

  # ---------------------------------------------------------------------------

  # Save all gen-level occurrences of eta -> mu+ mu- gamma
  def fillGen(self, prt):
    pid = prt.particleID().pid()
    mom = prt.momentum()
    pos = None
    key = self.key(prt)
    
    # Check for LHCb acceptance range
    eta = self.pseudorapidity(prt)
    if eta < 2.0 or eta > 4.5: return (None, None)

    # For MC data, we only want to save eta -> mu+ mu- (gamma) occurrences, but
    # only particles with pid 221 getting initially passed in, so this is a 
    # sufficient check.
    pre = 'gentag' if pid in [221] else 'genprt'
    # Prevent filling same particle more than once.
    if key in self.saved: return (pre, self.saved[key])

    # Daughters.
    if pre == 'gentag':
      dtrs = []
      idxs = []
      
      # Loop over decay vertices
      for vrt in prt.endVertices():
        # Loop over products of vertices (daughters)
        for dtr in vrt.products():
          if abs(dtr.particleID().pid()) in [13, 22]: 
            dtrs.append({'pid': dtr.particleID().pid(), 'dtr' : dtr})
        # Sort daughters by PID, i.e. always in {-13, 13, 22} order.
        dtrs = sorted(dtrs, key=lambda d: d['pid'])
      # Skip all etas which do not decay exactly to mu+ mu- gamma
      pids = [d['pid'] for d in dtrs]
      if pids != [-13, 13, 22]: return (None, None)
      for dtr in dtrs:
        # Skip all etas which do not have muons with p >= 3 GeV
        if dtr['pid'] in [-13, 13] and dtr['dtr'].p() < 3000: return (None, None)
        # Skip all eta which do not have daughters with pT >= 500 MeV
#        if dtr['dtr'].pt() < 500: return (None, None)
      # Otherwise recursively fill daughters
      for dtr in dtrs:
        (dtrPre, dtrIdx) = self.fillGen(dtr['dtr'])
        self.fill(dtrIdx, dtr['dtr'])
        idxs += [dtrIdx]
      
      # Fill daughter index
      for prtIdx, dtr in enumerate(idxs):
        self.fill('%s_idx_prt%i' % (pre, prtIdx), dtr)

    # Save current index of TTree array, get next row index of TTree.
    idx = self.ntuple['%s_px' % pre].size()
    self.saved[key] = idx

    # Momentum.
    self.fillMom(pre, mom)
    self.fill('%s_p' % pre, prt.p())
    self.fill('%s_pT' % pre, prt.pt())

    # Pseudorapidity
    self.fill('%s_eta' % pre, self.pseudorapidity(prt))

    # PID.
    self.fill('%s_q' % pre, float(prt.particleID().threeCharge()) / 3.0)
    self.fill('%s_pid' % pre, pid)
    # Fill PID of parent if possible
    try: self.fill('%s_pid_mom' % pre, prt.originVertex().mother().
                    particleID().pid())
    # Failure could be from primary particles or incorrect MC information
    except: self.fill('%s_pid_mom' % pre, 0)

    # Vertex.
    if not pos: pos = prt.originVertex().position()
    self.fill('%s_x' % pre, pos.X())
    self.fill('%s_y' % pre, pos.Y())
    self.fill('%s_z' % pre, pos.Z())

    # Primary vertex.
    pvr = prt.primaryVertex()
    if pvr:
      key = self.key(pvr)
      if not key in self.saved:
        self.saved[key] = self.ntuple['mcpvr_x'].size()
        self.fill('mcpvr_x', pvr.position().X())
        self.fill('mcpvr_y', pvr.position().Y())
        self.fill('mcpvr_z', pvr.position().Z())
      self.fill('%s_idx_pvr' % pre, self.saved[key])
    else: self.fill('%s_idx_pvr' % pre, -1)

    return (pre, idx)

  # ---------------------------------------------------------------------------

  # Save generator-level information for background particles around the invar-
  # -iant mass range of eta (547 MeV)
  def fillBkgd(self, prt, rec):
    pid = prt.particleID().pid()
    mom = prt.momentum()
    pos = None
    key = self.key(prt)

    pre = 'bgtag' if rec == 'tag' else 'bgprt'
    # Prevent filling same particle more than once.
    if key in self.saved: return (pre, self.saved[key])

    # Daughters.
    if pre == 'bgtag':
      dtrs = []
      idxs = []
      
      # Loop over decay vertices
      for vrt in prt.endVertices():
        # Loop over products of vertices (daughters)
        for dtr in vrt.products():
          # Save all daughters from composite particle
          dtrs.append({'pid': dtr.particleID().pid(), 'dtr' : dtr})
        # Sort daughters by PID, i.e. always in {-13, 13, 22} order.
        dtrs = sorted(dtrs, key=lambda d: d['pid'])
      pids = [d['pid'] for d in dtrs]
      # Save decays that weren't exactly eta -> mu+ mu- gamma (i.e. background)
      if pid == 221 and pids == [-13, 13, 22]: return
      # If check passed, recursively fill daughters
      for dtr in dtrs:
        (dtrPre, dtrIdx) = self.fillGen(dtr['dtr'])
        self.fill(dtrIdx, dtr['dtr'])
        idxs += [dtrIdx]
      
      # Fill daughter index
      for prtIdx, dtr in enumerate(idxs):
        self.fill('%s_idx_prt%i' % (pre, prtIdx), dtr)

    # Save current index of TTree array, get next row index of TTree.
    idx = self.ntuple['%s_px' % pre].size()
    self.saved[key] = idx

    # Momentum.
    self.fillMom(pre, mom)
    self.fill('%s_p' % pre, prt.p())
    self.fill('%s_pt' % pre, prt.pt())

    # PID.
    self.fill('%s_q' % pre, float(prt.particleID().threeCharge()) / 3.0)
    self.fill('%s_pid' % pre, pid)
    # Fill PID of parent if possible
    try: self.fill('%s_pid_mom' % pre, prt.originVertex().mother().
                    particleID().pid())
    # Failure could be from primary particles or incorrect MC information
    except: self.fill('%s_pid_mom' % pre, 0)

    # Vertex.
    if not pos: pos = prt.originVertex().position()
    self.fill('%s_x' % pre, pos.X())
    self.fill('%s_y' % pre, pos.Y())
    self.fill('%s_z' % pre, pos.Z())

    # Primary vertex.
    pvr = prt.primaryVertex()
    if pvr:
      key = self.key(pvr)
      if not key in self.saved:
        self.saved[key] = self.ntuple['mcpvr_x'].size()
        self.fill('mcpvr_x', pvr.position().X())
        self.fill('mcpvr_y', pvr.position().Y())
        self.fill('mcpvr_z', pvr.position().Z())
      self.fill('%s_idx_pvr' % pre, self.saved[key])
    else: self.fill('%s_idx_pvr' % pre, -1)

    return (pre, idx)


