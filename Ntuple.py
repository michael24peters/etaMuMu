###############################################################################
# Ntuple class for jet substructure analysis                                  #
# Author: Jake Pfaller                                                        #
###############################################################################

from collections import OrderedDict
import ROOT, array, GaudiPython
from GaudiPython.Bindings import gbl
STD  = gbl.std
LHCB = gbl.LHCb

# =============================================================================

class Ntuple:
  """
  Class to store an ntuple.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, name, tes, toolSvc, detSvc):
    """
    Initialize the ntuple with the needed tools.
    """

    ROOT.gInterpreter.ProcessLine(
      "#include \"/cvmfs/lhcb.cern.ch/lib/lhcb/HLT/"
      "HLT_v25r4/Hlt/HltDisplVertices/Kernel/"
      "IMatterVeto.h\"")
    self.pvrTool = toolSvc.create(
      'GenericParticle2PVRelator<_p2PVWithIPChi2, '
      'OfflineDistanceCalculatorName>/P2PVWithIPChi2',
      interface = 'IRelatedPVFinder')
    self.tagTool = toolSvc.create(
      'LoKi::BDTTag',
      interface = 'IJetTagTool')
    self.genTool = toolSvc.create(
      'DaVinciSmartAssociator',
      interface = 'IParticle2MCWeightedAssociator')
    self.mtrTool = toolSvc.create(
      'MatterVetoTool',
      interface = 'IMatterVeto')
    self.velTool = toolSvc.create(
      'VeloExpectation',
      interface = 'IVeloExpectation')
    self.dstTool = toolSvc.create(
      'LoKi::TrgDistanceCalculator',
      interface = 'IDistanceCalculator')
    self.trkTool = toolSvc.create(
      'TrackMasterExtrapolator',
      interface = 'ITrackExtrapolator')
    self.detTool = detSvc[
      '/dd/Structure/LHCb/BeforeMagnetRegion/Velo']
    self.name   = name
    self.tes    = tes
    self.saved  = {}
    self.ntuple = OrderedDict()
    self.tfile  = ROOT.TFile(self.name, 'RECREATE')
    self.ttree  = ROOT.TTree('data', 'data')
    self.vrs    = {}

  # ---------------------------------------------------------------------------

  def init(self, pre = '', vrs = []):
    """
    Initialize a set of variables, or intitialize the tuple
    """

    if ((pre == '') and (vrs == [])):
      for key, val in self.ntuple.items():
        if type(val) is array.array: self.ttree.Branch(key, val, key + '/D')
        else: self.ttree.Branch(key, val)
    else:
      self.saved[pre] = {}
      self.vrs[pre]   = vrs
      for v in vrs: self.ntuple['%s_%s' % (pre, v)] = ROOT.vector('double')()

  # ---------------------------------------------------------------------------

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

    self.tfile.Write(); self.tfile.Close()

  # ---------------------------------------------------------------------------

  def clear(self):
    """
    Clear the ntuple.
    """

    for key, val in self.saved.items(): val.clear()
    for key, val in self.ntuple.items():
      if type(val) is array.array: val[0] = -1
      else: val.clear()

  # ---------------------------------------------------------------------------

  def fill(self, key = None, val = None, idx = None, vrs = None):
    """
    Fill the ntuple for either an event or an object.
    """

    if key == None and val == None and idx == None and vrs == None:
      self.tfile.Cd(''); self.ttree.Fill()
    elif vrs != None and key != None:
      pre = key
      for key in self.vrs[pre]:
        val = vrs[key] if key in vrs else -1
        self.ntuple[pre + '_' + key].push_back(val)
    elif key in self.ntuple: 
      if idx == None: self.ntuple[key].push_back(val)
      elif idx < len(self.ntuple[key]): self.ntuple[key][idx] = val

  # ---------------------------------------------------------------------------

  def fillMom(self, obj, vrs):
    """
    Fill an object's momentum
    """
    if not obj: return
    vrs['px'] = obj.Px()
    vrs['py'] = obj.Py()
    vrs['pz'] = obj.Pz()
    vrs['e' ] = obj.E()

  # ---------------------------------------------------------------------------

  def fillPid(self, obj, vrs):
    """
    Fill an object's particle ID and charge
    """

    if not obj: return
    vrs['pid'] = obj.pid()

  # ---------------------------------------------------------------------------

  def fillPro(self, obj, vrs):
    """
    Fill protoparticle info
    """

    if not obj: return
    if obj.muonPID(): vrs['is_mu'] = obj.muonPID().IsMuon()
    vrs['pnn_e' ] = obj.info(700, -100)
    vrs['pnn_mu'] = obj.info(701, -100)
    vrs['pnn_pi'] = obj.info(702, -100)
    vrs['pnn_k' ] = obj.info(703, -100)
    vrs['pnn_p' ] = obj.info(704, -100)
    vrs['ecal']   = obj.info(332, -100)
    vrs['hcal']   = obj.info(333, -100)

  # ---------------------------------------------------------------------------

  def fillTrk(self, obj, pid, vrs):
    """
    Fill track info
    """

    if not obj or not pid: return
    vrs['prb_ghost'] = obj.ghostProbability()
    vrs['vid'] = 1.0
    vids = []
    for vid in obj.lhcbIDs():
      if vid.isVelo():
        vrs['vid'] *= float(vid.veloID().channelID())/1000000.0
        vids += [(self.detTool.sensor(vid.veloID()).z(), 
                  vid.veloID().channelID())]
    vids.sort()
    sta = LHCB.StateVector()
    if len(vids) > 0: self.trkTool.propagate(obj, vids[0][0], sta, pid)
    vrs['x'] = sta.x()
    vrs['y'] = sta.y()
    vrs['z'] = sta.z()

  # ---------------------------------------------------------------------------

  def fillDst(self, obj, pvr, dst, vrs):
    """
    Fill object distance
    """

    if not obj or not pvr: return
    val, valChi2 = ROOT.Double(-1), ROOT.Double(-1)
    self.dstTool.distance(obj, pvr, val, valChi2)
    vrs[dst] = val
    vrs[dst + '_chi2'] = valChi2

  # ---------------------------------------------------------------------------

  def fillPos(self, obj, vrs):
    """
    Fill object position
    """

    if not obj: return
    pos = obj.position()
    vrs['x'] = pos.X()
    vrs['y'] = pos.Y()
    vrs['z'] = pos.Z()
    vrs['in_mtr'] = self.mtrTool.isInMatter(pos)

  # ---------------------------------------------------------------------------

  def fillCov(self, obj, vrs):
    """
    Fill object covariance matrix
    """

    try:
      cov = obj.covMatrix()
      vrs['ndof'] = obj.nDoF()
      vrs['chi2'] = obj.chi2()
      vrs['dx']   = cov[0][0]
      vrs['dy']   = cov[1][1]
      vrs['dz']   = cov[2][2]
    except: pass
  
  # ---------------------------------------------------------------------------

  def fillPvr(self, obj, vrs, pre = 'gen_pvr'):
    """
    Fill object primary vertex index
    """

    if not obj: return # obj should be of type Vertex
    vrs['idx_pvr'] = self.addPvr(obj, pre)

  # ---------------------------------------------------------------------------
  
  def fillGen(self, obj, vrs):
    """
    Relate an MC particle to a Rec particle, then add the MC particle
    """

    if not obj: return
    vrs['idx_gen'] = self.addGen(gen)

  # ---------------------------------------------------------------------------

  def addGen(self, obj, pre = 'gen'):
    """
    Add an MC particle to the tuple
    """

    key = self.key(obj)
    if key in self.saved[pre]: return self.saved[pre][key]
    vrs = {}
    idx = len(self.saved[pre])
    self.fillPid(obj.particleID(), vrs)
    self.fillMom(obj.momentum(), vrs)
    self.saved[pre][key] = idx
    self.fill(pre, vrs = vrs)
    return idx

  # ---------------------------------------------------------------------------

  def addPvr(self, obj, pre):
    """
    Add primary vertex to the tuple
    """

    # Check if key already saved
    key = self.key(obj)
    if key in self.saved[pre]: return self.saved[pre][key]
    vrs = {}
    idx = len(self.saved[pre])
    # Fill primary vertex for generated event
    if pre == 'gen_pvr':
      self.fillPos(obj, vrs)
      self.fillCov(obj, vrs)
    # Fill primary vertex for reconstructed event
    elif pre == 'rec_pvr':
      obj = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary') # Get pv for associated rec_pvr Particle
      self.fillPos(obj, vrs)
      self.fillCov(obj, vrs)
    # Save key to prevent duplicates
    self.saved[pre][key] = idx
    # Fill primary vertex in ntuple
    self.fill(pre, vrs = vrs)
    return idx
  
  # ---------------------------------------------------------------------------

  def addTrk(self, obj, pre = 'trk'):
    """
    Add a track to the tuple
    """
    
    key = self.key(obj)
    if key in self.saved[pre]: return self.saved[pre][key]
    vrs = {}
    idx = len(self.saved[pre])
    pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
    self.fillMom(obj.momentum(), vrs)
    self.fillPid(obj.particleID(), vrs)
    self.fillPro(obj.proto(), vrs)
    self.fillDst(obj, pvr, 'ip', vrs)
    self.fillTrk(obj.proto().track(), obj.particleID(), vrs)
    self.fillPvr(pvr, vrs)
    self.fillGen(obj, vrs)
    self.saved[pre][key] = idx
    self.fill(pre, vrs = vrs)
    return idx

    