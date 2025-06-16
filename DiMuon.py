# Data.
from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles(['/data/dst/CL16.MD.Reco16MINIBIAS.0.00.dst'],
#                            clear = True) # CL16
#IOHelper('ROOT').inputFiles(['/data/dst/CL16.MD.Reco16DIMUON.0.00.dst'],
#                            clear=True)   # CL16
#IOHelper('ROOT').inputFiles(['/data/dst/CL16.MD.Turbo03LEPTONS.0.00.dst'],
#                            clear=True)   # Turbo03
#IOHelper('ROOT').inputFiles(['/data/dst/MC16.MD.40114002.0.00.dst'],
#                            clear=True)    # MC16
#IOHelper('ROOT').inputFiles(['/data/dst/MC16.MD.40114017.0.00.dst'],
#                            clear=True)   # MC16
#IOHelper('ROOT').inputFiles(['/data/dst/MC16.MD.40114007.0.00.dst'],
#                            clear=True)   # MC16
#IOHelper('ROOT').inputFiles(['/data/dst/MC16.MD.42122060.0.00.ldst'],
#                            clear=True)   # MC16

Type = 'Turbo03'
Mode = 'mue'

# DaVinci configuration.
from Configurables import DaVinci
DaVinci().DataType   = '2016'
DaVinci().Lumi = True
DaVinci().TupleFile = 'lumi.root'
DaVinci().DDDBtag    = 'dddb-20150724'
DaVinci().CondDBtag  = 'cond-20160522'
DaVinci().DQFLAGStag = 'dq-20150717'
if 'MC' in Type:
    DaVinci().Simulation = True
    DaVinci().Lumi = False
    DaVinci().DDDBtag = 'dddb-20150724'
    DaVinci().CondDBtag = 'sim-20160820-vc-md100'
if 'Turbo' in Type: 
    from Configurables import DstConf, TurboConf
    from PhysConf.Filters import LoKi_Filters
    hlt = LoKi_Filters(HLT2_Code = 
                       "HLT_PASS_RE('.*Hlt2Exotica.*TurboDecision.*')")
    DstConf().Turbo = True
    TurboConf().PersistReco = True
    DaVinci().EventPreFilters = hlt.filters('TriggerFilters')

# Tag configuration.
TrkCats = [('ve', 1), ('tt', 2), ('it', 3), ('ot', 4), ('mu', 7)]
TrgLocs = ['Hlt2ExoticaPrmptDiMuonSSTurbo', 
           'Hlt2ExoticaPrmptDiMuonTurbo',
           'Hlt2ExoticaDiMuonNoIPTurbo',
           'Hlt2ExoticaQuadMuonNoIP',
           'Hlt2ExoticaDisplDiMuon']

# Create the di-muons.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence
from StandardParticles import StdAllLooseMuons
cmbDimus = CombineParticles(
    'CmbDimus',
    DecayDescriptor = "KS0 -> mu+ mu-",
    CombinationCut = ("(AMAXDOCA('') < 0.2*mm) & "
                      "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) &"
                      "(AMINCHILD('mu-' == ABSID, P) > 10*GeV) &"
                      "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.5) &"
                      "(AMINCHILD('mu-' == ABSID, PT) > 0.5*GeV)"),
    MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")
selDimus = Selection(
    'SelDimus',
    Algorithm = cmbDimus, 
    RequiredSelections = [StdAllLooseMuons])
seqDimus = SelectionSequence('SeqDimus', TopSelection = selDimus)
if 'mu' in Mode: DaVinci().appendToMainSequence([seqDimus])

# Create the same-sign di-muons.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence
from StandardParticles import StdAllLooseMuons
cmbSsmus = CombineParticles(
    'CmbSsmus',
    DecayDescriptor = "[KS0 -> mu- mu-]cc",
    CombinationCut = ("(AMAXDOCA('') < 0.2*mm) & "
                      "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) &"
                      "(AMINCHILD('mu-' == ABSID, P) > 10*GeV) &"
                      "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.5) &"
                      "(AMINCHILD('mu-' == ABSID, PT) > 0.5*GeV)"),
    MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")
selSsmus = Selection(
    'SelSsmus',
    Algorithm = cmbSsmus, 
    RequiredSelections = [StdAllLooseMuons])
seqSsmus = SelectionSequence('SeqSsmus', TopSelection = selSsmus)
if Mode == 'mu' and 'MC' in Type: DaVinci().appendToMainSequence([seqSsmus])

# Create the di-electrons.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence
from StandardParticles import StdAllLooseElectrons
cmbDiels = CombineParticles(
    'CmbDiels',
    DecayDescriptor = "Lambda0 -> e+ e-",
    CombinationCut = ("(AMAXDOCA('') < 0.2*mm) & "
                      "(AMAXCHILD('e-' == ABSID, TRCHI2DOF) < 3) &"
                      "(AMINCHILD('e-' == ABSID, P) > 1*GeV) &"
                      "(AMINCHILD('e-' == ABSID, PROBNNe) > 0.5) &"
                      "(AMINCHILD('e-' == ABSID, PT) > 0.5*GeV)"),
    MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")
selDiels = Selection(
    'SelDiels',
    Algorithm = cmbDiels, 
    RequiredSelections = [StdAllLooseElectrons])
seqDiels = SelectionSequence('SeqDiels', TopSelection = selDiels)
if Mode == 'mue': DaVinci().appendToMainSequence([seqDiels])

# Create the quad-muons.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence
cmbQdmus = CombineParticles(
    'CmbQdmus',
    DecayDescriptor = "B0 -> KS0 KS0",
    MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")
selQdmus = Selection(
    'SelQdmus',
    Algorithm = cmbQdmus, 
    RequiredSelections = [selDimus])
seqQdmus = SelectionSequence('SeqQdmus', TopSelection = selQdmus)
if Mode == 'mu': DaVinci().appendToMainSequence([seqQdmus])

# Create the quad-electron/muons.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence, DataOnDemand
tesDimus = DataOnDemand(Location = '/Event/Turbo/Hlt2ExoticaPrmpt'
                            'DiMuonTurbo/Particles')
cmbQdlps = CombineParticles(
    'CmbQdlps',
    DecayDescriptor = "B0 -> KS0 Lambda0",
    MotherCut = "(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")
selQdlps = Selection(
    'SelQdlps',
    Algorithm = cmbQdlps, 
    RequiredSelections = [tesDimus, selDiels])
seqQdlps = SelectionSequence('SeqQdlps', TopSelection = selQdlps)
if Mode == 'mue': DaVinci().appendToMainSequence([seqQdlps])

# Create the Kss.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence
from StandardParticles import StdAllLoosePions
cmbKss = CombineParticles(
    'CmbKss',
    DecayDescriptor = "KS0 -> pi+ pi-",
    CombinationCut = ("(ADAMASS('KS0') < 50*MeV)"
                      " & (AMINCHILD(BPVIPCHI2(), ABSID == 'pi+') > 16)"
                      " & (AMINCHILD(PT,          ABSID == 'pi+') > 0.5*GeV)"),
    MotherCut = "(VFASPF(VCHI2PDOF) < 25)")
selKss = Selection(
    'SelKss',
    Algorithm = cmbKss, 
    RequiredSelections = [StdAllLoosePions])
seqKss = SelectionSequence('SeqKss', TopSelection = selKss)
if Mode == 'Ks': DaVinci().appendToMainSequence([seqKss])

# Create the di-pions.
from Configurables import CombineParticles
from PhysSelPython.Wrappers import Selection, SelectionSequence
from StandardParticles import StdAllLoosePions
cmbDipss = CombineParticles(
    'CmbDipss',
    DecayDescriptor = "[KS0 -> pi+ pi+]cc",
    CombinationCut = ("(AMINCHILD(PT, ABSID == 'pi+') > 0.5*GeV)"),
    MotherCut = "(VFASPF(VCHI2PDOF) < 25)")
selDipss = Selection(
    'SelDipss',
    Algorithm = cmbDipss, 
    RequiredSelections = [StdAllLoosePions])
seqDipss = SelectionSequence('SeqDipss', TopSelection = selDipss)
cmbDipos = CombineParticles(
    'CmbDipos',
    DecayDescriptor = "KS0 -> pi+ pi-",
    CombinationCut = ("(AMINCHILD(PT, ABSID == 'pi+') > 0.5*GeV)"),
    MotherCut = "(VFASPF(VCHI2PDOF) < 25)")
selDipos = Selection(
    'SelDipos',
    Algorithm = cmbDipos, 
    RequiredSelections = [StdAllLoosePions])
seqDipos = SelectionSequence('SeqDipos', TopSelection = selDipos)
if Mode == 'pi': DaVinci().appendToMainSequence([seqDipss, seqDipos])

# Create the jets.
from Configurables import HltParticleFlow, HltJetBuilder
pf = HltParticleFlow('pf')
pf.Inputs += [
    ['ProtoParticle',  'best',  'Rec/ProtoP/Charged'],
    ['ProtoParticle',  'gamma', 'Rec/ProtoP/Neutrals']]
if Mode == 'mue':
    pf.Inputs += [
        ['Particle', 'daughters', tesDimus.outputLocation()],
        ['Particle', 'daughters', seqDiels.outputLocation()],
        ['Particle', 'daughters', seqQdlps.outputLocation()]]
elif 'Turbo' in Type:
    pf.Inputs += [
        ['Particle', 'daughters', '/Event/Turbo/' + 
         loc + '/Particles'] for loc in TrgLocs]
elif Mode == 'mu':
    pf.Inputs += [
        ['Particle', 'daughters', seqDimus.outputLocation()],
        ['Particle', 'daughters', seqQdmus.outputLocation()]]
    if 'MC' in Type:
        pf.Inputs += [['Particle', 'daughters', seqSsmus.outputLocation()]]
elif Mode == 'Ks':
    pf.Inputs += [['Particle', 'daughters', seqKss.outputLocation()]]
elif Mode == 'pi':
    pf.Inputs += [
        ['Particle', 'daughters', seqDipss.outputLocation()],
        ['Particle', 'daughters', seqDipos.outputLocation()]]
pf.Output = 'Phys/PF/Particles'
pf.ProBestNames = ['mu+', 'e+', 'p+', 'K+', 'pi+']
pf.ProBestKeys  = [701,   700,  704,  703,  702]
pf.ProBestMins  = [0.5,   0.5,  0.5,  0.5,  0.5]
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
import GaudiPython, ROOT, IsoBdt
from Configurables import ApplicationMgr 
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

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
isoTool = IsoBdt.IsoBdt(dstTool)
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

# Simple ntuple class.
class Ntuple:
    def __init__(self, name, tes, genTool, rftTool, pvrTool, velTool, dstTool, 
                 isoTool, detTool, trkTool, l0Tool, hlt1Tool, hlt2Tool):
        from collections import OrderedDict
        import ROOT, array
        self.tes      = tes
        self.genTool  = genTool
        self.rftTool  = rftTool
        self.pvrTool  = pvrTool
        self.velTool  = velTool
        self.dstTool  = dstTool
        self.isoTool  = isoTool
        self.detTool  = detTool
        self.trkTool  = trkTool
        self.l0Tool   = l0Tool
        self.hlt1Tool = hlt1Tool
        self.hlt2Tool = hlt2Tool
        self.saved    = {}
        self.ntuple   = OrderedDict()
        self.tfile    = ROOT.TFile('output.root', 'RECREATE')
        self.ttree    = ROOT.TTree('data', 'data')
        vrsVrt = ['x', 'y', 'z', 'dx', 'dy', 'dz']
        vrsMom = ['pid', 'px', 'py', 'pz', 'e']
        if 'MC' in Type: vrsMom = ['idx_gen'] + vrsMom
        vrsPrt = ['iso', 'mu', 'loose_mu', 'tight_mu', 'pnn_mu', 'pnn_pi',
                  'pnn_k', 'pnn_p', 'pnn_ghost', 'prb_ghost', 'ip', 'ip_chi2',
                  'x0', 'y0', 'z0', 't0', 'p0', 'id0', 'z1', 'id1', 'id2',
                  'id3', 'xm2', 'ym2', 'zm2']
        vrsTag = ['m', 'dtf_m', 'dtf_dm', 'ip', 'ip_chi2', 'fd', 'fd_chi2', 
                  'doca', 'dtf_chi2', 'chi2', 've_ns', 'tt_ns', 'it_ns', 
                  'ot_ns', 'mu_ns']
        vrsTrg = ['l0_tos0', 'l0_tos1', 'l0_tis', 'hlt1_tos0', 'hlt1_tos1',
                  'hlt1_tos2', 'hlt1_tos3', 'hlt1_tis'] + [
            'hlt2_tos%i' % i for i in range(0, len(TrgLocs))] + ['hlt2_tis']
        if 'Turbo' in Type and Mode == 'pi':
            vrsPrt = ['mu', 'pnn_mu', 'prb_ghost', 'ip_chi2',
                      'xm2', 'ym2', 'zm2']
            vrsTag = ['m', 'ip_chi2', 'doca', 'chi2', 've_ns', 'mu_ns']
        self.vrsInit('pvr', vrsVrt)
        self.vrsInit('tag', ['idx_pvr', 'idx_prt0', 'idx_prt1'] +
                     vrsMom + vrsVrt + vrsTag + vrsTrg)
        if not 'Turbo' in Type: self.vrsInit('tag', [
                've_iso0', 've_iso1', 'ln_iso0', 'ln_iso1'])
        self.vrsInit('prt', ['idx_pvr'] + vrsMom + vrsPrt)
        if 'MC' in Type:
            vrsMcp = ['pid', 'q', 'px', 'py', 'pz', 'e', 'x', 'y', 'z']
            self.vrsInit('mcpvr', ['x', 'y', 'z'])
            self.vrsInit('mctag', ['idx_pvr', 'idx_prt0', 'idx_prt1', 
                                   'pid_mom'] + vrsMcp)
            self.vrsInit('mcprt', ['idx_pvr'] + vrsMcp)
        self.ntuple['pvr_n'] = array.array('d', [-1])
        self.ntuple['run_n'] = array.array('d', [-1])
        self.ntuple['evt_n'] = array.array('d', [-1])
        self.ntuple['evt_tck'] = array.array('d', [-1])
        self.ntuple['evt_spd'] = array.array('d', [-1])
        for key, val in self.ntuple.iteritems():
            if type(val) is array.array: self.ttree.Branch(key, val, key + '/D')
            else: self.ttree.Branch(key, val)
    def vrsInit(self, pre, vrs):
        import ROOT
        for v in vrs: self.ntuple['%s_%s' % (pre, v)] = ROOT.vector('double')()
    def keyPrt(self, obj):
        key = (obj.momentum().Px(), obj.momentum().Py(), obj.momentum().Pz())
        try: 
            trk = obj.proto().track()
            key = (trk.momentum().X(), trk.momentum().Y(), trk.momentum().Z())
        except:
            try:
                pos = obj.proto().calo()[0].position()
                key = (pos.x(), pos.y(), pos.z(), pos.e())
            except: pass
        return key
    def keyTrk(self, obj):
        key = float(1)
        try:
            trk = obj.proto().track()
            for lid in trk.lhcbIDs():
                if lid.isVelo():
                    key *= float(lid.veloID().channelID()/1000000.0)
        except: pass
        return key
    def keyVrt(self, obj):
        return (obj.position().X(), obj.position().Y(), obj.position().Z())
    def close(self):
        self.tfile.Write('data', ROOT.TObject.kOverwrite)
        self.tfile.Close()
    def clear(self):
        import array
        self.saved.clear()
        for key, val in self.ntuple.iteritems():
            if type(val) is array.array: val[0] = -1
            else: val.clear()
    def fill(self, key = None, val = None, idx = None):
        if key == None or val == None: self.tfile.Cd(''); self.ttree.Fill()
        elif key in self.ntuple: 
            if idx == None: self.ntuple[key].push_back(val)
            elif idx < len(self.ntuple[key]): self.ntuple[key][idx] = val
    def fillVrt(self, pre, vrt, cov = None, pos = None):
        from math import sqrt
        if cov == None: cov = vrt.covMatrix()
        if pos == None: pos = vrt.position() 
        self.fill('%s_x'  % pre, pos.X())
        self.fill('%s_y'  % pre, pos.Y())
        self.fill('%s_z'  % pre, pos.Z())
        self.fill('%s_dx' % pre, sqrt(abs(cov[0][0])))
        self.fill('%s_dy' % pre, sqrt(abs(cov[1][1])))
        self.fill('%s_dz' % pre, sqrt(abs(cov[2][2])))
    def fillMom(self, pre, mom):
        self.fill('%s_px' % pre, mom.Px())
        self.fill('%s_py' % pre, mom.Py())
        self.fill('%s_pz' % pre, mom.Pz())
        self.fill('%s_e'  % pre, mom.E())
        try:
            self.fill('%s_dtf_m'  % pre, mom.m().value())
            self.fill('%s_dtf_dm' % pre, mom.m().error())
        except:
            self.fill('%s_dtf_m'  % pre, -1)
            self.fill('%s_dtf_dm' % pre, -1)
    def fillPrt(self, prt, pvrs = None, prts = None):
        import ROOT
        pid = prt.particleID().pid()
        pro = prt.proto()
        try: prt = prt.data()
        except: pass
        try: pvr = self.pvrTool.relatedPV(prt, 'Rec/Vertex/Primary')
        except: pvr = None
        vrt = prt.endVertex()
        mom = prt.momentum()
        key = self.keyPrt(prt)
        pre = 'tag' if vrt else 'prt'
        if key in self.saved: return (pre, self.saved[key])
        # Daughters.
        if pre == 'tag':
            trks = []
            dtrs = prt.daughters()
            idxs = []
            for dtr in prt.daughters():
                (dtrPre, dtrIdx) = self.fillPrt(dtr)
                idxs += [dtrIdx]
                if dtr.proto() and dtr.proto().track():
                    trk, mu = dtr.proto().track(), None
                    hits = self.hits(trk)
	            if hits[7] == 0 and dtr.proto().muonPID():
 	                mu = dtr.proto().muonPID().muonTrack()
                    trks += [(trk, mu)]
            for prtIdx, dtr in enumerate(idxs):
                self.fill('%s_idx_prt%i' % (pre, prtIdx), dtr)
        idx = self.ntuple['%s_px' % pre].size()
        self.saved[key] = idx
        if pre == 'tag':
            # Shared hits.
            n = self.share(trks)
            # Vertex.
            from LoKiArrayFunctors.decorators import AMAXDOCA
            if pvr: 
                dtf = GaudiPython.gbl.DecayTreeFitter.Fitter(prt, pvr)
                dtf.fit()
            else: dtf = None
            if dtf and dtf.status() == 0:
                par = dtf.fitParams(prt)
                mom = par.momentum()
                self.fill('%s_dtf_chi2' % pre, dtf.chiSquare())
                self.fillVrt(pre, prt, par.posCovMatrix(), par.position())
            else:
                self.fill('%s_dtf_chi2' % pre, -1)
                self.fillVrt(pre, vrt)
            self.fill('%s_chi2' % pre, vrt.chi2())
            self.fill('%s_doca' % pre, AMAXDOCA('')(prt.daughters()))
            # Bs -> mu mu isolation.
            try: isos = self.isoTool.isolation(prt, pvr, pvrs, prts, 0)
            except: isos = [float('-inf'), float('-inf')]
            for dtrIdx, iso in enumerate(isos):
                self.fill('%s_ve_iso%i' % (pre, dtrIdx), iso)
            try: isos = self.isoTool.isolation(prt, pvr, pvrs, prts, 1)
            except: isos = [float('-inf'), float('-inf')]
            for dtrIdx, iso in enumerate(isos):
                self.fill('%s_ln_iso%i' % (pre, dtrIdx), iso)
        # Momentum and mass.
        self.fill('%s_m'  % pre, prt.measuredMass())
        self.fillMom(pre, mom)
        # Trigger.
        self.l0Tool.setOfflineInput(prt)
        self.hlt1Tool.setOfflineInput(prt)
        self.hlt2Tool.setOfflineInput(prt)
        self.l0Tool.setTriggerInput('L0DiMuonDecision')
        trg = l0Tool.tisTosTobTrigger()
        self.fill('%s_l0_tos0' % pre, trg.tos())
        self.l0Tool.setTriggerInput('L0MuonDecision')
        trg = l0Tool.tisTosTobTrigger()
        self.fill('%s_l0_tos1' % pre, trg.tos())
        self.l0Tool.setTriggerInput('L0HadronDecision')
        trg = self.l0Tool.tisTosTobTrigger()
        self.fill('%s_l0_tis' % pre, trg.tis())
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
        self.hlt2Tool.setTriggerInput('Hlt2Topo.*')
        trg = self.hlt2Tool.tisTosTobTrigger()
        self.fill('%s_hlt2_tis' % pre, trg.tis())
        for locIdx, loc in enumerate(TrgLocs):
            self.hlt2Tool.setTriggerInput(loc + 'Decision')
            trg = self.hlt2Tool.tisTosTobTrigger()
            self.fill('%s_hlt2_tos%i' % (pre, locIdx), trg.tos())
        # Particle ID.
        self.fill('%s_pid' % pre, pid)
        if pro:
            try:
                self.fill('%s_mu' % pre, pro.muonPID().IsMuon())
                self.fill('%s_loose_mu' % pre, pro.muonPID().IsMuonLoose())
                self.fill('%s_tight_mu' % pre, pro.muonPID().IsMuonTight())
            except: pass
            self.fill('%s_pnn_mu'    % pre, pro.info(701, -100))
            self.fill('%s_pnn_pi'    % pre, pro.info(702, -100))
            self.fill('%s_pnn_k'     % pre, pro.info(703, -100))
            self.fill('%s_pnn_p'     % pre, pro.info(704, -100))
            self.fill('%s_pnn_ghost' % pre, pro.info(705, -100))
            self.fill('%s_prb_ghost' % pre, pro.track().ghostProbability()
                      if pro.track() else -100)
            # Track.
            trk = pro.track()
            ids = []
            for i in trk.lhcbIDs():
                if i.isVelo():
                    d = self.detTool.sensor(i.veloID())
                    ids += [(d.z(), d, i)]
                ids.sort()
            for hit in range(0, 4):
                if hit >= len(ids):
                    self.fill('%s_x%i'  % (pre, hit), -1)
                    self.fill('%s_y%i'  % (pre, hit), -1)
                    self.fill('%s_z%i'  % (pre, hit), -1)
                    self.fill('%s_id%i' % (pre, hit), -1)
                    self.fill('%s_t%i' % (pre, hit),  -1)
                    self.fill('%s_p%i' % (pre, hit),  -1)
                    continue
                z, d, i = ids[hit]
                s = i.veloID()
                v = GaudiPython.gbl.LHCb.StateVector()
                self.trkTool.propagate(trk, z, v, prt.particleID())
                self.fill('%s_x%i'  % (pre, hit), v.x())
                self.fill('%s_y%i'  % (pre, hit), v.y())
                self.fill('%s_z%i'  % (pre, hit), v.z())
                self.fill('%s_id%i' % (pre, hit), i.lhcbID()*
                          (-1 if s.isPhiType() else 1))
                if (s.isPhiType()):
                    self.fill('%s_t%i' % (pre, hit), d.globalPhi(s.strip(), 0))
                    self.fill('%s_p%i' % (pre, hit), d.phiPitch(s.strip()))
                if (s.isRType()):
                    self.fill('%s_t%i' % (pre, hit), d.globalR(s.strip(), 0))
                    self.fill('%s_p%i' % (pre, hit), d.rPitch(s.strip()))
            # Muon extrapolation.
            z = 15270.0
            v = GaudiPython.gbl.LHCb.StateVector()
            self.trkTool.propagate(trk, z, v, prt.particleID())
            self.fill('%s_xm2' % pre, v.x())
            self.fill('%s_ym2' % pre, v.y())
            self.fill('%s_zm2' % pre, v.z())
        # Linked MC particle.
        if 'MC' in Type:
            gen = None; wgt = 0; rels = self.genTool.relatedMCPs(prt)
            for rel in rels: gen = rel.to() if rel.weight() > wgt else gen
            if gen: (genPre, genIdx) = self.fillMcp(gen) 
            else: genIdx = -1
            self.fill('%s_idx_gen' % pre, genIdx)
        # IP.
        ip, ipChi2 = ROOT.Double(-1), ROOT.Double(-1)
        if pvr: self.dstTool.distance(prt, pvr, ip, ipChi2)
        self.fill('%s_ip' % pre, ip)
        self.fill('%s_ip_chi2' % pre, ipChi2)
        fd, fdChi2 = ROOT.Double(-1), ROOT.Double(-1)
        if pvr and vrt: self.dstTool.distance(vrt, pvr, fd, fdChi2)
        self.fill('%s_fd' % pre, fd)
        self.fill('%s_fd_chi2' % pre, fdChi2)
        # Isolation.
        if pre == 'prt':
            jets = self.tes[jb.Output]
            iso  = None
            vid  = self.keyTrk(prt)
            for jet in jets:
                for dtr in jet.daughters():
                    if self.keyTrk(dtr) == vid: iso = jet; break
                if iso: break
            self.fill('%s_iso' % pre, float(prt.pt()/iso.pt()) if iso else -1)
        # Primary vertex.
        if pvr:
            key = self.keyVrt(pvr)
            if not key in self.saved:
                self.saved[key] = self.ntuple['pvr_x'].size()
                self.fillVrt('pvr', pvr)
            self.fill('%s_idx_pvr' % pre, self.saved[key])
        else: self.fill('%s_idx_pvr' % pre, -1)
        return (pre, idx)
    def deltaR(self, p1, p2):
        import ROOT
        v1 = ROOT.TVector3(p1.momentum().X(), p1.momentum().Y(), 
                           p1.momentum().Z())
        v2 = ROOT.TVector3(p2.momentum().X(), p2.momentum().Y(), 
                           p2.momentum().Z())
        return v1.DeltaR(v2)
    def fillMcp(self, prt):
        import ROOT
        pid = prt.particleID().pid()
        mom = prt.momentum()
        pos = None
        key = self.keyPrt(prt)
        pre = 'mctag' if pid in [23, 32, 36] else 'mcprt'
        if key in self.saved: return (pre, self.saved[key])
        # Daughters.
        if pre == 'mctag':
            idxs = []
            for vrt in prt.endVertices():
                for dtr in vrt.products():
                    if abs(dtr.particleID().pid()) != 13: continue
                    pos = vrt.position()
                    (dtrPre, dtrIdx) = self.fillMcp(dtr)
                    idxs += [dtrIdx]
            for prtIdx, dtr in enumerate(idxs):
                self.fill('%s_idx_prt%i' % (pre, prtIdx), dtr)
        idx = self.ntuple['%s_px' % pre].size()
        self.saved[key] = idx
        # Momentum.
        self.fill('%s_px' % pre, mom.Px())
        self.fill('%s_py' % pre, mom.Py())
        self.fill('%s_pz' % pre, mom.Pz())
        self.fill('%s_e'  % pre, mom.E())
        # PID.
        self.fill('%s_q' % pre, float(prt.particleID().threeCharge()) / 3.0)
        self.fill('%s_pid' % pre, pid)
        try: self.fill('%s_pid_mom' % pre, prt.originVertex().mother().
                       particleID().pid())
        except: self.fill('%s_pid_mom' % pre, 0)
        # Vertex.
        if not pos: pos = prt.originVertex().position()
        self.fill('%s_x' % pre, pos.X())
        self.fill('%s_y' % pre, pos.Y())
        self.fill('%s_z' % pre, pos.Z())
        # Primary vertex.
        pvr = prt.primaryVertex()
        if pvr:
            key = self.keyVrt(pvr)
            if not key in self.saved:
                self.saved[key] = self.ntuple['mcpvr_x'].size()
                self.fill('mcpvr_x', pvr.position().X())
                self.fill('mcpvr_y', pvr.position().Y())
                self.fill('mcpvr_z', pvr.position().Z())
            self.fill('%s_idx_pvr' % pre, self.saved[key])
        else: self.fill('%s_idx_pvr' % pre, -1)
        return (pre, idx)
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
    def share(self, trks):
        ns = {n:0 for n in range(0, 13)}
        for trk1 in trks:
            for trk2 in trks:
                if trk1[0] == trk2[0]: continue
                hits = self.hits(trk1[0], trk2[0], trk1[1], trk2[1])
                for key, val in hits.iteritems():
                    if val > ns[key]: ns[key] = val
        for cat in TrkCats:
            n = ns[cat[1]]
            self.fill('tag_%s_ns' % (cat[0]), n)
        return ns[TrkCats[0][1]]

# Run.
import sys, ROOT
try: evtmax = int(sys.argv[1])
except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, genTool, rftTool, pvrTool, None, dstTool,
                isoTool, None, trkTool, l0Tool, hlt1Tool, hlt2Tool)
while evtnum < evtmax:
    gaudi.run(1)
    if not bool(tes['/Event']): break
    evtnum += 1
    ntuple.clear()

    # Fill event info.
    daq = tes['DAQ/ODIN']
    try:
        ntuple.ntuple['run_n'][0] = daq.runNumber()
        ntuple.ntuple['evt_n'][0] = daq.eventNumber()
        ntuple.ntuple['evt_tck'][0] = daq.triggerConfigurationKey()
    except: continue
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
    if 'Turbo' in Type:
        if Mode == 'pi':
            for loc in [seqDipss.outputLocation(), seqDipos.outputLocation()]:
                prts = tes[loc]
                try: len(prts); run = True
                except: run = False
                if run:
                    for prt in prts: 
                        if prt.measuredMass() < 1300:
                            ntuple.fillPrt(prt); fill = True
        elif Mode == 'mue':
            prts = tes[seqQdlps.outputLocation()]
            try: len(prts); run = True
            except: run = False
            if run:
                for prt in prts: 
                    (prtPre, prtIdx) = ntuple.fillPrt(prt)
                    fill = True
        else:
            for locIdx, loc in enumerate(TrgLocs):
                prts = tes['/Event/Turbo/' + loc + '/Particles']
                try: len(prts); run = True
                except: run = False
                if run:
                    for prt in prts: 
                        (prtPre, prtIdx) = ntuple.fillPrt(prt)
                        fill = True
                        if prtPre == 'tag': ntuple.fill('%s_hlt2_tos%i' % (
                                prtPre, locIdx), 1, prtIdx)
    elif Mode == 'mu':
        if 'MC' in Type:
            prts = tes[seqSsmus.outputLocation()]
            try: len(prts); run = True
            except: run = False
            if run:
                for prt in prts: ntuple.fillPrt(prt, pvrs, trks); fill = True
        pvrs = tes['Rec/Vertex/Primary']
        trks = tes['Rec/Track/Best']
        prts = tes[seqDimus.outputLocation()]
        sigs = []
        try: len(prts); run = True
        except: run = False
        if run:
            for prt in prts:
                if (prt.measuredMass() < 1500 or 'CL' not in Type
                    or rndTool.Rndm() < 0.1):
                    sigs += [prt]; ntuple.fillPrt(prt, pvrs, trks); fill = True
        prts = tes[seqQdmus.outputLocation()]
        try: len(prts); run = True
        except: run = False
        if run:
            for prt in prts:
                for dtr in prt.daughters(): ntuple.fillPrt(dtr, pvrs, trks)
                fill = True
        prts = tes[pf.Output]
        try: len(prts); run = len(sigs) > 0
        except: run = False
        if run:
            for prt in prts:
                if abs(prt.charge()) != 1: continue
                try: pvr = pvrTool.relatedPV(prt, 'Rec/Vertex/Primary')
                except: pvr = None
                if not pvr: continue
                ip, ipChi2 = ROOT.Double(-1), ROOT.Double(-1)
                dstTool.distance(prt, pvr, ip, ipChi2)
                if ipChi2 > 9 and prt.momentum().Pt() > 200:
                    minDoca = 100
                    for sig in sigs:
                        for dtr in sig.daughters():
                            doca = docaTool.doca(prt, dtr)
                            if doca < minDoca: minDoca = doca
                    if minDoca < 0.2: ntuple.fillPrt(prt); fill = True
    elif Mode == 'Ks':
        prts = tes[seqKss.outputLocation()]
        try: len(prts); run = True
        except: run = False
        if run:
            for prt in prts: ntuple.fillPrt(prt); fill = True
    elif Mode == 'pi':
        for loc in [seqDipss.outputLocation(), seqDipos.outputLocation()]:
            prts = tes[loc]
            try: len(prts); run = True
            except: run = False
            if run:
                for prt in prts: ntuple.fillPrt(prt); fill = True
    if 'MC' in Type:
        prts = tes['Phys/StdAllLooseMuons/Particles']
        try: len(prts); run = True
        except: run = False
        if run:
            for prt in prts: ntuple.fillPrt(prt); fill = True

    # Fill MC.
    mcps = tes['MC/Particles']
    try: len(mcps); run = True
    except: run = False
    if run:
        for mcp in mcps:
            pid = mcp.particleID()
            if (abs(pid.pid()) in [13, 23, 32, 36] or pid.hasCharm() or 
                pid.hasBottom()): ntuple.fillMcp(mcp); fill = True
    if fill: ntuple.fill()

# Close and write the output.
ntuple.close()
