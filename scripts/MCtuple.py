
# Find every MC generator-level eta -> mu mu (gamma) and save it in ntuple.
def findMCs(self):
  try: mcs = self.mcParticles(); print("Found mcParticles()")
  except:
    try: mcprts = self.tes['/Event/MC/Particles']
    except: pass; print("Failed to MCParticles container")

  # for mcprt in mcprts:
  # if !221, continue
  # (if 221) find final state for mcpart
  # mup, mum, pho = 0
  # if -13->mups++, 13->mums++, 22->phos++
  # if mum,mup,phos=1,1,1, fillMcp(self,eta)->into new ntuple branch

# ---------------------------------------------------------------------------

## !! TODO: rewrite to fill into new ntuple branch "gen-level"
def fillMcp(self, prt):
  pid = prt.particleID().pid()
  mom = prt.momentum()
  pos = None
  key = self.key(prt)
  pre = 'mctag' if pid in [221] else 'mcprt'
  # Prevent filling same particle more than once.
  if key in self.saved: return (pre, self.saved[key])

  # Daughters.
  if pre == 'mctag':
    idxs = []
    # Loop over decay points (end vertices)
    for vrt in prt.endVertices():
      # Loop over final products (daughters)
      for dtr in vrt.products():
        # 13 == muon, 22 = photon
        # If daughters (of mc eta) are muon or photon, save particle
        if abs(dtr.particleID().pid()) not in [13, 22]: continue
        pos = vrt.position()
        (dtrPre,dtrIdx) = self.fillMcp(dtr) # Recursively loop thru daughters
        idxs += [dtrIdx]
    ## !! Index daughters by linking candiate to prtIdx-th daughter.
    ## Example: 'tag_idx_prt0 = 42' means the first daughter of the candidate
    ## is entry 42 in the tree. ??
    # Build truth-level decay chain
    for prtIdx, dtr in enumerate(idxs):
      self.fill('%s_idx_prt%i' % (pre, prtIdx), dtr)

  # Save current index of TTree array, get next row index of TTree.
  idx = self.ntuple['%s_px' % pre].size()
  self.saved[key] = idx

  # Momentum.
  self.fillMom(pre, mom)

  # PID.
  ### ?? WHAT IS THE PURPOSE OF THIS LINE ??
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
