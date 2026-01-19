###############################################################################
# Script to analyze check for correct filling of gen-level information.       #
# Author: Michael Peters                                                      #
###############################################################################

import ROOT

infile = 'eta2MuMuGamma_mc_20251121.root'
print(f'Reading from {infile}.')

tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# =============================================================================


def get_decays(gentag_pid, genprt_pid):
    """
    Groups gentag_pid and genprt_pid by event.
    Returns a list of (gentag_pids, genprt_pids) for each decay.
    """
    decays = []  # Each element: (idx_mom, gentag_pids, genprt_pids)
    idx_mom = None
    gentag_pids, genprt_pids = [], []
    for i, idx in enumerate(genprt_idx_mom):
        idx = int(idx)
        # For each event, categorize by mother index, each dtr particle is
        # at unique index i.
        # Set initial mother index
        if idx_mom is None: idx_mom = idx
        # If new mother index, (1) store previous decay info, (2) set new
        # idx, (3) reset lists
        if idx != idx_mom:
            decays.append((idx_mom, gentag_pids, genprt_pids))
            idx_mom = idx
            gentag_pids, genprt_pids = [], []
        # Otherwise, add gentag_pid and genprt_pid to current mother index lists
        if i < len(gentag_pid):
            gentag_pids.append(int(gentag_pid[i]))
        if i < len(genprt_pid):
            genprt_pids.append(int(genprt_pid[i]))
    # Add last decay not added in loop
    if idx_mom is not None:
        decays.append((idx_mom, gentag_pids, genprt_pids))
    return decays


# =============================================================================

nErrors, nMuErrors, nPhoErrors = 0, 0, 0
missingPrts = 0
for entryIdx in range(tree.GetEntries()):
    tree.GetEntry(entryIdx)

    gentag_pid = getattr(tree, 'gentag_pid')
    genprt_pid = getattr(tree, 'genprt_pid')
    genprt_idx_mom = getattr(tree, 'genprt_idx_mom')

    if len(gentag_pid) == 0: continue  # no gen-level tags, skip
    # print(f'Entry {entryIdx}, gentag_pid: {gentag_pid}, ' +
    #       f'genprt_pid: {genprt_pid}, genprt_idx_mom: {genprt_idx_mom}')  # debug 

    # Count number of each particle type in gentag_pid and genprt_pid
    nEta, nMup, nMum, nPho = 0, 0, 0, 0
    for pid in gentag_pid:
        if pid == 221: nEta += 1
    for pid in genprt_pid:
        if pid == -13: nMum += 1
        elif pid == 13: nMup += 1
        elif pid == 22: nPho += 1
    
    # There should be one eta for each muon pair + photon
    if nEta != nMup or nEta != nMum or nEta != nPho:
        nErrors += 1
        print(f'Error at entry {entryIdx}: nEta {nEta}, nMum {nMum}, ' +
              f'nMup {nMup}, nPho {nPho}')
        print(f'  - gentag_pid: {gentag_pid}')
        print(f'  - genprt_pid: {genprt_pid}')
        # If there are no muons, then the photon should also be missing
        if nEta != nMup or nEta != nMum:
            nPhoErrors += 1
            print('  - Photon error.')
        # If there is no photon, then the muon pair should also be missing
        if nEta != nPho:
            nMuErrors += 1
            print('  - Muon error.')
        missingPrts += (len(gentag_pid) * 3) - len(genprt_pid)

print('========================================')
print('Gen-level info check results:')
print(f'Total errors: {nErrors}')
print(f'  - Photon errors: {nPhoErrors}')
print(f'  - Muon pair errors: {nMuErrors}')
print(f'Missing particles: {missingPrts}')