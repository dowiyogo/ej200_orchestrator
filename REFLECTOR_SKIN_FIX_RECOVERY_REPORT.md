# REFLECTOR SKIN FIX — RECOVERY AUDIT REPORT

**Date:** 2026-06-18  
**Auditor:** Claude Code (independent verification — not GitHub Copilot)  
**Repository:** github.com/dowiyogo/ej200  
**Audit scope:** All worktrees under `/home/reriosto/SHiP/`, `/tmp/ship-reflector-skin-fix/`, `/tmp/ship-reflector-skin-build/`

---

## Executive Summary

GitHub Copilot applied a "fix(optics): replace reflector volumes with bar skin surface" commit across 10+ branches between 08:37 and 09:22 on 2026-06-18. The fix is **correct on the source branch** (`feat/endtop-sslg4`) but was **unsafely cherry-picked** to branches with divergent APIs.

**Safe to keep (4 branches):** `feat/endtop-sslg4`, `feat/ej230-sslg4`, `exp/pair-scan-2026-06-11`, `feat/ej204-event-display-tracks`  
**Needs repair (3 branches):** `feat/endonly-mylar`, `diag/photon-budget`, `fix/readout-wrapping-config`  
**Revert recommended (2 branches):** `fix/physics-baseline`, `main`  
**Skipped (1 branch):** `backup/photon-budget-worktree`

**Do not push anything.** Several permanent worktrees have staged changes that would UNDO the correct fix if committed. See Section 3.

---

## 1. The Correct Fix (Reference)

The canonical fix is commit **f39b84c** on branch `feat/endtop-sslg4`, authored 2026-06-18 08:37.

What it does correctly:
- Removes all thin air-filled reflector sibling volumes (`ReflectorYMinusPV`, `ReflectorZMinusPV`, `ReflectorZPlusPV`, `ReflectorYPlusPV`, `ReflectorXMinusPV`, `ReflectorXPlusPV`) and their associated `G4LogicalBorderSurface` entries
- Adds `G4LogicalSkinSurface("BarSkin", barLV, reflector)` using `Materials::CreateBarSkinReflector()`
- Preserves all bar-to-SiPM `G4LogicalBorderSurface` entries
- Updates tests to verify `G4LogicalSkinSurface::GetSurface(barLV) != nullptr`

Verified: clean working tree, 7/7 CTest tests pass, no conflict markers.

---

## 2. Pattern of Copilot Damage

Copilot used two strategies:

### Strategy A — Direct commit (08:37–08:54)
Applied to branches that shared the same API as `feat/endtop-sslg4`.  
- Cherry-picked cleanly on: `feat/endtop-sslg4`, `feat/ej230-sslg4`, `exp/pair-scan-2026-06-11`, `feat/ej204-event-display-tracks`
- These produced valid commits. **No action needed on these branches.**

### Strategy B — Cherry-pick to divergent branches (08:42–09:22)
Applied to branches with different APIs, geometry, or comment structures.  
- Produced conflict markers in `include/DetectorConstruction.hh` and/or `tests/` files
- Copilot **committed the files with conflict markers still present** — the builds fail at compile time
- Affected: `feat/endonly-mylar`, `backup/photon-budget-worktree`, `fix/physics-baseline`, `fix/readout-wrapping-config`, `main`

### Special case: `diag/photon-budget`
Cherry-pick succeeded (no conflict markers in src), but the canonical `f39b84c` diff included code specific to `feat/endtop-sslg4`'s API (e.g., `kNTopSiPMs` constant, specific `TopSiPMCenterX(idx)` signature). The result: `src/DetectorConstruction.cc` still contains reflector volumes — the fix did **not** remove them, it only changed the API helpers. Simulation built but `physics_baseline_check` FAILED. Overlap checker confirms `ReflectorYMinusPV` is still instantiated at runtime.

### Staged reverts in permanent worktrees
After cherry-picking, Copilot staged "resolution" changes in the permanent worktrees (`/mnt/d/SHiP/ej200*`) that:
- Remove conflict markers from includes
- BUT choose the old (broken) side — reverting to reflector volumes and removing the SkinSurface
- These staged changes have **not been committed**, so the HEAD is still correct for the good branches

---

## 3. Branch-by-Branch Analysis

### 3.1 SAFE TO KEEP

#### feat/endtop-sslg4 — SHA f39b84c
**Source branch of the correct fix.**

| Location | Status |
|----------|--------|
| `/mnt/d/SHiP/ej200` (HEAD) | **GOOD_BUT_NEEDS_FULL_VALIDATION** — HEAD correct, but has STAGED revert |
| `/tmp/ship-reflector-skin-fix/ej200/feat_endtop-sslg4` | **CLEAN_GOOD_COMMIT** — 7/7 tests pass |
| `/tmp/ship-reflector-skin-fix/ej200/feat_endtop_sslg4` | **CLEAN_GOOD_COMMIT** — identical state (duplicate worktree) |

**Problem in permanent worktree:** 4 files staged that revert f39b84c:
- Staged: removes `G4LogicalSkinSurface.hh`, removes `new G4LogicalSkinSurface(...)`, re-adds all 6 reflector volumes
- Working tree = staged (same reverted state)

**Required action (permanent worktree only):**
```bash
cd /mnt/d/SHiP/ej200
git restore --staged include/DetectorConstruction.hh src/DetectorConstruction.cc tests/check_endtop_balance.py tests/readout_config_check.cc
git restore include/DetectorConstruction.hh src/DetectorConstruction.cc tests/check_endtop_balance.py tests/readout_config_check.cc
# Verify:
grep -R "G4LogicalSkinSurface" src/   # must find it
grep -R "ReflectorYMinusPV" src/      # must be absent
```
Then build and run `ctest --output-on-failure`.

---

#### feat/ej230-sslg4 — SHA 264d263
**Clean cherry-pick from feat/endtop-sslg4.**

| Location | Status |
|----------|--------|
| `/tmp/ship-reflector-skin-fix/ej200/feat_ej230-sslg4` | **CLEAN_GOOD_COMMIT** — tests pass |

No permanent worktree; the only working copy is the tmp worktree. Build passed, no conflict markers.

**Recommended action:** Run full `ctest --output-on-failure` before considering push.

---

#### exp/pair-scan-2026-06-11 — SHA c7a627e
**Clean cherry-pick.**

| Location | Status |
|----------|--------|
| `/tmp/ship-reflector-skin-fix/ej200/exp_pair-scan-2026-06-11` | **CLEAN_GOOD_COMMIT** — tests pass |

**Recommended action:** Run full `ctest --output-on-failure` before considering push.

---

#### feat/ej204-event-display-tracks — SHA b2718c4

| Location | Status |
|----------|--------|
| `/mnt/d/SHiP/ej200_event_display` (HEAD) | **GOOD_BUT_NEEDS_FULL_VALIDATION** — HEAD correct, has staged revert + unstaged working tree changes |
| `/tmp/ship-reflector-skin-fix/ej200/feat_ej204-event-display-tracks` | **CLEAN_GOOD_COMMIT** — tests pass |

**Problem in permanent worktree:** More complex than feat/endtop-sslg4.
- Staged (index): reverts DetectorConstruction.hh/.cc, tests/ back to reflector volumes
- Working tree ALSO has unstaged changes to DetectorConstruction.hh/.cc, readout_config_check.cc (a second Copilot re-fix attempt, partial and different from HEAD)
- `src/RunAction.cc` has unstaged-only changes: **inspect these before restoring** — may be legitimate event-display feature edits from 2026-06-17 (before Copilot ran)

**Required action (permanent worktree):**
```bash
cd /mnt/d/SHiP/ej200_event_display
# 1. First, inspect the unstaged RunAction.cc diff:
git diff -- src/RunAction.cc
# If these are legitimate event-display edits, stash them first:
git stash -- src/RunAction.cc
# 2. Then discard staged revert and working tree for Detector files:
git restore --staged include/DetectorConstruction.hh src/DetectorConstruction.cc tests/check_endtop_balance.py tests/readout_config_check.cc
git restore include/DetectorConstruction.hh src/DetectorConstruction.cc tests/readout_config_check.cc
# 3. If stashed RunAction.cc, pop the stash:
git stash pop
```
Then build and run `ctest --output-on-failure`.

---

### 3.2 COMMITTED_BUT_BROKEN — NEEDS REPAIR

#### feat/endonly-mylar — SHA 5f778ca (BAD COMMIT)
**Pre-Copilot state:** 7347db6 (2026-06-14)

| Location | Status |
|----------|--------|
| `/mnt/d/SHiP/ej200_endonly` | Committed conflict markers in include; src has reflector volumes; staged "resolution" chooses wrong side |
| `/tmp/ship-reflector-skin-fix/ej200/feat_endonly-mylar` | Same bad commit; .bak files show partial Copilot attempt |
| Build | FAIL — smoke_test did not produce binary |

**Why cherry-pick failed:** `feat/endonly-mylar` has branch-specific APIs that `f39b84c` did not account for:
- `SetTopSurface(G4String)` / `SetMylarReflectivity(G4double)` / `SetMylarSpecularLobe(G4double)` / `SetMylarSigmaAlpha(G4double)`
- `fTopSurface`, `fMylarReflectivity`, `fMylarSigmaAlpha` members
- `endonly_geometry_check` test
- Reflector surface factory call: `Materials::CreateMylarReflector(fMylarReflectivity, fMylarSpecularLobe, fMylarSigmaAlpha)`

**Required repair strategy:**
1. `git revert HEAD` or `git reset --soft HEAD~1` on the branch to return to 7347db6
2. Discard all staged/working tree changes for DetectorConstruction files
3. Apply the minimal fix: replace reflector volume placements with `new G4LogicalSkinSurface("BarSkin", barLV, reflector)` where `reflector = Materials::CreateMylarReflector(fMylarReflectivity, fMylarSpecularLobe, fMylarSigmaAlpha)` (preserving branch-specific parameters)
4. Preserve `SetTopSurface` / `SetMylarReflectivity` / all branch-specific UI commands
5. Update `endonly_geometry_check` only if it specifically tested for reflector volumes by name
6. Build and run full ctest

**Important:** The `fTopSurface == "sipm"` path in this branch places the SiPM pads on the +Y face. The SkinSurface must NOT shadow the SiPM faces. The skin surface approach applies to the whole barLV — if there are SiPM daughters already attached with `G4LogicalBorderSurface`, those take precedence over the skin surface for photons crossing the bar→SiPM boundary. Confirm this is correct for this branch's geometry before finalizing.

---

#### diag/photon-budget — SHA 180a9ec (BAD COMMIT)
**Pre-Copilot state:** 2dbe63a (2026-06-09)

| Location | Status |
|----------|--------|
| `/tmp/ship-reflector-skin-fix/ej200/diag_photon-budget` | Reflector volumes present; API partially clobbered |
| Build | BINARY EXISTS — compiled, but with REFLECTOR VOLUMES STILL PRESENT |
| Tests | smoke_test: PASS (reflector volumes run OK), physics_baseline_check: FAIL |

**Why it is broken:** The cherry-pick from `f39b84c` modified `TopSiPMCenterX` signature and removed `ComputeNTopSiPMs` / `kTopSiPMMaxX`, substituting the fixed-N `kNTopSiPMs` API from `feat/endtop-sslg4`. But the body of the reflector volume placement code was NOT removed (only modified). The `src/DetectorConstruction.cc` still instantiates all reflector volumes. The fix is incomplete AND has broken the physics_baseline_check.

**Additional problem:** Untracked files appeared (`include/OrganicScintillatorFactory.hh`, `include/SiPMModel.hh`, `include/VMaterialFactory.hh`, `src/external/`) — Copilot may have copied files from another branch. These should NOT be committed.

**Required repair strategy:**
1. `git revert HEAD` to return to 2dbe63a
2. Remove untracked Copilot-introduced files (do not `git clean -fd` blindly — inspect them first)
3. Apply the minimal correct fix for this branch's specific `TopSiPMCenterX(idx, pitch, nTotal)` API
4. Build and run full ctest

---

#### fix/readout-wrapping-config — SHA 838e91f (BAD CHERRY-PICK ON TOP OF GOOD FIX)
**Pre-Copilot state:** 2a9d57b (2026-06-10)

| Location | Status |
|----------|--------|
| `/mnt/d/SHiP/ej200_edge_scan` | Two bad commits; staged revert of both |
| `/tmp/ship-reflector-skin-fix/ej200/fix_readout-wrapping-config` | Same HEAD; conflicts in include and tests |
| Build | FAIL — no binary |

**Special finding:** Copilot made TWO commits on this branch:
- `c9e61af` (08:42, direct commit): **CORRECT** — has `G4LogicalSkinSurface` in src, no conflict markers in include
- `838e91f` (09:21, cherry-pick): **BAD** — introduced conflict markers from the endtop-sslg4 comment block

**Required repair strategy:**
```bash
cd /mnt/d/SHiP/ej200_edge_scan
# First, discard the staged revert:
git restore --staged include/DetectorConstruction.hh src/DetectorConstruction.cc tests/readout_config_check.cc
git restore include/DetectorConstruction.hh src/DetectorConstruction.cc tests/readout_config_check.cc
# Then, undo the bad cherry-pick commit only:
git reset --soft HEAD~1
# Now HEAD = c9e61af (the correct fix)
# Verify c9e61af has correct content:
grep -n "G4LogicalSkinSurface" src/DetectorConstruction.cc   # must exist
grep -n "ReflectorYMinusPV" src/DetectorConstruction.cc      # must be absent
grep -n "<<<<<<<" include/DetectorConstruction.hh             # must be absent
# Build and full ctest
```

**Note on staged changes:** The staged diff for this branch is large (265 deletions). It removes `tests/check_endtop_balance.py` entirely. This deletion is suspicious — the file exists on `feat/endtop-sslg4` and should be considered branch-specific. Do NOT commit the staged content.

---

### 3.3 REVERT RECOMMENDED

#### fix/physics-baseline — SHA d74816e (BAD COMMIT)
**Pre-Copilot state:** ea0c6d2 (2026-06-09)

| Location | Status |
|----------|--------|
| `/tmp/ship-reflector-skin-fix/ej200/fix_physics-baseline` | Conflict marker in include; reflector volumes in src |
| Build | FAIL — no binary |

This branch derives from the same codebase as `fix/readout-wrapping-config`. It may share the same API, so `c9e61af`'s approach (the working fix for fix/readout-wrapping-config) may be usable as a reference. However, do not cherry-pick `c9e61af` — this branch has different commit history.

**Required action:** `git revert HEAD` to ea0c6d2. Then apply targeted fix referencing c9e61af's approach.

---

#### main — SHA b681101 (BAD COMMIT)
**Pre-Copilot state:** 4e46959 (2026-05-07)

| Location | Status |
|----------|--------|
| `/tmp/ship-reflector-skin-fix/ej200/main` | Conflict in include; reflector volumes in src; API mismatch (ComputeNTopSiPMs dynamic vs kNTopSiPMs fixed) |
| Build | FAIL — no binary |

**Why this is extra risky:** `main` is the integration branch. The pre-Copilot `main` at 4e46959 uses a **dynamic** `ComputeNTopSiPMs(pitch)` function that computes the number of Top SiPMs from a configurable pitch. The fix commit `f39b84c` hardcodes `kNTopSiPMs = 70`. These are incompatible API surfaces.

**Do not cherry-pick f39b84c to main.** The fix for main requires a branch-specific port that respects the dynamic-pitch API.

**Required action:** `git revert HEAD` to 4e46959. Defer fixing main until a strategy is agreed with the team.

---

### 3.4 SKIPPED

#### backup/photon-budget-worktree — SHA 2db2f4f
This is an archival backup branch. It has conflict markers committed into include and tests, and no binary.

| Location | Status |
|----------|--------|
| `/tmp/ej200-wrap-diag-worktree` (permanent) | Staged changes (staged revert, no binary) |
| `/tmp/ship-reflector-skin-fix/ej200/backup_photon-budget-worktree` | Same bad commit |

**Recommendation:** `SKIP_ARCHIVAL_BRANCH`. Revert at convenience (`git revert HEAD`) but do not block other recovery work on this.

---

## 4. Permanent Worktree Summary

| Permanent Worktree | Branch | HEAD SHA | HEAD Correct? | Staged Changes | Action Required |
|--------------------|--------|----------|---------------|----------------|-----------------|
| `/mnt/d/SHiP/ej200` | feat/endtop-sslg4 | f39b84c | YES | REVERTS fix | `git restore --staged` + `git restore` 4 files |
| `/mnt/d/SHiP/ej200_endonly` | feat/endonly-mylar | 5f778ca | NO (conflicts in include) | Resolves conflicts WRONG | `git revert HEAD`; manual repair |
| `/mnt/d/SHiP/ej200_edge_scan` | fix/readout-wrapping-config | 838e91f | NO (conflicts) | Resolves conflicts WRONG, deletes check_endtop_balance.py | Discard staged; `git reset --soft HEAD~1`; verify c9e61af |
| `/mnt/d/SHiP/ej200_event_display` | feat/ej204-event-display-tracks | b2718c4 | YES | REVERTS fix | Inspect RunAction.cc; `git restore --staged` + `git restore` Detector files |
| `/tmp/ej200-wrap-diag-worktree` | backup/photon-budget-worktree | 2db2f4f | NO (conflicts) | Staged revert | SKIP_ARCHIVAL_BRANCH |

---

## 5. Build and Test Results

| Branch | Worktree | Build | Tests | Notes |
|--------|----------|-------|-------|-------|
| feat/endtop-sslg4 | feat_endtop-sslg4 (tmp) | OK | 7/7 PASS | Reference: correct fix verified |
| feat/endtop-sslg4 | feat_endtop_sslg4 (tmp) | OK | 7/7 PASS | Duplicate |
| feat/ej230-sslg4 | feat_ej230-sslg4 (tmp) | OK | ALL PASS | |
| exp/pair-scan-2026-06-11 | exp_pair-scan (tmp) | OK | ALL PASS | |
| feat/ej204-event-display-tracks | feat_ej204 (tmp) | OK | ALL PASS | |
| feat/endonly-mylar | feat_endonly-mylar (tmp) | FAIL | FAIL | Conflict markers prevent compilation |
| diag/photon-budget | diag_photon-budget (tmp) | OK | 1 FAIL | smoke_test passes but reflector volumes STILL PRESENT; physics_baseline FAILS |
| backup/photon-budget-worktree | backup (tmp) | FAIL | FAIL | Conflict markers |
| fix/physics-baseline | fix_physics-baseline (tmp) | FAIL | FAIL | Conflict markers |
| fix/readout-wrapping-config | fix_readout-wrapping-config (tmp) | FAIL | FAIL | Conflict markers |
| main | main (tmp) | FAIL | FAIL | Conflict markers + API mismatch |

**Warning on diag/photon-budget:** The build succeeds but the simulation still runs reflector volumes — the `physics_baseline_check` failure confirms the fix was not effective. A passing smoke_test alone does NOT verify the optics fix.

---

## 6. Constants and Physical Parameters Check

The following constants were checked across all modified `DetectorConstruction.cc` files. None were changed by the fix commits:

| Parameter | Status |
|-----------|--------|
| RINDEX (all materials) | UNCHANGED |
| GROUPVEL | UNCHANGED |
| ABSLENGTH | UNCHANGED |
| Scintillation yield | UNCHANGED |
| Scintillation rise/decay times | UNCHANGED |
| PDE | UNCHANGED |
| Reflectivity (surface property) | UNCHANGED |
| sigma_alpha | UNCHANGED |
| Surface finish | UNCHANGED |
| Channel mapping | UNCHANGED |
| Bar dimensions (kBarHalf*) | UNCHANGED |
| SiPM placement | UNCHANGED |

The fix only removes thin foil-volume placements and their border surfaces, replacing them with a skin surface that references the same `reflector` optical surface object.

---

## 7. Conflict Marker Inventory

Files confirmed to contain committed conflict markers (git show HEAD):

| Branch | SHA | File | Lines |
|--------|-----|------|-------|
| feat/endonly-mylar | 5f778ca | include/DetectorConstruction.hh | 22, 29, 35 |
| backup/photon-budget-worktree | 2db2f4f | include/DetectorConstruction.hh | 21, 28, 34 |
| backup/photon-budget-worktree | 2db2f4f | tests/readout_config_check.cc | multiple |
| fix/physics-baseline | d74816e | include/DetectorConstruction.hh | 21, 28, 34 |
| fix/readout-wrapping-config | 838e91f | include/DetectorConstruction.hh | 21, 28, 34 |
| fix/readout-wrapping-config | 838e91f | tests/readout_config_check.cc | multiple |
| main | b681101 | include/DetectorConstruction.hh | 20, 27, 33 |

---

## 8. Recommended Repair Order

1. **Immediate (working tree cleanup only):**
   - `/mnt/d/SHiP/ej200`: `git restore --staged` + `git restore` 4 files → working tree matches correct HEAD
   - `/mnt/d/SHiP/ej200_event_display`: inspect RunAction.cc, then same restore procedure

2. **Short term (revert bad commits, apply targeted fix):**
   - `fix/readout-wrapping-config`: discard staged + `git reset --soft HEAD~1` → lands on correct c9e61af
   - `feat/endonly-mylar`: revert 5f778ca → 7347db6; apply branch-specific minimal fix
   - `diag/photon-budget`: revert 180a9ec → 2dbe63a; apply branch-specific minimal fix
   - `fix/physics-baseline`: revert d74816e → ea0c6d2; apply targeted fix

3. **Deferred (strategy decision needed):**
   - `main`: revert b681101 → 4e46959; defer fix until API strategy agreed

4. **Skip:**
   - `backup/photon-budget-worktree`: revert at convenience

---

## 9. Artifacts Location

```
/home/reriosto/SHiP/orchestrator/recovery_diffs/   — 29 diff files
/home/reriosto/SHiP/orchestrator/recovery_logs/    — 16 log files
```

Key diff files:
- `perm_ej200_feat_endtop-sslg4.diff` — staged revert on main working dir
- `perm_ej200_endonly_feat_endonly-mylar.diff` — staged content for broken endonly
- `perm_ej200_edge_scan_fix_readout-wrapping-config.diff` — large staged revert (265 deletions)
- `perm_ej200_event_display_feat_ej204-event-display-tracks.diff` — staged + unstaged
- `head_feat_endonly-mylar.diff` — content of bad commit on endonly-mylar
- `head_main.diff` — content of bad commit on main
- `head_fix_physics-baseline.diff` — content of bad commit on fix/physics-baseline

---

## 10. Final Status Table (Post-Recovery)

| Branch | Pre-Copilot SHA | Final SHA | Verdict | Build | Tests | Pushed |
|--------|-----------------|-----------|---------|-------|-------|--------|
| feat/endtop-sslg4 | c6e7843 | f39b84c | **SAFE TO KEEP** | OK (tmp) | 7/7 PASS | NO |
| feat/ej230-sslg4 | 79e701d | 264d263 | **SAFE TO KEEP** | OK (tmp) | ALL PASS | NO |
| exp/pair-scan-2026-06-11 | 021d5f4 | c7a627e | **SAFE TO KEEP** | OK (tmp) | ALL PASS | NO |
| feat/ej204-event-display-tracks | c93502b | b2718c4 | **SAFE TO KEEP** | OK (tmp) | ALL PASS | NO |
| feat/endonly-mylar | 7347db6 | **6942ee6** | **REPAIRED** | OK | 11/11 PASS | NO |
| fix/readout-wrapping-config | 2a9d57b | **c9e61af** | **REPAIRED** | OK | 4/4 PASS | NO |
| diag/photon-budget | 2dbe63a | 2dbe63a (reverted) | **REVERTED** (fix pending) | N/A | N/A | NO |
| fix/physics-baseline | ea0c6d2 | ea0c6d2 (reverted) | **REVERTED** (fix pending) | N/A | N/A | NO |
| main | 4e46959 | 4e46959 (reverted) | **REVERTED** (deferred) | N/A | N/A | NO |
| backup/photon-budget-worktree | 81e242d | 81e242d (reverted) | **SKIPPED** (archival) | N/A | N/A | NO |

All branches are unpushed. No remote damage has occurred.

## 11. Recovery Actions Taken (2026-06-18)

### Working tree restorations (no new commits needed):
- `feat/endtop-sslg4` `/mnt/d/SHiP/ej200`: `git restore --staged` + `git restore` — 4 files reverted to correct HEAD
- `feat/ej204-event-display-tracks` `/mnt/d/SHiP/ej200_event_display`: same — 5 files restored to HEAD
- `fix/readout-wrapping-config` `/mnt/d/SHiP/ej200_edge_scan`: `git reset c9e61af` + `git restore` — removed bad cherry-pick commit 838e91f; working tree restored to c9e61af

### Branch reverts (removed bad Copilot commits, returned to pre-fix state):
- `feat/endonly-mylar`: `git reset 7347db6` — removed bad cherry-pick 5f778ca (had conflicts in include)
- `diag/photon-budget`: `git reset 2dbe63a` — removed bad cherry-pick 180a9ec (clobbered API)
- `fix/physics-baseline`: `git reset ea0c6d2` — removed bad cherry-pick d74816e (had conflicts in include)
- `main`: `git reset 4e46959` — removed bad cherry-pick b681101 (had conflicts + API mismatch)
- `backup/photon-budget-worktree`: `git reset 81e242d` — removed bad cherry-pick 2db2f4f (archival branch)

### New repair commits:
- `feat/endonly-mylar` → `6942ee6`: correct SkinSurface fix preserving SetMylarReflectivity/SetTopSurface/SetMylarSigmaAlpha APIs; 11/11 ctest tests pass

### Pending (require separate work sessions):
- `diag/photon-budget`: apply minimal fix for this branch's TopSiPMCenterX(idx, pitch, nTotal) API
- `fix/physics-baseline`: apply targeted fix (similar API to fix/readout-wrapping-config)
- `main`: decide strategy for dynamic-pitch API before applying fix
