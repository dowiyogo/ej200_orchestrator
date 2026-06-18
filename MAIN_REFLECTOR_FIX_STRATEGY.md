# main Branch — Reflector-Volume-to-SkinSurface Fix Strategy

**Date:** 2026-06-18  
**Starting SHA:** 4e46959  
**Branch:** main

---

## 1. How main computes the number of top SiPMs

`main` uses a **dynamic pitch** mechanism:

```cpp
// In DetectorConstruction.hh:
static G4double TopSiPMCenterX(G4int idx, G4double pitch, G4int nTotal);
static G4int    ComputeNTopSiPMs(G4double pitch);  // private
void     SetTopSiPMPitch(G4double pitchMm);
G4int    GetNTopSiPMs() const { return fNTopSiPMs; }
```

`fNTopSiPMs` is recomputed at each `Construct()` call via `ComputeNTopSiPMs(fTopSiPMPitch)`. Default pitch = 70 mm → 20 Top SiPMs.

This is the same API as `diag/photon-budget` and `fix/physics-baseline`. It is **not** the same as `feat/endtop-sslg4`'s fixed-N API (`kNTopSiPMs = 70`). The dynamic pitch must be preserved.

---

## 2. How reflector panels are currently constructed (4e46959)

Located in `src/DetectorConstruction.cc` lines ~152–196:

```cpp
auto* reflector = Materials::CreateBarSkinReflector();
const G4double foilHalfT = 0.5 * um;
// 5 solids: reflYMinus, reflZ (shared), reflX (shared)
// 5 placements: ReflectorYMinusPV, ReflectorXMinusPV, ReflectorXPlusPV,
//               ReflectorZMinusPV, ReflectorZPlusPV
// 5 border surfaces: BarReflector_YMinus, _XMinus, _XPlus, _ZMinus, _ZPlus
```

**No +Y reflector panel** — the +Y face is left open for Top SiPMs.  
**No X panels in instrumented condition** — all 5 panels are unconditional (X panels always placed regardless of SiPM configuration).

---

## 3. Bar-to-SiPM border surfaces that must be preserved

Located after the reflector block (lines ~209–260):

- `fSiPMSurfaces[leftId]` → BarPV→EndSiPMLeft_PV border surfaces (8 per side)
- `fSiPMSurfaces[rightId]` → BarPV→EndSiPMRight_PV border surfaces
- `fSiPMSurfaces[globalId]` → BarPV→TopSiPMPV border surfaces (fNTopSiPMs entries)

These are in `barLV` daughters and are addressed via `fBarPhys`. The fix touches only the reflector block; all SiPM surfaces are downstream and unaffected.

---

## 4. Where to attach the G4LogicalSkinSurface

After `barLV` is created (line ~143) and after `fBarPhys` is placed (line ~149), immediately replacing lines 152–196:

```cpp
// Apply reflector properties as a skin on the bar logical volume.
auto* reflector = Materials::CreateBarSkinReflector();
auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV, reflector);
(void)barSkin;
```

Add `#include "G4LogicalSkinSurface.hh"` to includes.

---

## 5. Tests that need updating

**Tracked tests at 4e46959:**
- `smoke_test` → runs `ej200_bar_sim -m macros/test.mac` — does not check reflector geometry ✓
- `edge_scan_smoke` → runs `ej200_bar_sim -m macros/scan_edge.mac` — does not check reflector geometry ✓

**Untracked Copilot additions:**
- `tests/check_endtop_balance.py` — not in CMakeLists.txt; not compiled by ctest
- `tests/readout_config_check.cc` — not in CMakeLists.txt; not compiled by ctest

**Conclusion:** No tracked test needs updating. The fix is purely a source change.

---

## 6. Decision: Fix now or defer?

**Fix now.** The analysis confirms:
- `main` uses the same API as `diag/photon-budget` and `fix/physics-baseline`
- The fix is minimal and bounded (exactly the same change as those branches)
- The two tracked tests are smoke-tests that do not reference reflector volumes
- No API conflict with the SkinSurface approach

The earlier "defer" recommendation was based on the mistaken assumption that `main` used a different API (`ComputeNTopSiPMs` dynamic) from `feat/endtop-sslg4` (fixed `kNTopSiPMs`). While that difference is real, it is irrelevant to the reflector fix — `ComputeNTopSiPMs` is used for Top SiPM placement, not for reflector volume placement, and is untouched by the fix.

---

## Summary

Apply the same targeted fix as `fix/physics-baseline` and `diag/photon-budget`:

1. Add `#include "G4LogicalSkinSurface.hh"`
2. Replace the 5-panel reflector block with `new G4LogicalSkinSurface("BarSkin", barLV, reflector)`
3. Update edgeWrap message string
4. Update comment block in `include/DetectorConstruction.hh`
5. Build and run `ctest --output-on-failure`
6. Commit with message `fix(optics): replace reflector volumes with bar skin surface`
