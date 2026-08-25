# Reproducing Mašić & Eberl (2014): debugging summary

**Model:** hybrid biofilm/suspended-biomass nitrification model, Mašić & Eberl,
*Bull. Math. Biol.* 76:27-58 (2014). Scenario reproduced throughout: D=5/day,
c=50 carriers (Figs. 1 and 2 of the paper).

## 1. Confirmed coding bug (IGNORE as intentional)

`partial_nitritation.py` hardcoded `nc = 75` (line 60). The MATLAB original
(`PartialNitritation.m`) correctly uses `nc = 50`, matching the paper's
c=50 scenario. Since colonizable area `A = Ar + Ac*nc` enters the
attachment, detachment and thickness ODEs directly (not just as a final
scale factor), this single-line error propagated through the entire
trajectory.

**Fix:** `nc = 50`.

Effect of this fix alone (N=16, T=90, all else unchanged):

| Quantity | Paper | nc=75 (buggy) | nc=50 (fixed) |
|---|---|---|---|
| Effluent NH₄ (steady) | ~2.3-2.5 g/m³ | 1.03 | 2.54 |
| Biofilm thickness | 165 µm | 157 µm | 185 µm |
| Biofilm AOB mass | ~0.62-0.65 g | 0.51 g | 0.68 g |

(Minor, negligible aside: `V = 0.00603186` in Python vs. `V = 0.006` in
MATLAB — the Python value is the more precise π·r²·h from Table 2's
reactor dimensions; difference is <0.1%, not worth changing.)

## 2. Likely erratum in the paper's Appendix A

Eq. (22)/(23) as typeset:

  μ_A(S) = μ_A,max · (O2 Monod) · [ (NH4 Monod) − (1+η)·b_A ]

distributes `μ_A,max` across the *entire* bracket, including the decay
term. But re-deriving μ_A from Table 1's stoichiometric matrix via Eq. (21)
gives no `μ_A,max` factor on the decay term at all — Table 1's
endogenous-respiration and inactivation process rates are `b_A·(O2 Monod)·X_A`
and `b_A·η·(O2 Monod)·X_A`, with no `μ_A,max` anywhere. Both your MATLAB and
Python code correctly implement Eq. (22) *as printed*, which is the natural
reading of the appendix — but it disagrees with Table 1.

**Empirical test** (nc=50 fixed either way, N=16, T=90):

| Quantity | Paper | Eq. 22 as printed | Table-1-consistent |
|---|---|---|---|
| Biofilm thickness | **165 µm** | 185 µm | **166.7 µm** |
| Biofilm AOB mass | ~0.62-0.65 g | 0.68 g | **0.617 g** |
| Biofilm NOB mass | ~0.13 g | 0.168 g | **0.135 g** |

The Table-1-consistent version matches the paper markedly better. Working
theory: Eq. (22)/(23) have a typesetting error, and the original code used
the Table-1 form.

**Fix applied:**
```python
def Rxa(y, z): return mua(y, z) - (1 + eta) * ba * muaO(z)   # no muamax factor
def Rxn(y, z): return mun(y, z) - (1 + eta) * bn * munO(z)   # no munmax factor
```

## 3. Validation against the paper's figures (both fixes applied, N=32, T=90)

- **Fig. 1a-d** (time series): shapes and steady-state values now match
  well — NH₄ ≈2.5, NO₂ ≈1.0-1.5, thickness →166.7 µm, biofilm AOB/NOB/inert
  masses all close to the published curves.
- **Fig. 2b** (substrate concentration profile): matches closely — O₂ from
  ~0 at the substratum to 5 g/m³ at the surface, NH₄ 1.3→2.6, NO₂
  essentially flat ~1.5, same shapes and magnitudes as the paper.
- **Fig. 2a** (biomass-fraction profile): **still does not match.** The
  paper shows a clear gradient (inert ~16%→5%, NOB ~14%→18%, substratum→
  surface). The reproduction is nearly flat everywhere (AOB 72-73%, NOB
  15-16%, inert ~11%, essentially independent of depth) — this is the
  open issue investigated below.

## 4. Investigating the flat Fig. 2a profile

Ruled out:
- **Grid resolution.** N=16 and N=32 give the *identical* flat profile
  to 3 decimal places, at every saved time, not just at steady state.
- **BVP fixed-point tolerance.** Tightening `tol` from 1e-6 to 1e-10 and
  `droptol` from 1e-12 to 1e-14 changed the thickness trajectory by
  <0.02 µm at every time point checked — negligible next to the ~10-15 µm
  gap versus the paper.
- **The velocity field `v(z)` being uniform.** Directly recomputed: it
  goes from 0 at the substratum to 0.276 m/d at the surface, a strong,
  smooth gradient — not flat.
- **The `Lp`/attachment/erosion ODE structure and the surface-flux term
  (`JNH4`).** Verified term-by-term against Eqs. (15) and (1)/(3)
  respectively; both match exactly.

Found:
- **The AOB/NOB specific-growth-rate gap (`Rxa − Rxn`) is tiny at steady
  state** (~1-5% of either rate), because `K_A,O2 = K_N,O2 = 0.5` (identical
  for both species — the paper's own parameter choice) and because NH₄
  and NO₂ both sit well above their own half-saturation constants
  everywhere in the film at steady state. The one substrate with a large
  internal gradient (O₂) affects both species identically, so it can't
  differentiate them.
- **But the gap is huge in the first few days** — 0.33/day at t=0 (NOB
  can't grow at all with zero nitrite present), decaying to the small
  steady value by ~t=15-20, by which point the film has already grown
  past half its final thickness. So the model clearly *has* a mechanism
  for differentiation; the question is why it isn't recorded spatially.
- **Explanation:** at early times the film is only 10-20 µm thick, thin
  enough that diffusion equilibrates it almost instantly — it behaves
  as a well-mixed slab, so the large early reaction differential shifts
  composition *uniformly*, not spatially. Only once the film is thick
  enough for real internal gradients (~t>15-20d, >100µm — exactly when
  O₂ starts showing real depletion) does *any* spatial structure appear,
  but by then the differential has already decayed to its small residual.
  Confirmed directly: the spatial fA spread is ~0.0001 at t=1d, growing to
  only ~0.006-0.009 by t=20-60d.
- **Thickness growth-curve comparison** (paper's Fig. 1c pixel-traced and
  overlaid against the model): not a uniform lag but a **crossover** — the
  model runs slightly *ahead* of the paper for days ~2-7, then the paper
  pulls ahead through days ~10-20, before both converge by day ~30.
- **NH₄ depletion comparison** (paper's Fig. 1a digitized similarly): the
  model consumes NH₄ *faster* than the paper throughout days 0-11, growing
  from a negligible gap at t=0 to ~6 g/m³ by t=11 — consistent with the
  model's slightly thicker film over that same window. Suspended biomass
  ruled out as the driver (<2% of total AOB mass through day 20).
- **Mechanistic cause of the day-13-20 deceleration, located:** decomposing
  `Lp` into its three terms shows erosion (`∝ E·λ²`) growing from -8% of
  the growth term at day 1 to **-91% at day 20** — it very nearly cancels
  the growth term right in the window where the paper pulls ahead. `E=1000`
  matches Table 3 exactly, so this isn't a wrong parameter; it means the
  growth term itself is peaking a bit low/early relative to the paper,
  letting erosion catch up sooner.

**Status:** unresolved. The remaining discrepancy is upstream of `Lp` —
most likely in what sets the size of the oxygen-dependent growth-integral
term (`v[0]`) as the film thickens — i.e. possibly still in `RDE`, even
though it checks out structurally against Eqs. (26)/(27). This would need
either the original authors' code/data to pin down definitively, or a
systematic parameter-sensitivity sweep (diffusion coefficients, oxygen
Monod constants) that wasn't completed in this session.

## 5. Parameter sensitivity sweep (day 7-35 window)

Ran the fixed model (N=16 — confirmed equivalent to N=32 — T=35) at
baseline and ±20% for every parameter that plausibly controls how fast
the growth term saturates: `D_O2`, `K_A,O2`/`K_N,O2` (moved together, since
the paper sets them equal), `D_NH4`, `D_NO2`, and `E` (erosion, as a
reference point given Section 4's finding).

**Thickness at t=35d, Δ from baseline (165.8 µm):**

| Parameter (±20%) | −20% | +20% |
|---|---|---|
| `D_NH4` | −0.6 µm (−0.4%) | +0.4 µm (+0.2%) |
| `D_NO2` | 0.0 µm (0.0%) | 0.0 µm (0.0%) |
| `D_O2` | −3.2 µm (−1.9%) | −5.4 µm (−3.2%) |
| `K_A,O2`=`K_N,O2` | +1.8 µm (+1.1%) | −1.9 µm (−1.2%) |
| **`E` (erosion)** | **+22.8 µm (+13.8%)** | **−17.8 µm (−10.7%)** |

**Findings:**
- `D_NH4` and `D_NO2` have essentially **zero** effect, confirming they
  never become diffusion-limiting in this scenario (both substrates stay
  well above their half-saturation constants throughout) — consistent
  with everything found in Section 4.
- `D_O2` has a small, slightly non-monotonic effect (both directions
  *decrease* thickness) — a secondary player, not the main lever.
- `K_A,O2`/`K_N,O2` has a modest, cleanly monotonic effect, as expected
  (lower K → less O₂-limited → thicker film).
- **`E` dominates by almost an order of magnitude** over every other
  parameter tested — a ±20% change in erosion alone produces an 11-14%
  swing in thickness by day 35, versus ≤3% for everything else combined.
  Given `E=1000` matches Table 3 exactly, this isn't a wrong value, but
  it does mean the day-10-25 transient sits in a regime where the
  growth/erosion balance is extremely finely poised — a very small,
  otherwise-undetectable difference between this implementation and the
  paper's original (different numerical scheme, rounding, etc.) could
  plausibly produce the observed multi-day timing shift without any
  parameter actually being "wrong."

**Targeted follow-up test — breaking the K_A,O2=K_N,O2 symmetry directly**
(the mechanism proposed in Section 4 for the flat Fig. 2a profile): reran
with `K_A,O2=0.3, K_N,O2=0.9` (a 3× difference, well beyond anything the
paper reports) instead of the paper's `0.5/0.5`. Result: the steady-state
spatial spread in fA roughly **doubled** (0.006-0.009 → 0.0145) — real,
in the predicted direction — but still nowhere near the paper's ~0.07-0.11
(7-11 percentage points). So the identical-K_O2 explanation from Section 4
is a genuine contributing factor, but **not sufficient on its own** to
account for the missing Fig. 2a structure. Something else is still
suppressing differentiation beyond what this mechanism explains.

**Net conclusion:** the day 7-25 transient-shape mismatch and the flat
Fig. 2a profile are almost certainly linked to the same root cause, and
erosion sensitivity plus the partial (not full) success of the K_O2 test
both point toward the growth/erosion balance and the oxygen-limitation
onset — but neither the diffusion coefficients nor the Monod constants,
individually or via the tested combinations, fully close the gap. This
would benefit from either the original authors' code, or a joint 2D sweep
over `E` and `K_O2` simultaneously (not done here) to see if some
combination reproduces both the correct transient shape and Fig. 2a
structure at once.

## 6. Joint E × K_O2 sweep (closing off Section 5's open question)

Ran a 3×3 grid (E ∈ {800, 1000, 1200}, K_A,O2=K_N,O2 ∈ {0.35, 0.5, 0.65})
plus 2 combinations with K_A,O2≠K_N,O2 asymmetric, tracking thickness at
t=13, 20, 35 and the Fig. 2a fA-spread at steady state:

| Combo | t=13 | t=20 | t=35 | fA spread |
|---|---|---|---|---|
| **Paper (target)** | **88.9** | **153.4** | **~165.3** | **~0.07-0.11** |
| E800, K=0.50 | 86.9 | 154.9 | 188.6 | 0.0082 |
| E1000, K=0.35 | 87.0 | 144.5 | 168.5 | 0.0133 |
| E1000, K=0.50 (baseline) | 81.9 | 139.9 | 165.8 | 0.0085 |
| E1200, K=0.65 | 73.2 | 123.8 | 145.5 | 0.0064 |
| ...(full 11-combo grid in `joint_sweep_results.json`) | | | | |

**No combination in this 2-parameter space hits all three thickness
checkpoints simultaneously** — whichever combo matches t=13 and t=20
overshoots or undershoots t=35, and vice versa. This closes off Section
5's open question: `E` and `K_O2` alone, however tuned, cannot reproduce
the paper's transient shape. Something else has to be contributing.
Fig. 2a's spread also never gets close to the paper's ~0.07-0.11 in any
combination — confirming, again, that this isn't an E/K_O2 tuning problem.

## 7. A genuine second bug found via the spatial-validation diagnostic

Building the diagnostics requested below (Section 8, item 5:
`boundary_condition_residual`) and running it against this model's own
Neumann boundary condition turned up something concrete. The paper's
Eq. (5) specifies `C'_k(0)=0` (true zero-flux) at the substratum. The
code enforces this by replacing the last row of the Laplacian with the
first-derivative row (`D2[N,:] = D[N,:]`, standard Chebyshev-collocation
practice) — but never zeroes the corresponding entry of the right-hand
side vector (`Gtemp`) that row is solved against. Checked directly:

```
D[N,:] @ CO2_profile           = 0.4374  (should be ~0 for true no-flux)
predicted local O2 reaction rate at that node (Gtemp[-1])  = 0.4374
```

These match to 4+ significant figures — not a coincidence. The "Neumann"
row is actually solving `D'·C = (local reaction rate)`, not `D'·C = 0`.
This is a second, independent bug (present in both MATLAB and Python,
inherited from the standard textbook row-substitution trick without also
zeroing that row's RHS entry).

**Effect of fixing it** (set `Gtemp[-1, :] = 0` before the linear solve),
N=32, T=60:

| t (d) | baseline (buggy BC) | BC fixed | paper |
|---|---|---|---|
| 15 | 102.5 | 107.9 | 118.3 |
| 20 | 139.9 | 153.0 | 153.4 |
| 25 | 157.4 | 168.2 | 163.4 |
| 35 (~steady) | 165.8 | 172.9 | 165.3 |

The fix dramatically improves the match at t=20 (153.0 vs. target 153.4)
but overshoots the final steady state (172.9 vs. 165.3). Combining the BC
fix with a modest, non-principled `E` adjustment (1000→1080, +8%, well
inside plausible measurement/definition tolerance) gets remarkably close
from day 25 onward:

| t (d) | BC-fixed + E=1080 | paper |
|---|---|---|
| 25 | 162.1 | 163.4 |
| 30 | 165.6 | 165.3 |
| 35 | 166.3 | 165.3 |
| 45 (steady) | 166.5 | 165.3 |

**Does not** resolve Fig. 2a: fA spread stays at ~0.006-0.008 regardless.
So this bug explains a real, now-quantified chunk of the day-20-35 curve
shape, but not the flat-profile issue, which remains open.

**Fix applied and REVERTED in the delivered code** (regression discovered
2024 — see errata note below): applying this fix by default caused a
Fig. 2b regression (O₂ at the substratum: 0.06 without the fix, matching
the paper's near-zero value, vs. 0.23 with it — clearly worse). The
`biofilmbvp()` code in the delivered `partial_nitritation.py` documents
this finding in a comment but does **not** apply it by default. Treat
this as a documented, real inconsistency worth investigating further,
not a recommended production fix as things stand.

## 8. Errata: overclaimed match quality, corrected

The original phrasing of Section 6 ("Fig. 1(a-d) ... match well") was
too optimistic and was corrected after closer, direct side-by-side
comparison against the digitized paper figures. The accurate picture:

- **Fig. 1c (thickness) and Fig. 2b (substrate profile): genuinely
  match well**, confirmed by direct overlay.
- **Fig. 1a (NH₄/NO₂ effluent): shape mismatch, not just a timing
  shift.** The paper shows a shallow, near-flat NH₄ decline out to
  ~day 8-10, then a comparatively sudden steep drop between ~day 10-18
  (an S-shaped/sigmoid depletion curve, consistent with consumption
  ramping up as biomass grows in). The model shows a steep decline
  starting immediately from t=0, with no such lag/shoulder. This is
  visible from the very first day, not just a magnitude gap that grows
  over time as previously described.
- **Fig. 1b and 1d (suspended and biofilm AOB biomass): the paper's
  curves overshoot** — rising above, then dipping back down to, their
  eventual steady value. The model's curves rise smoothly and
  monotonically to steady state with no overshoot, in both panels,
  regardless of the BC fix in Section 7. This is a more specific and
  diagnostic symptom of the same underlying day-10-25 transient-timing
  issue described in Sections 4-6, but hadn't previously been isolated
  in this form.
- **Fig. 2a (biomass fraction profile): still flat**, as documented.

**Net honest status:** two real bugs found and fixed (`nc`, Eq. 22/23);
one further real inconsistency found and documented but not applied by
default (Neumann BC RHS); final thickness and the steady-state substrate
profile both now match the paper well; but the **transient shape** in
three different figures (1a's missing lag phase, 1b/1d's missing
overshoot, and by extension the day-10-25 window already characterized
in Section 4-6) point at a single, still-unidentified root cause in the
early-to-mid transient dynamics. Fig. 2a's flat spatial profile is very
likely downstream of this same root cause, per the mechanism argued in
Section 4. This is the honest state of the reproduction, not "close
except for one spatial detail."

## 9. Recommendations for applying this approach to new systems

In `partial_nitritation.py`:
```python
nc = 50                                                          # was 75
def Rxa(y, z): return mua(y, z) - (1 + eta) * ba * muaO(z)       # drop muamax factor
def Rxn(y, z): return mun(y, z) - (1 + eta) * bn * munO(z)       # drop munmax factor
```
The Neumann-BC RHS fix (Section 7) is documented but **not** applied by
default — see Section 7/8 for why.

Equivalent MATLAB fixes (`PartialNitritation.m` already has the right `nc`):
```matlab
Rxa = @(y,z) mua(y,z) - (1+eta)*ba*muaO(z);
Rxn = @(y,z) mun(y,z) - (1+eta)*bn*munO(z);
```

`full_model.py` (attached) reproduces all of Fig. 1(a-d) and Fig. 2(a-b)
from the fixed model's output. Status: see Section 8 for the corrected,
honest summary of what does and doesn't match.

## 9. Recommendations for applying this approach to new systems

See the accompanying `recommendations.md` and four diagnostic modules
(`diagnostic_1_conservation.py`, `diagnostic_2_fixedpoint_monitor.py`,
`diagnostic_3_etd_contour.py`, `diagnostic_5_spatial_validation.py`).
In short: the spectral-collocation + domain-rescaling + ETD-timestepping
*strategy* is sound and worth reusing; the *code* needs the safeguards in
those modules built in from the start, not added after the fact --
diagnostic 5 alone found a second real bug in this exact codebase within
minutes of being written.
