# Next Validation Plan for HSR v3.x

## Scope of this plan

HSR v3.x is a reduced-order, pre-engineering transport-transition and safety-diagnostic framework for ultra-deep REY mud lifting. It is not a sea-trial-validated mining system, not a final riser design, not certified offshore equipment, and not a commercial-deployment package.

This file is a staged validation roadmap. It is **not** a deployment plan, not a sea-trial plan, and not an offshore operating procedure.

The purpose is to identify what must be tested before HSR v3.x can support stronger engineering claims.

## Validation posture

The next validation stages should preserve the following principles:

- Do not relax thresholds to rescue a preferred conclusion.
- Do not reinterpret surrogate results as field validation.
- Keep diagnostic claims separate from operational control claims.
- Treat deposition management as unresolved unless tested directly.
- Treat every stronger practical statement as conditional on additional evidence.
- Preserve negative, null, and failure results.

## V0. Repository-level verification

### Goal

Confirm that the public repository is internally consistent and reviewable.

### Tasks

- Check that `README.md`, `README_ja.md`, `CLAIM_BOUNDARY.md`, `LIMITATIONS.md`, and outreach files use consistent claim language.
- Confirm that result summaries match the manuscript tables.
- Confirm that file manifests and source materials correspond to the included archive.
- Ensure that browser-only visual orientation is presented as orientation, not simulation.
- Add a short reviewer path in both README files.

### Exit criterion

Repository can be reviewed without mistaking it for field-validated equipment or final riser design.

## V1. Closure and Ucrit refinement

### Goal

Reduce ambiguity in the critical-velocity and closure interpretation.

### Tasks

- Add an explicit `Ucrit` or equivalent critical-transport column to the closure table if possible.
- Document how the current `Ucrit` estimate was reconstructed.
- Test sensitivity to wall concentration, center concentration, particle radius, and wall shear.
- Compare nominal, conservative, and pessimistic closure variants.
- Identify whether the `uss ≈ 1.25–1.26` transition estimate is robust or model-specific.

### Exit criterion

Critical transport estimates are traceable, reproducible, and not dependent on hidden closure assumptions.

## V2. Diagnostic robustness expansion

### Goal

Test whether the residual pressure-gradient derivative remains useful under more realistic measurement conditions.

### Current status

The included summary table tests additive Gaussian noise, simple delay, and smoothing. This is not sufficient for field instrumentation.

### Tasks

- Add sensor drift.
- Add calibration bias.
- Add missing samples.
- Add quantization error.
- Add irregular sampling.
- Add time-delay sensitivity beyond one sample.
- Add non-Gaussian noise.
- Add false-positive and false-negative audits across no-blackout and blackout cases.
- Compare residual pressure-gradient derivative against simpler baseline alarms.

### Exit criterion

The diagnostic trigger can be described as robust only if it survives realistic sensor degradation without unacceptable false alarms or missed detections.

## V3. Long-duration deposition and restart fragility

### Goal

Move from short-duration sustained transport toward long-duration deposition-risk characterization.

### Tasks

- Extend no-blackout duration beyond 120 s.
- Quantify deposition growth rates.
- Test sensitivity to deposition thresholds.
- Test partial blockage, restart, and re-fluidization scenarios.
- Separate hydraulic arrest from deposition-limited degradation.
- Identify conditions where short-term sustained transport becomes long-term deposition failure.

### Exit criterion

The model can identify deposition-limited horizons without claiming that deposition management is solved.

## V4. Parameter realism

### Goal

Connect the surrogate to plausible REY mud and riser conditions.

### Tasks

- Introduce realistic mud density and solids concentration ranges.
- Introduce particle-size distribution rather than a single representative particle radius.
- Include viscosity and rheology uncertainty.
- Include riser length, pipe diameter, and pressure-envelope scaling.
- Include pump-curve constraints.
- Compare sensitivity across plausible ultra-deep operating regimes.

### Exit criterion

Results can be discussed as plausible parameterized scenarios, while still avoiding field-validation claims.

## V5. Higher-fidelity comparison

### Goal

Compare representative HSR cases with higher-fidelity models.

### Tasks

- Select a minimal benchmark set:
  - no-blackout near-critical case,
  - no-blackout margin case,
  - blackout collapse case,
  - long-duration deposition case.
- Compare with 2D or axisymmetric models if available.
- Compare selected cases with 3D CFD.
- Add FSI checks for riser response if feasible.
- Identify where the reduced model fails.

### Exit criterion

HSR v3.x can be positioned as a screening surrogate only where higher-fidelity comparison does not contradict the qualitative transition classification.

## V6. Laboratory or controlled-loop validation

### Goal

Move from numerical comparison to controlled experimental evidence.

### Tasks

- Define a laboratory slurry-loop analogue.
- Choose measurable channels corresponding to the surrogate states.
- Test sustained, near-critical, and collapse-like regimes.
- Test sensor diagnostic channels.
- Test deposition growth and restart behavior.
- Compare experimental regimes with HSR-classified regimes.

### Exit criterion

The framework gains controlled experimental support, but still does not imply sea-trial validation.

## V7. Field-test design preparation

### Goal

Prepare the measurement logic that would be needed for future field testing.

### Tasks

- Define required sensors and sampling rates.
- Define blackout/coastdown safety scenarios.
- Define stop/restart observation metrics.
- Define criteria for safe termination.
- Define environmental and operational exclusion conditions.
- Identify review requirements by offshore engineering and safety specialists.

### Exit criterion

A future field-test protocol can be discussed with specialists without presenting the current repository as a field-tested system.

## Claim-upgrade ladder

| Evidence stage | Maximum safe claim |
|---|---|
| Current repository | Reduced-order pre-engineering framework |
| V1–V2 complete | Improved surrogate and diagnostic robustness |
| V3–V4 complete | More realistic long-duration and parameterized scenarios |
| V5 complete | Surrogate compared with selected higher-fidelity models |
| V6 complete | Controlled experimental support |
| V7 plus external review | Candidate field-test planning support |
| Actual sea trials and certification | Only then discuss field validation or certified deployment |

## Near-term recommended priority

The most valuable next steps are:

1. Update README reading path and add field-value files.
2. Add explicit claim-to-evidence mapping.
3. Expand diagnostic robustness beyond Gaussian noise/delay/smoothing.
4. Extend deposition duration and restart-fragility analysis.
5. Prepare a small set of benchmark cases for CFD/FSI or laboratory review.

## Final validation statement

HSR v3.x should move forward by increasing evidence quality, not by increasing rhetorical strength.  
The current repository is useful precisely because it can identify what remains unverified.
