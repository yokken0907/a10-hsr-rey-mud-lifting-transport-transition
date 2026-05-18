# Reviewer Start Here

## Purpose of this repository

This repository is a paper companion and review package for **Hybrid Smart Riser / HSR v3.x**.

HSR v3.x is a reduced-order, pre-engineering transport-transition and safety-diagnostic framework for ultra-deep REY mud lifting. It is not a sea-trial-validated mining system, not a final riser design, not certified offshore equipment, and not a commercial-deployment package.

The intended review posture is therefore:

> Review this repository as a pre-engineering surrogate and diagnostic framework for transport-transition risk in ultra-deep REY mud lifting, not as a completed offshore mining technology.

## One-paragraph summary

HSR v3.x studies ultra-deep rare-earth-element and yttrium (REY) mud lifting as a transport-transition problem. Within the current reduced-order model, it separates sustained slurry transport, blackout/coastdown-induced collapse, and deposition-limited long-duration constraints. Its practical value is to provide an early diagnostic and screening framework before high-cost field trials, detailed 3D CFD/FSI, environmental assessment, safety certification, or hardware design.

## What the repository is allowed to support

The repository supports the following limited claims:

- HSR v3.x is a reduced-order pre-engineering transport-transition framework.
- It examines sustained transport, fault-induced collapse, blackout/coastdown response, and deposition-limited behavior in a tested surrogate.
- In the included summary tables, the no-blackout critical inlet-velocity region is approximately `uss ≈ 1.25–1.26`.
- In the included no-blackout comparison, `uss = 1.6` remains in sustained transport for 40 s, 80 s, and 120 s within the tested reduced-order model.
- Deposition growth remains a next-stage constraint rather than a solved long-duration management problem.
- The residual pressure-gradient time derivative is treated only as a **candidate diagnostic trigger** in the tested model.

## What the repository does not claim

This repository does **not** claim:

- sea-trial validation,
- field-validated deep-sea mining,
- final riser hardware design,
- certified offshore equipment,
- production-ready or commercial-ready mining technology,
- environmental permitting readiness,
- completed long-duration deposition management,
- replacement of 3D CFD/FSI, full offshore engineering review, laboratory slurry-loop tests, or field trials.

## Recommended reading order

For a first technical review, read in this order:

1. `CLAIM_BOUNDARY.md`  
   Establishes the allowed and forbidden claim language.

2. `FIELD_VALUE_en.md` or `FIELD_VALUE_ja.md`  
   Explains what problem the theory could help with in practical evaluation.

3. `INDUSTRY_RELEVANCE_ja.md`  
   Gives a Japanese industry-facing interpretation for public-agency or company discussions.

4. `NEXT_VALIDATION_PLAN.md`  
   Shows what must be tested next before any stronger engineering claim is made.

5. `results_summary/summary.json` and the CSV files in `results_summary/`  
   Provide compact result tables used for repository-level orientation.

6. `manuscript/`  
   Contains the manuscript PDF and related materials.

7. `source_materials/`  
   Preserves selected scripts, logs, JSON files, and closure/source evidence used to support the paper companion archive.

8. `docs/technical_visual_orientation/index.html`  
   Browser-only visual orientation page for the transport-transition logic. This is an orientation aid, not an executable engineering simulator.

## Evidence hierarchy

When evaluating claims, use the following hierarchy:

| Level | Material | Role |
|---|---|---|
| 1 | `CLAIM_BOUNDARY.md` / `LIMITATIONS.md` | Defines the maximum allowed interpretation |
| 2 | Manuscript PDF | Main narrative and technical synthesis |
| 3 | `results_summary/*.csv` / `summary.json` | Compact numerical summaries |
| 4 | `source_materials/` | Selected scripts/logs/source evidence |
| 5 | `docs/technical_visual_orientation/index.html` | Visual orientation only |
| 6 | Outreach and field-value files | Reader-specific explanation, not stronger evidence |

If an outreach document sounds stronger than `CLAIM_BOUNDARY.md`, the claim-boundary file takes priority.

## Quick technical orientation

The current reduced-order result package should be read as follows:

- The no-blackout sweep indicates a transition around `uss ≈ 1.25–1.26`.
- The `uss = 1.6` no-blackout cases remain positive-flow sustained cases over 40 s, 80 s, and 120 s in the reduced-order model.
- Longer duration increases deposition thickness; therefore, deposition creep is not solved and remains a validation target.
- Blackout/coastdown cases show a much larger residual pressure-gradient derivative than the no-blackout reference in the included diagnostic table.
- The diagnostic trigger is promising only within the tested assumptions; sensor drift, missing data, calibration error, time delay, and harsh offshore measurement conditions require future testing.

## Reviewer checklist

Before accepting any practical interpretation, check:

- Is the statement explicitly restricted to the reduced-order model?
- Does it avoid field-validation language?
- Does it avoid final-design or certification language?
- Does it distinguish diagnostic screening from operational control?
- Does it acknowledge deposition management as incomplete?
- Does it point to the next validation stage rather than implying deployment readiness?

## Safe citation-style summary

> HSR v3.x is a reduced-order pre-engineering framework for analyzing transport-transition risk in ultra-deep REY mud lifting. The included results distinguish sustained transport, blackout/coastdown collapse, and deposition-limited constraints within a tested surrogate. The repository does not claim sea-trial validation, final riser design, safety certification, environmental permitting readiness, or commercial deployment readiness.
