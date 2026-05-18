# Claim-to-Evidence Map

## Purpose

This file maps each safe claim to the evidence category currently available in the repository.

It is intended to prevent overclaiming. If a claim is not listed here, it should not be made without additional evidence.

## Claim boundary priority

If there is any conflict between this file and `CLAIM_BOUNDARY.md`, `CLAIM_BOUNDARY.md` takes priority.

## Safe claim map

| Safe claim | Current evidence location | Maximum interpretation |
|---|---|---|
| HSR v3.x is a reduced-order pre-engineering transport-transition framework | README, manuscript, claim-boundary documents | Framework-level claim only |
| It studies ultra-deep REY mud lifting as a transport-transition problem | Manuscript, README, field-value files | Problem-framing claim |
| It separates sustained transport, blackout/coastdown collapse, and deposition-limited behavior | Manuscript, `results_summary/summary.json`, summary CSV files | Surrogate-level classification |
| The no-blackout critical inlet-velocity region is approximately `uss ≈ 1.25–1.26` | `results_summary/table1_no_blackout_40s_sweep.csv` | Current reduced-order model only |
| `uss = 1.6` remains sustained over 40 s, 80 s, and 120 s in no-blackout cases | `results_summary/table2_long_no_blackout_comparison.csv` | Current reduced-order model only |
| Deposition thickness grows with longer duration | `results_summary/table2_long_no_blackout_comparison.csv` | Deposition remains next-stage constraint |
| Conservative Ucrit sensitivity does not overturn the `uss = 1.6` sustained result in included cases | `results_summary/table3_ucrit_conservative_sensitivity.csv` | Included sensitivity only |
| Residual pressure-gradient derivative is a candidate blackout diagnostic trigger | `results_summary/table4_diagnostic_noise_delay_robustness.csv` | Candidate trigger, not operational alarm |
| The repository is suitable for specialist review and next-stage validation planning | README, claim boundary, next validation plan | Review package only |

## Forbidden claim map

| Forbidden claim | Why it is not supported |
|---|---|
| Field-validated deep-sea mining system | No sea-trial evidence is included |
| Final riser hardware design | No certified engineering design package is included |
| Certified offshore equipment | No safety certification evidence is included |
| Production-ready mining equipment | No production design, manufacturing, deployment, or operations evidence is included |
| Commercial deployment readiness | No field validation, safety certification, environmental permitting, or operational approval is included |
| Environmental permitting readiness | No environmental impact assessment or permitting process is included |
| Complete long-duration deposition management | Deposition growth remains a next-stage constraint |
| Replacement for full CFD/FSI | The model is reduced-order and requires higher-fidelity comparison |
| Approved blackout alarm logic | Diagnostic robustness is not yet validated under realistic field instrumentation |
| Safe stop/restart procedure | Restart and operational safety procedures are not validated |

## Recommended wording

### English

> In the current reduced-order model, HSR v3.x distinguishes sustained transport, blackout/coastdown-induced collapse, and deposition-limited constraints. These results support pre-engineering screening and validation planning, not field deployment or certified equipment claims.

### Japanese

> 現在のreduced-order modelにおいて、HSR v3.xは持続輸送、blackout/coastdownに伴う崩壊、および堆積制約を区別する。これらの結果は前工学的な候補選別と検証計画を支えるものであり、実海域導入や認証済み装置を主張するものではない。

## Review rule

A practical statement is acceptable only if it passes all three checks:

1. It states or implies the reduced-order nature of the evidence.
2. It avoids field-validation, certification, final-design, and commercial-readiness language.
3. It identifies what additional validation is required.
