# PUBLIC-GATE-0: A10-HSR / REY Mud Lifting Transport-Transition Surrogate

## Gate status

`PASS-WITH-MINOR-PUBLICATION-FIXES-A10-HSR-PUBLIC-GATE-0`

After applying the v0.1.1-public-gate repository patch, the GitHub archive is suitable for conservative public release, subject to replacing the placeholder GitHub account in `CITATION.cff` and leaving Zenodo DOI as pending until release.

## Repository identity

- Repository name: `a10-hsr-rey-mud-lifting-transport-transition`
- Public version: `v0.1.1-public-gate`
- Japanese classification: `超深海REY泥揚鉱・スマートライザー輸送遷移理論`
- Archive type: paper companion archive, not a certified engineering design package

## Positive findings

- `CLAIM_BOUNDARY.md` and `LIMITATIONS.md` are already strong.
- HSR v3.3 is framed as a reduced-order pre-engineering transport-transition theory.
- Sea-trial validation, final riser design, production mining equipment, environmental permitting readiness, and safety certification are explicitly excluded.
- No Bitcoin/donation sentence was detected in repository text or extracted PDF text.
- The manuscript PDF contains a limitations section and does not present the work as a field-validated offshore system.

## Issues found before patching

1. `FILE_MANIFEST.csv/json` used Japanese source paths while the uploaded zip extracted the corresponding directory with escaped `#U...` path segments.
2. `CITATION.cff` still contained `version: v0.1`, `date-released: 2026-XX-XX`, and a placeholder repository URL.
3. Manuscript-facing terms such as `practically actionable`, `完成判定`, and `完成域` are somewhat strong. They are acceptable for GitHub only with explicit non-claim language; Jxiv should receive an additional wording review.
4. Jxiv metadata is still draft-level.

## Repairs applied in v0.1.1-public-gate

- Updated public version to `v0.1.1-public-gate`.
- Added PUBLIC-GATE-0 wording to README / README_ja.
- Added conservative public wording to `CLAIM_BOUNDARY.md`.
- Normalized selected source-material directory to `source_materials/selected_from_uploaded_hsr_project/hsr_simulation_rareearth_selected/`.
- Regenerated `FILE_MANIFEST.csv/json` using manifest-excluding-self policy.
- Added `MANIFEST_NOTE.md` and `tools/verify_manifest_excluding_self.py`.
- Added PUBLIC-GATE documentation under `docs/public_gate_0/`.

## Allowed public claims

A10-HSR / HSR v3.3 may be described as a reduced-order pre-engineering transport-transition framework for ultra-deep REY mud lifting. It may report tested reduced-model behaviors: sustained-transport sector, fault-induced collapse sector, deposition-limited horizon, Ucrit sensitivity checks, and residual pressure-gradient derivative as a candidate diagnostic trigger in the tested reduced model.

## Forbidden public claims

Do not claim sea-trial validation, full 3D offshore validation, final riser hardware design, field-ready operational procedure, certified offshore deployment, environmental permitting readiness, industrial safety certification, production mining equipment readiness, commercial deployment readiness, complete long-duration deposition management, or replacement of full CFD/FSI, field testing, and expert offshore engineering review.

## Standard public wording

> A10-HSR / HSR v3.3 is presented as a reduced-order, pre-engineering transport-transition framework for ultra-deep REY mud lifting. The repository preserves a paper companion archive, compact result summaries, selected source materials, and deferred-material inventory. It does not claim sea-trial validation, final riser hardware design, certified offshore deployment, environmental permitting readiness, production mining equipment readiness, commercial deployment readiness, or complete long-duration deposition management.

## Japanese standard public wording

> A10-HSR / HSR v3.3は、超深海REY泥揚鉱を対象とするreduced-orderの前工学的輸送遷移frameworkである。本リポジトリは、論文補助アーカイブ、compact result summaries、selected source materials、およびdeferred-material inventoryを保存するものであり、実海域試験による検証、最終ライザー・ハードウェア設計、認証済み海洋実装、環境許認可準備完了、量産採鉱装置、商用展開準備完了、または長時間堆積管理の完成を主張するものではない。

## Publication decision

| Target | Decision |
|---|---|
| GitHub | GO after replacing GitHub account placeholder |
| Zenodo | GO after GitHub release |
| Jxiv | Conditional GO after manuscript metadata/title/conclusion wording review |
| Engineering deployment claim | NO |
| Sea-trial / field validation claim | NO |
| Safety / certification claim | NO |
| Commercial deployment claim | NO |

## Gate decision

`PASS-WITH-MINOR-PUBLICATION-FIXES-A10-HSR-PUBLIC-GATE-0`
