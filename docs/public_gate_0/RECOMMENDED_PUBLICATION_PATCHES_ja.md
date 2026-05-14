# 推奨公開前修正: A10-HSR v0.1.1-public-gate

## 反映済み

- versionを `v0.1.1-public-gate` に更新。
- README / README_jaにPUBLIC-GATE-0 noteを追加。
- CLAIM_BOUNDARY.mdに公開用標準文言を追加。
- selected source-material directoryをASCII-safe pathへ正規化。
- FILE_MANIFEST.csv/jsonをmanifest-excluding-self方式で再生成。
- MANIFEST_NOTE.mdを追加。
- tools/verify_manifest_excluding_self.pyを追加。
- docs/public_gate_0/に監査資料を同梱。
- CITATION.cffのrelease dateを2026-05-07に更新。

## GitHub公開前に手作業で必要

`CITATION.cff` の以下を実際のGitHubアカウント名に置換する。

```text
https://github.com/<your-github-account>/a10-hsr-rey-mud-lifting-transport-transition
```

## Zenodo後に必要

Zenodo DOI発行後、`CITATION.cff`の `doi: "pending"` を実DOIに置換する。

## Jxiv前に推奨

論文PDF本文では、`practically actionable`、`完成判定`、`完成域` という表現がある。現状でも限界節により境界は示されているが、Jxiv正式提出前には、`completion` を `pre-engineering readiness` へ弱めるか、タイトル・結論に非主張注記を追加することを推奨する。
