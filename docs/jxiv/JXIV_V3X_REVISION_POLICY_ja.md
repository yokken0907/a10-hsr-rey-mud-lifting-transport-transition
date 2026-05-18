# Jxiv v3.0 / v3.x 改訂方針メモ

## 目的

このメモは、HSR v3.x論文をJxiv上で改訂する場合の方針を整理する。  
目的は論文を増やすことではなく、既存理論の主張境界を明確化し、実務価値を誤読されにくくすることである。

## 基本方針

HSR v3.xは、超深海REY泥揚鉱を対象とするreduced-orderの前工学的輸送遷移・安全診断フレームワークであり、実海域検証済み採鉱装置、最終ライザー設計、安全認証済み技術、または商用化準備完了技術を主張するものではない。

Jxiv改訂では、理論の価値を弱めるのではなく、以下を明確化する。

- 実海域検証済み装置ではない。
- 最終ライザー設計ではない。
- 安全認証済み技術ではない。
- 商用化準備完了ではない。
- reduced-order modelによる前工学的診断・候補選別・検証計画整理である。
- 現場価値は「完成」ではなく「次に検証すべき危険相・診断量・堆積制約を示すこと」にある。

## タイトル方針

現行タイトルに「完成判定」「practically actionable」「completion」などの語が含まれる場合、読者に完成装置・実用準備完了と誤読される可能性がある。

### 推奨日本語タイトル案

> 超深海REY泥揚鉱のためのHybrid Smart Riser v3.x：  
> 前工学的評価に向けた輸送遷移・ブラックアウト診断・堆積制約フレームワーク

### 推奨英語タイトル案

> Hybrid Smart Riser v3.x for Ultra-Deep REY Mud Lifting:  
> A Pre-Engineering Framework for Transport Transition, Blackout Diagnostics, and Deposition Constraints

## サブタイトル方針

### 日本語候補

> 実海域検証済み採鉱装置ではなく、reduced-order modelに基づく安全診断・候補選別理論としての位置づけ

### 英語候補

> A reduced-order diagnostic and screening framework, not a field-validated mining system or final riser design

## 抄録改訂方針

抄録では、以下の順序にすると安全である。

1. 背景：超深海REY泥揚鉱では持続輸送、障害時崩壊、堆積制約が問題になる。
2. 方法：reduced-order / pre-engineering surrogateを用いる。
3. 結果：no-blackout transition、blackout diagnostic candidate、deposition-limited horizonを報告する。
4. 限界：実海域試験、3D CFD/FSI、安全認証、最終設計ではない。
5. 価値：次段階検証へ進むための診断・候補選別フレームワークである。

## 置き換え推奨語

| 避けたい表現 | 推奨表現 |
|---|---|
| 実務接続可能性の完成判定 | 前工学的評価に向けた整理 |
| 完成域 | reduced-order model上の余裕域 |
| 実用可能性を示した | 実務検証へ接続可能な候補条件を整理した |
| practically actionable | pre-engineering actionable / validation-oriented |
| completion | framework consolidation / validation roadmap |
| safety trigger | candidate diagnostic trigger |
| mining system | transport-transition surrogate |
| operational procedure | diagnostic logic / validation target |

## 結論部の推奨文

> 本研究は、超深海REY泥揚鉱における輸送遷移、blackout/coastdown collapse、およびdeposition-limited horizonをreduced-order model上で整理した前工学的フレームワークである。  
> これは実海域検証済み採鉱装置、最終ライザー設計、安全認証済み技術、または商用化準備完了技術を主張するものではない。  
> 本研究の実務的意義は、実鉱泥試験、CFD/FSI、センサ診断検証、長時間堆積管理、および専門家レビューへ進むための候補条件と検証課題を明確化した点にある。

## Jxiv改訂で追加すべき節

可能であれば、以下の短い節を追加する。

### Practical Interpretation and Claim Boundary

内容：

- 本研究の実務的読み方
- 現時点で主張しないこと
- 次段階検証で必要なこと
- GitHubリポジトリ内の `CLAIM_BOUNDARY.md`、`FIELD_VALUE_ja.md`、`NEXT_VALIDATION_PLAN.md` への参照

## 改訂の優先順位

1. タイトル・抄録・結論の過強表現を弱める。
2. `pre-engineering`、`reduced-order`、`diagnostic candidate` を明示する。
3. 実海域検証・最終設計・認証済み技術ではないことを明記する。
4. GitHub側の field-value / validation-plan 文書との整合性を取る。
5. 数値結果そのものは変えず、解釈境界を保守化する。

## 最終方針

Jxiv v3.x改訂は、新規理論の追加ではなく、**同じ理論をより安全に、より専門家が評価しやすい形へ再配置する改訂**として行う。
