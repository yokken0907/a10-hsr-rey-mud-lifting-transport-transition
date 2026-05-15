# A10-HSR REY Mud Lifting Transport-Transition Surrogate

このリポジトリは、以下のAI支援独立研究論文に対応するGitHub配置用フォルダである。

**超深海 REY 泥揚鉱のための Hybrid Smart Riser v3.3  
実務接続可能性の完成判定に向けた輸送遷移理論の再構成と検証**

英題: **Hybrid Smart Riser v3.3 for Ultra-Deep REY Mud Lifting: Toward a Practically Actionable Transport-Transition Theory**

著者: 吉村圭司（Independent Researcher）  
状態: GitHub-ready paper companion archive v0.1.1-public-gate

## 簡易分類フォルダ名

**超深海REY泥揚鉱・スマートライザー輸送遷移理論**

## 位置づけ

HSR v3.3は、超深海REY泥揚鉱を対象とする reduced-order pre-engineering transport-transition theory である。  
通常揚水下の sustained-transport sector、blackout/coastdown下の collapse sector、長時間運転における deposition-limited horizon を同一の枠組みで区別する。

## 中心的解釈

この理論は「安全に止まれるか」だけではなく、「通常運転相・障害遷移相・堆積制約相を識別できるか」を扱う。  
したがって、単なるフェイルセーフ停止理論ではなく、実務設計へ接続可能な前工学的な輸送遷移理論として扱うのが安全である。

## 技術的ビジュアル案内

初めて本リポジトリを見る技術的関心のある読者向けに、ブラウザだけで開ける技術的ビジュアル案内ページを同梱しています。

`docs/technical_visual_orientation/index.html`

このページは、A10-HSR REY mud-lifting transport-transition の構造をプロジェクト固有の観点から整理する補助資料です。本リポジトリにおける mission variable は commercial extraction yield、seabed mining deployment、または certified riser engineering ではなく、continuous transport、deposition-limited loss、blackout/coastdown collapse、restart fragility、pump/flow-margin exhaustion を reduced lifting/transport surrogate 上で識別する sustained slurry-transport viability under transition risk です。

また、このページでは reduced transport-surrogate state channels、critical transport velocity、solids loading、deposition thresholds、blackout/coastdown and restart regimes、evidence hierarchy、リポジトリ閲覧順、および claim boundary を短く整理しています。主要な図解セクションには replay control を付けており、静的テンプレートではなく診断ロジックを段階的に確認できます。

このページは説明補助であり、offshore engineering simulation を実行するものではありません。海底採鉱システム、riser / pump-system design、operational / restart procedure、environmental clearance、安全認証、または商用抽出可能性を示すものでもなく、論文本体、source materials、figures、または専門家による独立評価を置き換えるものでもありません。

## 含まれるもの

- 論文PDF
- README日本語版・英語版
- claim boundary
- limitations
- AI支援開示
- 実務位置づけ表
- アップロードされたHSR project archiveから抽出したselected scripts / figures / logs / JSON / closure CSV
- 論文表から再構成したsummary CSV
- virtual environment / binary dependenciesを除外したclassification index

## 主張しないこと

本リポジトリは、以下を主張しない。

- 実海域導入済み採鉱技術
- full 3D offshore validation
- sea-trial validation
- production-ready mining equipment
- final riser design
- environmental permitting readiness
- industrial safety certification
- 長時間堆積管理の完成

## 現在の状態

これはpaper companion archiveであり、海洋鉱業用の認証済み実装設計、実海域検証済みシステム、または最終操業手順ではない。

## PUBLIC-GATE-0 note

このv0.1.1-public-gate版では、公開前監査に基づき、保守的なclaim boundaryを再固定している。本リポジトリはreduced-orderの前工学的輸送遷移アーカイブとしてのみ説明し、実海域検証済み採鉱システム、最終ライザー設計、認証済み海洋実装、環境許認可準備完了、量産採鉱装置、商用展開準備完了、または長時間堆積管理の完成を主張しない。
