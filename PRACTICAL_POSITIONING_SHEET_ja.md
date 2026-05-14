# 実務位置づけシート

## 技術名

Hybrid Smart Riser v3.3 for Ultra-Deep REY Mud Lifting  
または  
A10-HSR REY Mud Lifting Transport-Transition Theory

## 簡易分類フォルダ名

超深海REY泥揚鉱・スマートライザー輸送遷移理論

## 対象産業

- 深海採鉱
- REY泥・レアアース泥揚鉱
- 海洋資源輸送
- スラリー輸送
- ライザー流体・構造連成
- フェイルセーフ輸送制御
- offshore pre-engineering assessment

## 現場課題

超深海REY泥揚鉱では、平均輸送性能だけでなく、blackout、coastdown、急減速、堆積成長、偏析、残差圧力スパイク、構造応答を同時に考える必要がある。  
「安全に止まれる」だけでは不十分であり、通常運転で持続輸送できるか、障害時にcollapseを診断できるか、長時間運転でdeposition creepが支配制約になるかを分ける必要がある。

## A10/HSRの役割

HSR v3.3は、通常運転相、障害遷移相、堆積制約相を同一のreduced-order model上で分類する。  
A10的には、単一性能最大化ではなく、相境界・診断量・安全遷移を整理するstructured pre-engineering frameworkとして機能する。

## 期待効果

- sustained-transport sectorとcollapse sectorを分離できる。
- no-blackout条件での臨界入口速度を評価できる。
- blackout/coastdownによるfault-induced collapseを診断できる。
- 長時間運転でdeposition creepが支配制約に移ることを検出できる。
- residual pressure-gradient derivativeを安全トリガ候補として評価できる。

## 検証済み範囲

論文内では以下が報告されている。

- no-blackout 40 s sweepで、臨界入口速度はおおよそ uss ≈ 1.25–1.26。
- uss=1.6 no-blackout条件では、40 s, 80 s, 120 sの全てで正の流速を維持。
- 堆積厚は40 sから120 sにかけて増加し、hydraulic arrestではなくdeposition creepが長時間制約となる。
- Ucrit再構成の保守側感度試験でも、uss=1.6 sustained resultは崩れない。
- residual pressure-gradient time derivativeは、2–10%ノイズ、1サンプル遅れ、5点平滑でもblackout検出率1.000、no-blackout偽陽性0.000。

## 未検証範囲

- 実海域試験
- full 3D CFD/FSI
- 実鉱泥・実粒度分布・実ライザー長での検証
- 長時間deposition management
- 実センサのドリフト・欠測・校正誤差
- offshore safety certification
- environmental permitting
- production deployment

## 実装への次ステップ

1. paper companion archiveとして公開する。
2. HSR v3.3を「採鉱装置」ではなく「前工学的輸送遷移理論」として説明する。
3. Ucritをclosure tableに明示列として追加する。
4. 長時間堆積管理則と再流動化戦略を追加検証する。
5. 実センサドリフト・欠測・バイアスに対する診断ロジックを検証する。
6. full 3D CFD/FSIおよび実海域条件への段階的接続を検討する。

## 想定読者

- 海洋資源工学研究者
- 深海採鉱技術者
- スラリー輸送研究者
- 流体・構造連成研究者
- フェイルセーフ制御研究者
- INPIT / 知財・共同研究相談員
- レアアース資源政策・技術評価関係者

## 誇張しない一文の結論

HSR v3.3は、超深海REY泥揚鉱における持続輸送相・障害遷移相・堆積制約相をreduced-order model上で区別する前工学的輸送遷移理論であり、実海域検証済み採鉱システムや最終ライザー設計を主張するものではない。
