# 再現性メモ

## 現在の状態

このリポジトリはpaper companion archiveである。

## あるもの

- 論文PDF
- 論文表に対応するsummary CSV
- selected source scripts
- closure_table.csv
- selected figures
- logs / JSON
- claim boundary / limitations

## あえて入れていないもの

- venv
- site-packages
- compiled dependencies
- pycache
- monolithic code dump

## 再現性の意味

このサロゲートはreduced-order transport-transition modelであるため、再現性とは実海域試験との一致ではなく、同じscripts、closure table、input assumptions、scenario conditionsで同等のphase behaviorを再生成できることを意味する。

## 現在の限界

full 3D CFD/FSIや実海域検証の再現パッケージではない。


## v0.1.1-public-gate path normalization

元アップロード由来の日本語フォルダ名は、zip展開環境によって `#U....` 形式へエスケープされる場合がある。v0.1.1-public-gateでは、manifest検証性を優先し、selected source directoryをASCII-safeな `hsr_simulation_rareearth_selected` に正規化した。
