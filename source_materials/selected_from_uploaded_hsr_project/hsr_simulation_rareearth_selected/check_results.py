import csv
from collections import Counter

filename = "closure_table.csv"
status_counts = Counter()
tau_w_list = []
beta_eq_list = []
r_p_list = []

try:
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            status_counts[row['status']] += 1
            if row['status'] == 'ok':
                tau_w_list.append(float(row['tau_w_Pa']))
                beta_eq_list.append(float(row['beta_eq']))
                r_p_list.append(float(row['r_p_m']))

    print("="*40)
    print(" 🏆 HSR v2.1 クロージャーテーブル生成 最終リザルト 🏆")
    print("="*40)
    print(f"総ケース数: {sum(status_counts.values())} 件")
    for status, count in status_counts.items():
        print(f" - {status.upper()}: {count} 件")
    
    if tau_w_list:
        print("-" * 40)
        print("【抽出された物理パラメータの範囲】")
        print(f"壁面剪断応力 (tau_w): {min(tau_w_list):.2f} 〜 {max(tau_w_list):.2f} Pa")
        print(f"偏析振幅 (beta_eq)  : {min(beta_eq_list):.6f} 〜 {max(beta_eq_list):.6f}")
        print(f"プラグ半径 (r_p)    : {min(r_p_list):.4f} 〜 {max(r_p_list):.4f} m")
    print("="*40)

except FileNotFoundError:
    print(f"エラー: {filename} が見つかりません。同じフォルダで実行してください。")
