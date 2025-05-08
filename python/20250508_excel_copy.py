import os
import shutil
import pandas as pd

"""ファイルコピー操作

original ディレクトリ内のすべてのサブディレクトリを再帰的に探索し、ファイル名に "knock" を含み拡張子が .xlsx のファイルだけを新しく作成する同階層の copy ディレクトリにまとめてコピーする。

"""

# ベースとなるディレクトリ名
BASE_DIR = 'original'
# コピー先ディレクトリ名
DEST_DIR = 'copy'

# カレントディレクトリがスクリプトと同じ階層であることを前提
# copy ディレクトリがなければ作成
if not os.path.exists(DEST_DIR):
    os.makedirs(DEST_DIR)

# original 以下を再帰的に探索
for root, dirs, files in os.walk(BASE_DIR):
    for filename in files:
        # ファイル名に "knock" を含み、拡張子が .xlsx ならコピー対象
        if 'knock' in filename and filename.lower().endswith('.xlsx'):
            src_path = os.path.join(root, filename)
            dest_path = os.path.join(DEST_DIR, filename)
            # 同名ファイルは上書き
            shutil.copy(src_path, dest_path)
            print(f'Copied: {src_path} -> {dest_path}')


"""Excelデータ操作

copy ディレクトリ内のすべての Excel ファイルについて、それぞれのファイルの Sheet1 から

C2:D11 → 新しいシートの B,C 列

I2:K6 → 新しいシートの F,G,H 列

を順番に縦に連結して書き出す。

"""
# 処理対象ディレクトリ
input_dir = 'copy'
# 出力ファイル名
output_file = 'outputs.xlsx'
# 出力用の DataFrame を初期化
# B,C 列用と F,G,H 列用にそれぞれカラム名を設定
df_result = pd.DataFrame(columns=['B', 'C', 'F', 'G', 'H'])

# ディレクトリ内の .xlsx ファイルをループ
for fname in os.listdir(input_dir):
    if not fname.lower().endswith('.xlsx'):
        continue
    path = os.path.join(input_dir, fname)
    # sheet_name='sheet1' を明示。C2:D11 は列 C,D を 2 行目から 11 行目まで
    df_cd = pd.read_excel(
        path,
        sheet_name='Sheet1',
        usecols='C:D',
        skiprows=0,     # 1 行ヘッダを飛ばして、2 行目から
        nrows=10        # 2～11 行目の計 10 行
    )
    # I2:K6 は列 I,J,K を 2 行目から 6 行目まで
    df_ijk = pd.read_excel(
        path,
        sheet_name='Sheet1',
        usecols='I:K',
        skiprows=0,     # 1 行ヘッダを飛ばして
        nrows=5         # 2～6 行目の計 5 行
    )

    # 列名を B,C および F,G,H にリネーム
    df_cd.columns = ['B', 'C']
    df_ijk.columns = ['F', 'G', 'H']

    # 行数をそろえるため、足りない行には NaN が入る
    df_tmp = pd.concat([df_cd.reset_index(drop=True),
                        df_ijk.reset_index(drop=True)],
                       axis=1)

    # 結果に追加
    df_result = pd.concat([df_result, df_tmp], ignore_index=True)

# 結果を新しい Excel ファイルに書き出し
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    df_result.to_excel(writer, sheet_name='Sheet1', index=False)

print(f'処理完了。{output_file} を出力しました。')
