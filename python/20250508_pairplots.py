import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# load csv file
df = pd.read_csv("data.csv")

# 特徴量リスト
features = ["A", "B", "C", "D", "E", "F", "G", "H"]

"""plotlyバージョン"""
# 2×4のサブプロットを作成
fig = make_subplots(
    row = 2, cols = 4, 
    subplot_titles = features,
    horizontal_spacing = 0.05, 
    vertical_spacing = 0.1
)

# 各特徴量の散布図を作成
for idx, feat in enumerate(features):
    # 相関係数
    r = df["x"].corr(df[feat])
    
    # subplotの位置を決定
    row = idx // 4 + 1
    col = idx % 4 + 1
    
    # 散布図を追加
    fig.add_trace(
        go.Scatter(
            x = df["x"],
            y = df[feat],
            mode = "markers",
            marker = dict(size = 5, opacity = 0.5),
            showlegend = True 
        ),
        row = row, col = col
    )
    
    # 相関係数を追記
    fig.add_anotation(
        text = f"r = {r:.3f}",
        xref = f"x{idx + 1} domain",
        yref = f"y{idx + 1} domain",
        x = 0.05, 
        y = 0.90
        xanchor = "left",
        yanchor = "top",
        showarrow = False,
        font = dixt(size = 12, color = "black"),
        row = row, col = col
    ) 
    
# レイアウトの設定
fig.update_layout(
    title = "Pairplots",
    height = 600,
    width = 1200,
    title_x = 0.5
)

# 保存
fig.save("20250508_pairplots_plotly.png")
fig.show()



"""seabornバージョン"""
# 1. 相関係数を描画用に計算・注釈する関数
def corrfunc(x, y, **kws):
    """散布図上にPearson相関係数を表示する関数"""
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate(f"r = {r:.2f}",
                xy = (0.1, 0.9), 
                xycoords = ax.transAxes,
                fontsize = 12)

# 2. PairGridの設定
g = sns.PairGrid(df, x_vars = ['X'], y_vars = features[1:], height = 3, aspect = 1)

# 3. 散布図のプロットと相関注釈のマッピング
g.map(sns.scatterplot)        
g.map_lower(corrfunc)      

# 4. レイアウト調整と表示
plt.tight_layout()
plt.savefig("20250508_pairplots_seaborn.png", dpi = 300)
plt.show()
