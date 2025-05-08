import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# load csv
df = pd.read_csv("20250508_pairplots.csv")

# 特徴量リスト
features = ["A", "B", "C", "D", "E", "F", "G", "H"]

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
        xref = f"x{idx + 1}",
        yref = f"y{idx + 1}",
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
fig.save("20250508_pairplots.png")
fig.show()