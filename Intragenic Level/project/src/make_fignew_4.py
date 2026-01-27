import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
out_part_date_3='output_part.csv'
savefile ="output_part.csv4.png"
pos = "bot"     #"bot" 底部	"left" 左侧 	"right" 右侧	"top" 顶部

current_vars=['out_part_date_3','savefile','pos']
try:
    # 加载TOML配置文件
    with open('config.toml', "r", encoding="utf-8") as f:
        config = toml.load(f)
    for var_name in current_vars:# 遍历所有需要检查的变量
        if var_name in config:
            globals()[var_name] = config[var_name]
except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")

date = pd.read_csv(
    out_part_date_3,
    sep='\t',
    header=0,
)
min_query_start = date['query_start'].min()
min_reference_start = date['reference_start'].min()
max_query_end = date['query_end'].max()
max_reference_end = date['reference_end'].max()
min_values = min(min_query_start, min_reference_start)  # 或者使用其中一个
max_values = max(max_query_end, max_reference_end)
max_dot_new = max_values - min_values
new_max_dot = round(max_dot_new / 1000000)
dot_c = date['byevents'].tolist()  # 提取 转换为列表-怕【
dot_red_1new = date['query_start'].tolist()
dot_red_2new = date['reference_start'].tolist()
min_val = min(dot_c)
max_val = max(dot_c)
num_bins = int(np.ceil(max_val) - np.floor(min_val)) + 1

rainbow_colors = ['#6A0DAD', '#1E90FF', '#00FA9A', '#FF8C00', '#FF0000']
positions = [0.0, 0.35, 0.5, 0.65, 1.0]

cmap = LinearSegmentedColormap.from_list('uneven_rainbow', list(zip(positions, rainbow_colors)))
norm = plt.Normalize(vmin=70, vmax=max_val)
result = [(a + b) / 2 - min_values for a, b in zip(dot_red_1new, dot_red_2new)]
resultn = [(a + b) / 2 - min_values - min(a, b) for a, b in zip(dot_red_1new, dot_red_2new)]

fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
if pos == "bot":  # 底部
    ax_dot = fig.add_axes([0.1, 0.1, 0.75, 0.7])  # 几乎全屏    左边开始 底部开始 宽度 高度
    ax = fig.add_axes([0.15, 0.7, 0.15, 0.15])

    xs = result
    ys = resultn
if pos == "left":  # 左侧
    ax_dot = fig.add_axes([0.15, 0.1, 0.75, 0.7])  # 几乎全屏    左边开始 底部开始 宽度 高度
    ax = fig.add_axes([0.7, 0.7, 0.15, 0.15])
    xs = resultn
    ys = result

if pos == "right":  # 右侧
    ax_dot = fig.add_axes([0.1, 0.1, 0.75, 0.7])  # 几乎全屏    左边开始 底部开始 宽度 高度
    ax = fig.add_axes([0.15, 0.7, 0.15, 0.15])
    xs = resultn
    ys = result
    ax_dot.invert_xaxis()

if pos == "top":  # 顶部
    ax_dot = fig.add_axes([0.1, 0.05, 0.75, 0.6])  # 几乎全屏    左边开始 底部开始 宽度 高度
    ax = fig.add_axes([0.15, 0.75, 0.15, 0.15])
    ax_dot.invert_yaxis()
    xs = result
    ys = resultn
# 左边开始 底部开始 宽度 高度
# 先创建大图，后创建小图 小图会覆盖在大图上面
ax.spines['top'].set_visible(False)  # 隐藏上边框
ax.spines['right'].set_visible(False)  # 隐藏右边框
ax_dot.spines['top'].set_visible(False)  # 隐藏上边框
ax_dot.spines['right'].set_visible(False)  # 隐藏右边框
ax_dot.spines['left'].set_visible(False)  # 隐藏左边框
ax_dot.spines['bottom'].set_visible(False)  # 隐藏下边框
# ax_dot.set_xticks([])  # 隐藏x轴刻度
ax_dot.set_yticks([])  # 隐藏y轴刻度
ax_dot.scatter(xs, ys,  # x轴数据  y轴数据
               c=dot_c, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
               norm=norm, s=1,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
               alpha=1, marker='o',  # 透明度    点的形状
               edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
               )  # 画出点
n, bins, patches = ax.hist(dot_c, bins=num_bins, edgecolor=None, alpha=1)
ax.set_xlim(left=70, right=max_val)
xticks_positions = np.arange(70, max_val + 5, 5)  # +5确保包含max_val
ax.set_xticks(xticks_positions)

ax.set_xlabel("%identity", fontsize=10, labelpad=8)
ax.set_ylabel("Quantity", fontsize=10, labelpad=8)

for i, patch in enumerate(patches):
    bin_center = (bins[i] + bins[i + 1]) / 2

    patch.set_facecolor(cmap(norm(bin_center)))
# ax_dot.set_title(f"Chr1: {int(min_values)} - {int(max_values)}", fontsize=17, pad=20)
ax_dot.text(0.5, -0.08, f"Chr1: {int(min_values)} - {int(max_values)}",
            transform=ax_dot.transAxes,
            fontsize=18,
            ha='center',
            va='bottom',
            style='italic')

extended_max = np.ceil(new_max_dot / 5) * 5
if pos == "bot" or pos == "top":
    ax3 = ax_dot.twiny()
    # 将刻度线和标签移动到下方
    # 简单修改1：设置ax_dot的范围为bp单位
    ax_dot.set_xlim(-1 * 1000000, extended_max * 1000000)  # 延长到最近的5的倍数
    # 简单修改2：设置ax3的范围为Mbp单位
    ax3.set_xlim(-1, extended_max)  # 延长到最近的5的倍数
    # 简单修改3：同步spine位置确保0对齐
    ax3_ticks = np.arange(0, extended_max + 5, 5)  # 生成刻度位置
    ax3_ticks = ax3_ticks[ax3_ticks < extended_max]  # 过滤超出范围的刻度
    ax3.set_xticks(ax3_ticks)  # 设置刻度位置
if pos == "right" or pos == "left":
    ax3 = ax_dot.twinx()
    ax_dot.set_ylim(-1 * 1000000, extended_max * 1000000)  # 延长到最近的5的倍数
    ax3.set_ylim(-1, extended_max)  # 延长到最近的5的倍数
    ax3_ticks = np.arange(0, extended_max + 5, 5)  # 生成刻度位置
    ax3_ticks = ax3_ticks[ax3_ticks < extended_max]  # 过滤超出范围的刻度
    ax3.set_yticks(ax3_ticks)  # 设置刻度位置

if pos == "bot":  # 底部

    ax3.spines['bottom'].set_position(('axes', 0.025))
    ax3.xaxis.set_ticks_position('bottom')  # 刻度线在下方
    ax3.xaxis.set_label_position('bottom')  # 标签在下方
    ax3.spines['top'].set_visible(False)  # 隐藏上边框
    ax3.spines['right'].set_visible(False)  # 隐藏右边框
    ax3.spines['left'].set_visible(False)  # 隐藏左边框
    ax3.text(extended_max, 0, 'Mbp', transform=ax3.transData,
             ha='left', va='top', fontsize=10)
if pos == "top":  # 顶部

    ax3.spines['top'].set_position(('axes', 0.98))
    ax3.xaxis.set_ticks_position('top')
    ax3.xaxis.set_label_position('top')
    # ax3.spines['top'].set_visible(False)  # 隐藏上边框
    ax3.spines['right'].set_visible(False)  # 隐藏右边框
    ax3.spines['left'].set_visible(False)  # 隐藏左边框
    ax3.spines['bottom'].set_visible(False)  # 隐藏下边框
    ax3.text(extended_max, 0, 'Mbp', transform=ax3.transData,
             ha='left', va='bottom', fontsize=10)  # x 坐标   y 坐标

if pos == "right":
    ax3.spines['right'].set_position(('axes', 0.98))
    ax3.yaxis.set_ticks_position('right')
    ax3.yaxis.set_label_position('right')
    ax3.spines['top'].set_visible(False)  # 隐藏上边框
    # ax3.spines['right'].set_visible(False)  # 隐藏右边框
    ax3.spines['left'].set_visible(False)  # 隐藏左边框
    ax3.spines['bottom'].set_visible(False)  # 隐藏下边框
    ax3.text(extended_max * 0.5, extended_max, 'Mbp', transform=ax3.transData,
             ha='left', va='bottom', fontsize=10, rotation=90)

if pos == "left":
    ax3.spines['left'].set_position(('axes', 0.02))
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')
    ax3.spines['top'].set_visible(False)  # 隐藏上边框
    ax3.spines['right'].set_visible(False)  # 隐藏右边框
    # ax3.spines['left'].set_visible(False)  # 隐藏左边框
    ax3.spines['bottom'].set_visible(False)  # 隐藏下边框
    ax3.text(0, extended_max, 'Mbp', transform=ax3.transData,
             ha='right', va='bottom', fontsize=10, rotation=-90)

# 移除默认的刻度标签
ax_dot.set_xticks([])
ax_dot.set_yticks([])


plt.savefig(savefile, dpi=1000, format='png', bbox_inches='tight')
plt.close()
# plt.show()