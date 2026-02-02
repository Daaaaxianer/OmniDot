import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
import numpy as np
import toml
import statistics
from matplotlib.colors import LinearSegmentedColormap, Normalize
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
import click
"""
gene_show
"chr","tnew_name"
  1	  coca01g00001

gff: six column
"chr","tnew_name","tpos_start","tpos_end","tchain","torder","told_name"
  1	  coca01g00001	   390	       1820	      -         1	 Cc01_g00010

lens: three line
"chr",   "tpos",     "torder"
  1       23021474     1399

color_dates
"tnew_name"   , "color"
coca01g00001    #F09201

Collinearity
# Alignment 1: score=1451 pvalue=0.2214 N=49 1&1 plus
vivi01g00514 514 vivi01g00515 515 1

KS
id1	id2	ka_NG86	ks_NG86	ka_YN00	ks_YN00
vivi01g00514	vivi01g00515	0.0686	0.1182	0.0678	0.1258
vivi01g00517	vivi01g00516	0.0796	0.137	0.0835	0.1197
vivi01g00518	vivi01g00517	0.0278	0.0306	0.0299	0.0248
"""
mods = "blast"  # "blast" "collinearity" 画blast 还是 collinearity
mod = "torder"  # "torder"或"tpos"根据gff文件的 tpos 或 torder 画图 "torder" 根据基因在染色体总的基因上的相对顺序 "tpos" 根据基因在染色体上的真实位置
name_1 = "NAN"  # 第1列的物种名称
name_2 = "NAN"  # 第2列的物种名称
name_3 = "NAN"#在mod_show=rate下需要填写
name_4 = "NAN"#在mod_show=rate下需要填写
collinearity_ks_mod = "False"  # "False" "True"   是否输入Collinearity_ks文件  False情况下只显示Collinearity点图无ks
Ks = "False"  # "False" "True"是否显示ks的值      collinearity情况下
ks_mod = "Median"  # "Average" 或 "Median"   "Average"显示ks的平均数"Median"显示ks的中位数       collinearity情况下
name_id = "False"  # "True" "False" 是否去除物种自身比较时的对角线 (当并非是物种自身比较时 应为  "False"  )   collinearity情况下
N = 5  # Alignment  ＞  几   就可以合并数据    collinearity情况下
NSD = 10  # Alignment  ＞  几   就可以显示数据     collinearity情况下
Distance = 250  # 合并数据的距离            collinearity情况下
half = "False"  # "self_col_2" "False" "True" "self_blco_2"  "False" "True"是否单个物种为自身比对且画一半"self_col_2"为2个物种的自身比对且画一半,并合并 "self_blco_2"collinearity与blast分别画
block = "False"  # 在collinearity 与  "self_col_2"情况下  "False" "True"  是否显示 block 的大小
top_block = "False"  # block情况下  "False" "True"  是否突出显示最近的一次加倍 需要ks文件
ks_dem_1 = 1.5  # Collinearity文件 top_block"True"情况下 分割最近的一次加倍的ks
ks_dem_1_1=0.8
ks_dem_2 = 1.8  # Collinearity_2文件 top_block"True"情况下 分割最近的一次加倍的ks
ks_dem_2_2= 1.2
pos = "bot_left"  # 显示在"bot_left" 左下   "bot_right" 右下  "top_left" 左上   "top_right"   右上
mod_show = "nan"  # "rate" "tree" "density" "nan"        collinearity情况下"rate" 显示基因的比例  "tree"  显示特定基因的名称  "density"显示基因的密度  "nan" 都不显示   blast情况下 "tree" "density"
bin_sizes = 20  # "density"显示基因的密度  定义分箱大小  默认单位为百万
color_date = "False"  # "False" "True"是否有显示基因名称的连接线的颜色文件       collinearity情况下     tree情况下
color_dates = 'NAN'  # 写入"tree"下显示的基因名称对应颜色的文件路径     tree情况下
rcolor_dates = 'NAN'
wide = 400  # 显示的名称上下的距离差         collinearity情况下     tree情况下
longs = 350  # 显示的名称左右的距离差         collinearity情况下     tree情况下
name_1_4 = "NAN"  # 第1列的物种名称
name_2_4 = "NAN"  # 第2列的物种名称
posd = "bot_left"  # "bot_left"   "bot_right"
Collinearity_ks = "NAN"
Collinearity = "NAN"
blast_name = "NAN"
blast_gff_1 ="NAN"
blast_gff_2 = "NAN"
blast_lens_1 = "NAN"
blast_lens_2 = "NAN"
Collinearity_2 = "NAN"
Collinearity_ks_2 = "NAN"
blast_name_2 ="NAN"
blast_gff_1_2 = "NAN"
blast_gff_2_2 = "NAN"
blast_lens_1_2 ="NAN"
blast_lens_2_2 = "NAN"
gene_show = "NAN"
savefile='NAN'
current_vars=['color_dates','gene_show','Collinearity_ks','Collinearity','blast_name','blast_gff_1','blast_gff_2','blast_lens_1','blast_lens_2','mods','mod',
              'name_1','name_2','name_3','name_4','collinearity_ks_mod','Ks','ks_mod','name_id','N','NSD','Distance','half','block','top_block','ks_dem_1','ks_dem_1_1','ks_dem_2',
              'ks_dem_2_2','pos','mod_show','bin_sizes','color_date','wide','longs','name_1_4','name_2_4','Collinearity_2','Collinearity_ks_2','blast_name_2','blast_gff_1_2',
              'blast_gff_2_2','blast_lens_1_2','blast_lens_2_2','posd','savefile']
config_1 = "Chromosomal"    #3种模式 "Chromosomal" "Fragment" "Intragenic"
config_2 = "blast"          # 在"Chromosomal" 下的哪一个板块 "blast","Collinearity","blast_blast","Collinearity_Collinearity","blast_Collinearity"

current_vard=['config_1','config_2']
try:
    # 加载TOML配置文件
    with open('../config.toml', "r", encoding="utf-8") as f:
        config = toml.load(f)
    if "Config" in config:
        intragenic_config = config["Config"]
        for var_name in current_vard:# 遍历所有需要检查的变量
            if var_name in intragenic_config:
                globals()[var_name] = intragenic_config[var_name]
except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")

try:
    # 加载TOML配置文件
    with open('config.toml', "r", encoding="utf-8") as f:
        config = toml.load(f)
    updated = 0# 统计更新情况
    missing = 0
    if config_2 in config:
        new_config = config[config_2]
    else:
        print(f"❌ 错误: 配置文件加载失败 ")
    for var_name in current_vars:# 遍历所有需要检查的变量
        if var_name in new_config:# 更新全局变量
            globals()[var_name] = new_config[var_name]
            updated += 1
        else:
            missing += 1
except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")
dot_mod_1 = dict(c='grey', s=0.7, alpha=0.8, marker='o', edgecolors=None, linewidths=0)  # 定义画blast图得到的点样式
dot_mod_2 = dict(c='blue', s=0.7, alpha=0.8, marker='o', edgecolors=None, linewidths=0)  # 定义画blast图得到的点样式
dot_mod_3 = dict(c='red', s=0.7, alpha=0.8, marker='o', edgecolors=None, linewidths=0)  # 定义画blast图得到的点样式
class AlignmentResult:  # "torder" x     "torder" y
    def __init__(self, alignment, ns, red_dict_x_3, red_dict_y_3,
                 red_dict_x_4, red_dict_y_4, ks_1, ks_2,  # "tpos" x    "tpos" y      平均数  中位数
                 x_start1, x_end1, y_start1, y_end1,  # "torder"
                 x_start2, x_end2, y_start2, y_end2, ks  # "tpos"
                 ):
        self.alignment = str(alignment)
        self.ns = int(ns)
        self.x1 = float(red_dict_x_3)
        self.y1 = float(red_dict_y_3)
        self.x2 = float(red_dict_x_4)
        self.y2 = float(red_dict_y_4)
        self.ks1 = float(ks_1)
        self.ks2 = float(ks_2)
        self.x_start1 = float(x_start1)
        self.x_end1 = float(x_end1)
        self.y_start1 = float(y_start1)
        self.y_end1 = float(y_end1)
        self.x_start2 = float(x_start2)
        self.x_end2 = float(x_end2)
        self.y_start2 = float(y_start2)
        self.y_end2 = float(y_end2)
        self.ks = ks
class gene_id:
    def __init__(self, gene_name,
                 xy,
                 genes,
                 xys,
                 date_xy
                 ):
        self.gene_name = str(gene_name)
        self.xy = int(xy)
        self.genes = genes
        self.xys = xys
        self.date_xy = date_xy
# 画图
# 设置全局字体为 Times New Roman 将 Matplotlib 的默认字体设置为 Times New Roman，影响所有文本
mpl.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'custom'  # 使用自定义数学字体
plt.rcParams['mathtext.rm'] = 'Times New Roman'  # 常规数学字体
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'  # 斜体数学字体
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'  # 加粗数学字体
# 所有数学公式（如 $x$、$\alpha$）会使用 Times New Roman 字体渲染。
def read_gff(gff):
    new_gffs = pd.read_csv(
        gff,
        sep='\t',
        header=None,
        names=["chr", "tnew_name", "tpos_start", "tpos_end", "tchain", "torder", "told_name"],
        dtype={
            "chr": int,
            "tnew_name": str,
            "tpos_start": int,
            "tpos_end": int,
            "tchain": str,
            "torder": int,
            "told_name": str
        }
    )
    return new_gffs
def read_lens(lens):
    new_lens = pd.read_csv(
        lens,
        sep='\t',
        header=None,
        names=["chr", "tpos", "torder"],
        dtype={
            "chr": int,
            "tpos": int,
            "torder": int
        }
    )
    return new_lens
def read_blast(blasts,gff_1s,gff_2s,gff_chr_1s,chr_lens_1s,gff_chr_2s,chr_lens_2s):
    new_blasts = {}
    new_blast=[]
    blastsd = pd.read_csv(
        blasts,
        sep='\t',
        header=None,
        names=["query_id", "subject_id", "percent_identity", "alignment_length","mismatches",
               "gap_opens", "query_start", "query_end", "subject_start","subject_end", "evalue", "bit_score"],
        dtype={
            "query_id": str,
            "subject_id": str,
            "percent_identity": float,
            "alignment_length": int,
            "mismatches": int,
            "gap_opens": int,
            "query_start": int,
            "query_end": int,
            "subject_start": int,
            "subject_end": int,
            "evalue": float,
            "bit_score": float
        }
    )
    for blast in blastsd.itertuples(index=False):
        if not (blast.bit_score >= 100 and blast.evalue < 1e-5):
            continue
        if not (blast.query_id in gff_1s['tnew_name'].values and blast.subject_id in gff_2s['tnew_name'].values):# 检查query_id和subject_id是否在对应的gff文件的tnew_name列中
            continue
        if gff_chr_1s.get(blast.query_id) not in chr_lens_1s or gff_chr_2s.get(blast.subject_id) not in chr_lens_2s: # 确保基因位于两个物种lens中的染色体上
            continue
        if blast.query_id not in new_blasts: # 如果query_id不在new_blasts中，则初始化列表
            new_blasts[blast.query_id] = []
        if len(new_blasts[blast.query_id]) < 5:# 只保留前5个匹配
            new_blasts[blast.query_id].append(blast.subject_id)
        else:
            continue
        new_blast.append(blast)
    new_blasts_x = pd.DataFrame(new_blast,
                                columns=["query_id", "subject_id", "percent_identity", "alignment_length","mismatches", "gap_opens",
                                         "query_start", "query_end", "subject_start","subject_end", "evalue", "bit_score"])
    return new_blasts_x
def bla_dot(dot_x,gff_chr_x1,gff_torder_x1,chr_torder_x1,gff_chr_x2,gff_torder_x2,chr_torder_x2):
    dot_red_1s = []
    dot_red_2s = []
    for red in dot_x.itertuples():
        chr_red_1 = gff_chr_x1[red.query_id]
        dict_red_1 = gff_torder_x1[red.query_id]
        red_dict_x = chr_torder_x1[chr_red_1] + dict_red_1
        dot_red_1s.append(red_dict_x)# 得到红色x
        chr_red_2 = gff_chr_x2[red.subject_id]
        dict_red_2 = gff_torder_x2[red.subject_id]
        red_dict_y = chr_torder_x2[chr_red_2] + dict_red_2
        dot_red_2s.append(red_dict_y)# 得到红色y
    return dot_red_1s,dot_red_2s
def read_collinearity(Collinearity):
    blast_ld = pd.read_csv(
        Collinearity,
        sep=' ',
        header=None,
        names=["id1", "n", "id2", "x", "c", "u", "i", "t"],
    )
    return blast_ld
def order_col_dot(blast_ld,gff_chr_x1,gff_chr_x2,gff_torder_x1,gff_torder_x2,chr_torder_x1,chr_torder_x2,name_id):
    dot_len_1s = []
    dot_star_1s = []
    dot_end_1s = []
    dot_red_1s = []
    dot_red_2s = []
    doc_1s = []
    doc_2s = []
    iiss = 0
    for index, red in blast_ld.iterrows():
        if red.id1 == "#":
            if iiss == 0:
                dot_red_1s.extend(doc_1s)
                dot_red_2s.extend(doc_2s)
                if len(doc_1s) != 0:
                    dot_star_1s.append(doc_1s[0])
                    dot_end_1s.append(doc_1s[-1])
                    dot_len_1s.append(len(doc_1s))
            doc_1s = []  # dot_red_1=[]
            doc_2s = []  # dot_red_2=[]
            iiss = 0
            continue
        if name_id == "True":  # 根据  name_id  的赋值决定是否去除一条直线（自身与自身比对）的区块
            if gff_chr_x1[red["id1"]] == gff_chr_x2[red["id2"]]:
                if (abs(gff_torder_x1[red["id1"]] - gff_torder_x2[red["id2"]]) < 150):
                    iiss = 1
        chr_red_1 = gff_chr_x1[red["id1"]]
        dict_red_1 = gff_torder_x1[red["id1"]]
        red_dict_x = chr_torder_x1[chr_red_1] + dict_red_1
        doc_1s.append(red_dict_x)# 得到红色x
        chr_red_2 = gff_chr_x2[red["id2"]]
        dict_red_2 = gff_torder_x2[red["id2"]]
        red_dict_y = chr_torder_x2[chr_red_2] + dict_red_2
        doc_2s.append(red_dict_y) # 得到红色y
    if iiss == 0:
        dot_red_1s.extend(doc_1s)
        dot_red_2s.extend(doc_2s)
        if len(doc_1s) != 0:
            dot_star_1s.append(doc_1s[0])
            dot_end_1s.append(doc_1s[-1])
            dot_len_1s.append(len(doc_1s))
    return dot_red_1s,dot_red_2s,dot_star_1s,dot_end_1s,dot_len_1s
def tpos_col_dot(blast_ld,gff_chr_x1,gff_chr_x2,gff_torder_x1,gff_torder_x2,gff_dict_x1,chr_dict_x1,gff_dict_x2,chr_dict_x2,name_id):
    dot_len_1s = []
    dot_star_1s = []
    dot_end_1s = []
    dot_red_1s = []
    dot_red_2s = []
    doc_1s = []
    doc_2s = []
    iiss = 0
    for index, red in blast_ld.iterrows():
        if red.id1 == "#":
            if iiss == 0:
                dot_red_1s.extend(doc_1s)
                dot_red_2s.extend(doc_2s)
                if len(doc_1s) != 0:
                    dot_star_1s.append(doc_1s[0])
                    dot_end_1s.append(doc_1s[-1])
                    dot_len_1s.append(len(doc_1s))
            doc_1s = []  # dot_red_1=[]
            doc_2s = []  # dot_red_2=[]
            iiss = 0
            continue
        if name_id == "True":  # 根据  name_id  的赋值决定是否去除一条直线（自身与自身比对）的区块
            if gff_chr_x1[red["id1"]] == gff_chr_x2[red["id2"]]:
                if (abs(gff_torder_x1[red["id1"]] - gff_torder_x2[red["id2"]]) < 150):
                    iiss = 1
        chr_red_1 = gff_chr_x1[red["id1"]]
        dict_red_1 = gff_dict_x1[red["id1"]]
        red_dict_x = chr_dict_x1[chr_red_1] + dict_red_1
        doc_1s.append(red_dict_x)# 得到红色x
        chr_red_2 = gff_chr_x2[red["id2"]]
        dict_red_2 = gff_dict_x2[red["id2"]]
        red_dict_y = chr_dict_x2[chr_red_2] + dict_red_2
        doc_2s.append(red_dict_y)# 得到红色y
    if iiss == 0:
        dot_red_1s.extend(doc_1s)
        dot_red_2s.extend(doc_2s)
        if len(doc_1s) != 0:
            dot_star_1s.append(doc_1s[0])
            dot_end_1s.append(doc_1s[-1])
            dot_len_1s.append(len(doc_1s))
    return dot_red_1s, dot_red_2s, dot_star_1s, dot_end_1s, dot_len_1s
def read_ks(Collinearity_ksd,blast_ld):
    collinearity = pd.read_csv(
        Collinearity_ksd,
        sep='\t',
        dtype={
            "id1": str,
            "id2": str,
            "ka_NG86": float,
            "ks_NG86": float,
            "ka_YN00": float,
            "ks_YN00": float,
        }
    )
    collinearity_reverse = collinearity.copy()# 创建反向匹配的数据副本
    collinearity_reverse.columns = ["id2", "id1", "ka_NG86", "ks_NG86", "ka_YN00","ks_YN00"]  # 将列名中的id1和id2交换，使得反向的基因对（id2-id1）也能被表示。
    collinearity_combined = pd.concat([collinearity, collinearity_reverse],ignore_index=True)  # 将正向和反向的基因对合并为一个新表collinearity_combined，方便后续匹配。
    result = pd.merge(      # 先尝试正向合并
        blast_ld,  # 左表（主表），保留所有行。
        collinearity[['id1', 'id2', 'ks_NG86']],  # 右表，仅提取需要的列。
        on=['id1', 'id2'],  # 按照id1和id2列匹配。
        how='left'  # 左连接，保留blast_l的所有行，未匹配的ks_NG86填充为NaN。
    )
    unmatched = result[result['ks_NG86'].isna() & ~blast_ld['id1'].str.startswith('#')]# 找出未匹配的行
    # result['ks_NG86'].isna()：ks_NG86为NaN的行。# ~blast_l['id1'].str.startswith('#')：排除以#开头的注释行。
    # 对这些未匹配的行尝试反向匹配
    for idx, row in unmatched.iterrows():
        reverse_match = collinearity_combined[
            (collinearity_combined['id1'] == row['id2']) &
            (collinearity_combined['id2'] == row['id1'])
            ]
        if not reverse_match.empty:
            result.at[idx, 'ks_NG86'] = reverse_match['ks_NG86'].values[0]
    result.loc[blast_ld['id1'].str.startswith('#'), 'ks_NG86'] = None# 确保注释行（以#开头）的ks_NG86为None，避免干扰数据分析。
    final_result = result[['id1', 'n', 'id2', 'x', 'c', 'u', 'i', 't', 'ks_NG86']]# 最终结果
    return final_result
def order_ks_dot(final_result,gff_chr_x1,gff_chr_x2,gff_torder_x1,gff_torder_x2,chr_torder_x1,chr_torder_x2,top_block,ks_mod,ks_dem_1,ks_dem_1_1,name_id):
    dot_block_1,dot_block_2,dot_block_ks_1,dot_len_1,dot_star_1,dot_end_1,dot_red_1,dot_red_2,dot_c ,doc,doc_1,doc_2= [],[],[],[],[],[],[],[],[],[],[],[]
    iiss = 0
    for index, red in final_result.iterrows():
        if red.id1 == "#":
            if iiss == 0:
                dot_red_1.extend(doc_1)
                dot_red_2.extend(doc_2)
                dot_c.extend(doc)
                if top_block == "True":
                    # "Average"平均数"Median"中位数
                    if ks_mod == "Average":
                        if len(dot_block_ks_1) != 0:
                            average_block = sum(dot_block_ks_1) / len(dot_block_ks_1)  # 平均数
                            if ks_dem_1_1 <= average_block <= ks_dem_1:
                                dot_block_1.extend(doc_1)
                                dot_block_2.extend(doc_2)
                    if ks_mod == "Median":
                        if len(dot_block_ks_1) != 0:
                            median_block = statistics.median(dot_block_ks_1)  # 中位数
                            if ks_dem_1_1 <= median_block <= ks_dem_1:
                                dot_block_1.extend(doc_1)
                                dot_block_2.extend(doc_2)
                if len(doc_1) != 0:
                    dot_star_1.append(doc_1[0])
                    dot_end_1.append(doc_1[-1])
                    dot_len_1.append(len(doc_1))
            doc = []  # dot_c=[]
            doc_1 = []  # dot_red_1=[]
            doc_2 = []  # dot_red_2=[]
            dot_block_ks_1 = []
            iiss = 0
            continue
        if name_id == "True":
            if gff_chr_x1[red["id1"]] == gff_chr_x2[red["id2"]]:
                if (abs(gff_torder_x1[red["id1"]] - gff_torder_x2[red["id2"]]) < 150):
                    iiss = 1
        if red['ks_NG86'] == -1:
            continue
        if top_block == "True":
            dot_block_ks_1.append(red['ks_NG86'])
        chr_red_1 = gff_chr_x1[red["id1"]]
        dict_red_1 = gff_torder_x1[red["id1"]]
        red_dict_x = chr_torder_x1[chr_red_1] + dict_red_1
        doc_1.append(red_dict_x)# 得到红色x
        chr_red_2 = gff_chr_x2[red["id2"]]
        dict_red_2 = gff_torder_x2[red["id2"]]
        red_dict_y = chr_torder_x2[chr_red_2] + dict_red_2
        doc_2.append(red_dict_y)# 得到红色y
        doc.append(red.ks_NG86)
    if iiss == 0:
        dot_red_1.extend(doc_1)
        dot_red_2.extend(doc_2)
        dot_c.extend(doc)
        if top_block == "True":# "Average"平均数"Median"中位数
            if ks_mod == "Average":
                if len(dot_block_ks_1) != 0:
                    average_block = sum(dot_block_ks_1) / len(dot_block_ks_1)  # 平均数
                    if ks_dem_1_1 <= average_block <= ks_dem_1:
                        dot_block_1.extend(doc_1)
                        dot_block_2.extend(doc_2)
            if ks_mod == "Median":
                if len(dot_block_ks_1) != 0:
                    median_block = statistics.median(dot_block_ks_1)  # 中位数
                    if ks_dem_1_1 <= median_block <= ks_dem_1:
                        dot_block_1.extend(doc_1)
                        dot_block_2.extend(doc_2)
        if len(doc_1) != 0:
            dot_star_1.append(doc_1[0])
            dot_end_1.append(doc_1[-1])
            dot_len_1.append(len(doc_1))
    if top_block == "True":
        return dot_red_1,dot_red_2,dot_c,dot_block_1,dot_block_2,dot_star_1,dot_end_1,dot_len_1
    if top_block == "False":
        return dot_red_1, dot_red_2, dot_c, dot_star_1, dot_end_1, dot_len_1
def tpos_ks_dot(final_result,gff_chr_x1,gff_chr_x2,gff_torder_x1,gff_torder_x2,gff_dict_x1,chr_dict_x1,gff_dict_x2,chr_dict_x2,top_block,ks_mod,ks_dem_1,ks_dem_1_1,name_id):
    dot_block_1,dot_block_2,dot_block_ks_1,dot_len_1,dot_star_1,dot_end_1,dot_red_1,dot_red_2,dot_c ,doc,doc_1,doc_2= [],[],[],[],[],[],[],[],[],[],[],[]
    iiss = 0
    for index, red in final_result.iterrows():
        if red.id1 == "#":
            if iiss == 0:
                dot_red_1.extend(doc_1)
                dot_red_2.extend(doc_2)
                dot_c.extend(doc)
                if top_block == "True":
                    # "Average"平均数"Median"中位数
                    if ks_mod == "Average":
                        if len(dot_block_ks_1) != 0:
                            average_block = sum(dot_block_ks_1) / len(dot_block_ks_1)  # 平均数
                            if ks_dem_1_1 <= average_block <= ks_dem_1:
                                dot_block_1.extend(doc_1)
                                dot_block_2.extend(doc_2)
                    if ks_mod == "Median":
                        if len(dot_block_ks_1) != 0:
                            median_block = statistics.median(dot_block_ks_1)  # 中位数
                            if ks_dem_1_1 <= median_block <= ks_dem_1:
                                dot_block_1.extend(doc_1)
                                dot_block_2.extend(doc_2)
                if len(doc_1) != 0:
                    dot_star_1.append(doc_1[0])
                    dot_end_1.append(doc_1[-1])
                    dot_len_1.append(len(doc_1))
            doc = []  # dot_c=[]
            doc_1 = []  # dot_red_1=[]
            doc_2 = []  # dot_red_2=[]
            dot_block_ks_1 = []
            iiss = 0
            continue
        if name_id == "True":
            if gff_chr_x1[red["id1"]] == gff_chr_x2[red["id2"]]:
                if (abs(gff_torder_x1[red["id1"]] - gff_torder_x2[red["id2"]]) < 150):
                    iiss = 1
        if red['ks_NG86'] == -1:
            continue
        if top_block == "True":
            dot_block_ks_1.append(red['ks_NG86'])
        chr_red_1 = gff_chr_x1[red["id1"]]
        dict_red_1 = gff_dict_x1[red["id1"]]
        red_dict_x = chr_dict_x1[chr_red_1] + dict_red_1
        doc_1.append(red_dict_x)
        # 得到红色x
        chr_red_2 = gff_chr_x2[red["id2"]]
        dict_red_2 = gff_dict_x2[red["id2"]]
        red_dict_y = chr_dict_x2[chr_red_2] + dict_red_2
        doc_2.append(red_dict_y)# 得到红色y
        doc.append(red.ks_NG86)
    if iiss == 0:
        dot_red_1.extend(doc_1)
        dot_red_2.extend(doc_2)
        dot_c.extend(doc)
        if top_block == "True":# "Average"平均数"Median"中位数
            if ks_mod == "Average":
                if len(dot_block_ks_1) != 0:
                    average_block = sum(dot_block_ks_1) / len(dot_block_ks_1)  # 平均数
                    if ks_dem_1_1 <= average_block <= ks_dem_1:
                        dot_block_1.extend(doc_1)
                        dot_block_2.extend(doc_2)
            if ks_mod == "Median":
                if len(dot_block_ks_1) != 0:
                    median_block = statistics.median(dot_block_ks_1)  # 中位数
                    if ks_dem_1_1 <= median_block <= ks_dem_1:
                        dot_block_1.extend(doc_1)
                        dot_block_2.extend(doc_2)
        if len(doc_1) != 0:
            dot_star_1.append(doc_1[0])
            dot_end_1.append(doc_1[-1])
            dot_len_1.append(len(doc_1))
    if top_block == "True":
        return dot_red_1,dot_red_2,dot_c,dot_block_1,dot_block_2,dot_star_1,dot_end_1,dot_len_1
    if top_block == "False":
        return dot_red_1, dot_red_2, dot_c, dot_star_1, dot_end_1, dot_len_1
def get_block(final_result,gff_chr_x1,gff_chr_x2,gff_torder_x1,gff_torder_x2,chr_torder_x1,gff_dict_x1,chr_dict_x1,chr_torder_x2,gff_dict_x2,chr_dict_x2,name_id,mod,Distance,N):
    alignment_id = []
    for index, Alignment in final_result.iterrows():  # 根据文件将不同区块分割 初步处理
        if Alignment["id1"] == "#":
            alignment = str(Alignment["id2"].split(":")[0])
            ns = int(Alignment["u"].split("=")[1])
            i = 1
            iis,ks_1,nan,nan_1 = 0,0,0,0
            ks = []
        else:
            if Alignment['ks_NG86'] == -1:
                if i == 1:
                    nan_1 += 1
                i += 1
                nan += 1# 当处理到最后基因对为-1时触发区块保存
                if i > ns:
                    if (ns - nan) == 0:
                        continue
                    avg_ks = ks_1 / (ns - nan)  # 计算平均KS值（防止除零）
                    median_ks = statistics.median(ks)
                    # 计算中点坐标（使用第1个和最后有效的位置）
                    mid_x_torder = (red_dict_x + red_dict_x_2) / 2
                    mid_y_torder = (red_dict_y + red_dict_y_2) / 2
                    mid_x_tpos = (red_dict_x_1 + red_dict_x_2_1) / 2
                    mid_y_tpos = (red_dict_y_1 + red_dict_y_2_1) / 2
                    align = AlignmentResult(            # 创建结果对象
                        alignment,  # 比对名称
                        (ns - nan),  # 有效基因对数
                        mid_x_torder, mid_y_torder,  # 中点坐标
                        mid_x_tpos, mid_y_tpos,
                        avg_ks, median_ks,  # KS统计
                        red_dict_x, red_dict_x_2,  # X起点和终点（第1和第5基因）
                        red_dict_y, red_dict_y_2,  # Y起点和终点
                        red_dict_x_1, red_dict_x_2_1,  # X的tpos起点和终点
                        red_dict_y_1, red_dict_y_2_1,  # Y的tpos起点和终点
                        ks  # KS值列表（有效值）
                    )
                    if iis == 0:
                        alignment_id.append(align)  # print(f"{align.alignment} {align.ns}  {align.ks1}")
                continue  # 跳过无效基因对的后续处理
            if name_id == "True":
                if gff_chr_x1[Alignment["id1"]] == gff_chr_x2[Alignment["id2"]]:
                    if (abs(gff_torder_x1[Alignment["id1"]] - gff_torder_x2[Alignment["id2"]]) < 150):
                        iis = 1
            ks_1 += float(Alignment["ks_NG86"])
            ks.append(float(Alignment["ks_NG86"]))
            if i == 1 + nan_1:
                chr_red_1 = gff_chr_x1[Alignment["id1"]]
                dict_red_1 = gff_torder_x1[Alignment["id1"]]
                red_dict_x = chr_torder_x1[chr_red_1] + dict_red_1  # "torder"  x1
                chr_red_1_1 = gff_chr_x1[Alignment["id1"]]
                dict_red_1_1 = gff_dict_x1[Alignment["id1"]]
                red_dict_x_1 = chr_dict_x1[chr_red_1_1] + dict_red_1_1  # "tpos"   x1
                chr_red_2 = gff_chr_x2[Alignment["id2"]]
                dict_red_2 = gff_torder_x2[Alignment["id2"]]
                red_dict_y = chr_torder_x2[chr_red_2] + dict_red_2  # "torder"  y1
                chr_red_2_1 = gff_chr_x2[Alignment["id2"]]
                dict_red_2_1 = gff_dict_x2[Alignment["id2"]]
                red_dict_y_1 = chr_dict_x2[chr_red_2_1] + dict_red_2_1  # "tpos"   y1
            chr_red_1_2 = gff_chr_x1[Alignment["id1"]]
            dict_red_1_2 = gff_torder_x1[Alignment["id1"]]
            red_dict_x_2 = chr_torder_x1[chr_red_1_2] + dict_red_1_2  # X终点 torder
            red_dict_y_2 = chr_torder_x2[gff_chr_x2[Alignment["id2"]]] + gff_torder_x2[Alignment["id2"]]  # Y终点 torder
            red_dict_x_2_1 = chr_dict_x1[chr_red_1_2] + gff_dict_x1[Alignment["id1"]]  # X终点 tpos
            red_dict_y_2_1 = chr_dict_x2[gff_chr_x2[Alignment["id2"]]] + gff_dict_x2[Alignment["id2"]]
            if i == ns:         # 所有有效基因对都会更新 处理了最后一个基因对ks值为-1的情况，确保区块能被正确处理
                chr_red_1_2 = gff_chr_x1[Alignment["id1"]]
                dict_red_1_2 = gff_torder_x1[Alignment["id1"]]
                red_dict_x_2 = chr_torder_x1[chr_red_1_2] + dict_red_1_2  # "torder"  x2
                chr_red_1_1_1 = gff_chr_x1[Alignment["id1"]]
                dict_red_1_1_1 = gff_dict_x1[Alignment["id1"]]
                red_dict_x_2_1 = chr_dict_x1[chr_red_1_1_1] + dict_red_1_1_1  # "tpos"   x2
                chr_red_2_2 = gff_chr_x2[Alignment["id2"]]
                dict_red_2_2 = gff_torder_x2[Alignment["id2"]]
                red_dict_y_2 = chr_torder_x2[chr_red_2_2] + dict_red_2_2  # "torder"  y2
                chr_red_2_1_1 = gff_chr_x2[Alignment["id2"]]
                dict_red_2_1_1 = gff_dict_x2[Alignment["id2"]]
                red_dict_y_2_1 = chr_dict_x2[chr_red_2_1_1] + dict_red_2_1_1  # "tpos"   y2
                red_dict_x_3 = (red_dict_x + red_dict_x_2) / 2  # "torder"  x
                red_dict_y_3 = (red_dict_y + red_dict_y_2) / 2  # "torder"  y
                red_dict_x_4 = (red_dict_x_1 + red_dict_x_2_1) / 2  # "tpos"     x
                red_dict_y_4 = (red_dict_y_1 + red_dict_y_2_1) / 2  # "tpos"     y
                if ns - nan == 0:
                    continue  # 当区块全为-1时跳过区块
                ks_1 = ks_1 / (ns - nan)  # ks的平均数ks
                ks_2 = statistics.median(ks)  # ks的中位数
                align = AlignmentResult(alignment, (ns - nan),
                                        red_dict_x_3,red_dict_y_3,red_dict_x_4,
                                        red_dict_y_4,ks_1, ks_2,red_dict_x,
                                        red_dict_x_2,red_dict_y,red_dict_y_2,
                                        red_dict_x_1,red_dict_x_2_1,red_dict_y_1,
                                        red_dict_y_2_1,ks
                                        )
                if iis == 0:
                    alignment_id.append(align)
            i += 1
    i = 0
    if mod == "torder":
        while i < len(alignment_id) - 1:
            current_line = alignment_id[i]# print(f"{current_line.alignment}    {current_line.ns}  {current_line.ks1}")
            if current_line.ns < N:
                i += 1
                continue
            next_line = alignment_id[i + 1]# print(f"{current_line.alignment}    {current_line.ns}  {current_line.ks1}")
            dx = (next_line.x_start1 - current_line.x_end1)
            dy = (next_line.y_start1 - current_line.y_end1)
            distance = math.sqrt(dx ** 2 + dy ** 2) # 计算当前线终点和下一条线起点的距离
            if distance < Distance: # 合并两条线（取第一条线的起点和第二条线的终点）
                merged_line = AlignmentResult(current_line.alignment, (current_line.ns + next_line.ns),
                                              ((current_line.x_start1 + next_line.x_end1) / 2),
                                              ((current_line.y_start1 + next_line.y_end1) / 2),
                                              ((current_line.x_start2 + next_line.x_end2) / 2),
                                              ((current_line.y_start2 + next_line.y_end2) / 2),
                                              (((current_line.ks1 * current_line.ns) + (
                                                      next_line.ks1 * next_line.ns)) / (
                                                       current_line.ns + next_line.ns)),
                                              (statistics.median(current_line.ks + next_line.ks)),
                                              current_line.x_start1, next_line.x_end1,
                                              current_line.y_start1, next_line.y_end1,
                                              current_line.x_start2, next_line.x_end2,
                                              current_line.y_start2, next_line.y_end2,
                                              (current_line.ks + next_line.ks)
                                              )
                alignment_id[i] = merged_line
                alignment_id.pop(i + 1)  # 跳过下一条线，因为已经合并
            else:
                i += 1  # 不合并，保留当前线
    elif mod == "tpos":
        while i < len(alignment_id) - 1:
            current_line = alignment_id[i]
            if current_line.ns < N:
                i += 1
                continue
            # 如果是最后一条线，直接加入
            next_line = alignment_id[i + 1]
            # 计算当前线终点和下一条线起点的距离
            dx = (next_line.x_start2 - current_line.x_end2)
            dy = (next_line.y_start2 - current_line.y_end2)
            distance = math.sqrt(dx ** 2 + dy ** 2)
            if distance < Distance:      # 合并两条线（取第一条线的起点和第二条线的终点）
                merged_line = AlignmentResult(current_line.alignment, (current_line.ns + next_line.ns),
                                              ((current_line.x_start1 + next_line.x_end1) / 2),
                                              ((current_line.y_start1 + next_line.y_end1) / 2),
                                              ((current_line.x_start2 + next_line.x_end2) / 2),
                                              ((current_line.y_start2 + next_line.y_end2) / 2),
                                              (((current_line.ks1 * current_line.ns) + (
                                                      next_line.ks1 * next_line.ns)) / (
                                                       current_line.ns + next_line.ns)),
                                              (statistics.median(current_line.ks + next_line.ks)),
                                              current_line.x_start1, next_line.x_end1,
                                              current_line.y_start1, next_line.y_end1,
                                              current_line.x_start2, next_line.x_end2,
                                              current_line.y_start2, next_line.y_end2,
                                              (current_line.ks + next_line.ks)
                                              )
                alignment_id[i] = merged_line
                alignment_id.pop(i + 1)  # 跳过下一条线，因为已经合并
            else:# 不合并，保留当前线
                i += 1  # 不合并，保留当前线
    return alignment_id
def rate_date(blast_ld,gff_chr_x1,gff_chr_x2,gff_torder_x1,gff_torder_x2,name_id):
    alignment_new = []
    alignment_new_1 = []
    iss = 0
    for index, Alignment in blast_ld.iterrows():
        if Alignment["id1"] == "#":
            if iss == 0:
                alignment_new = alignment_new + alignment_new_1
            iss = 0
            alignment_new_1 = []
            continue
        else:
            if name_id == "True":       # 是否去除自身与自身比对
                if gff_chr_x1[Alignment["id1"]] == gff_chr_x2[Alignment["id2"]]:
                    if (abs(gff_torder_x1[Alignment["id1"]] - gff_torder_x2[Alignment["id2"]]) < 150):
                        iss = 1
            alignment_new_1.append(Alignment)
    if iss == 0:
        alignment_new = alignment_new + alignment_new_1
    alignment_df = pd.DataFrame(alignment_new)
    alignment_id1 = alignment_df.groupby('id1')['id2'].size().reset_index(
        name='count')  # .groupby() 按 () 列的值对 DataFrame 进行分组 [] 从分组后的数据中选择 [] 列
    alignment_id2 = alignment_df.groupby('id2')['id1'].size().reset_index(
        name='count')  # .size()计算每个分组中有多少个 [] 列值 返回一个 Pandas Series，其中索引是 id1 的值，值是对应的计数
    # reset_index() 将 id1 从索引变回普通列 name='count' 指定计数列的列名为 "count"  得到1个 DataFrame 第1列为基因/元素 第2列为与多少个其他基因/元素有关联
    alignment_id1_date = {}
    alignment_id2_date = {}
    id1_date = 0
    id2_date = 0
    for i in alignment_id1["count"]:  # 计算基因对比第1列比第2列
        if i > 5:
            continue
        if i not in alignment_id1_date.keys():
            alignment_id1_date[i] = 1
        if i in alignment_id1_date.keys():
            alignment_id1_date[i] = alignment_id1_date[i] + 1
        id1_date += 1
    for i in alignment_id2["count"]:  # 计算基因对比第2列比第1列
        if i > 5:
            continue
        if i not in alignment_id2_date.keys():
            alignment_id2_date[i] = 1
        if i in alignment_id2_date.keys():
            alignment_id2_date[i] = alignment_id2_date[i] + 1
        id2_date += 1
    for key, value in alignment_id1_date.items():
        alignment_id1_date[key] = int(value / id1_date * 100 + 0.5)
    for key, value in alignment_id2_date.items():
        alignment_id2_date[key] = int(value / id2_date * 100 + 0.5)
    alignment_id1_date == dict(sorted(alignment_id1_date.items()))  # 默认按 key 升序    得到第1列比第2列的比例
    alignment_id2_date == dict(sorted(alignment_id2_date.items()))  # 默认按 key 升序    得到第2列比第1列的比例
    return alignment_id1_date,alignment_id2_date
def get_density(gff_1,long1,gff_2,long1_1,chr_dict_1,gff_dict_1,chr_dict_2,gff_dict_2,bin_sizes,n,m,n_2,m_2,ax_top,ax_right,name_1,name_2):
    gff_1_new = []
    gff_2_new = []
    for gff in gff_1.itertuples(index=False):
        if gff.chr not in long1.chr.values:
            continue
        else:
            gff_1_new.append(gff)
    new_gff_1 = pd.DataFrame(gff_1_new,
                             columns=["chr", "tnew_name", "tpos_start", "tpos_end", "tchain", "torder", "told_name","tops"]
                             )
    for gff in gff_2.itertuples(index=False):
        if gff.chr not in long1_1.chr.values:
            continue
        else:
            gff_2_new.append(gff)
    new_gff_2 = pd.DataFrame(gff_2_new,
                             columns=["chr", "tnew_name", "tpos_start", "tpos_end", "tchain", "torder", "told_name", "tops"]
                             )
    new_gff_1["ops"] = new_gff_1["chr"].map(chr_dict_1) + new_gff_1["tnew_name"].map(gff_dict_1)
    new_gff_2["ops"] = new_gff_2["chr"].map(chr_dict_2) + new_gff_2["tnew_name"].map(gff_dict_2)
    new_gff_1_i = []
    new_gff_2_i = []
    new_gff_1_le = []
    new_gff_2_le = []
    ture_gff_1 = []
    bin_size = bin_sizes * 1000000  # 定义分箱大小
    i = 0
    while i < n:
        ture_gff_1.append(i / 1000000)
        i += 50 * 1000000
    i = 0
    while i < n:
        k = i
        i += bin_size
        current_bin_end = min(i, n)  # 确保不超过最大值m
        bin_center = (k + current_bin_end) / 2 / 1000000  # 转换为百万单位 # 计算分箱中点
        count_in_bin = ((new_gff_1["ops"] > k) & (new_gff_1["ops"] <= current_bin_end)).sum()# 统计当前分箱内的数据点数量
        # 计算百分比
        if current_bin_end < n:
            percentage = (count_in_bin * 100) / n_2
        else:# 处理最后一个不完整的分箱
            actual_bin_size = n - k
            scaling_factor = actual_bin_size / bin_size
            percentage = (count_in_bin * 100) / (n_2 * scaling_factor)
        new_gff_1_i.append(bin_center)
        new_gff_1_le.append(percentage)
    ture_gff_2 = []
    bin_size = bin_sizes * 1000000  # 定义分箱大小
    i = 0
    while i < m:
        ture_gff_2.append(i / 1000000)
        i += 50 * 1000000
    i = 0
    while i < m:
        k = i
        i += bin_size
        current_bin_end = min(i, m)  # 确保不超过最大值m
        # 计算分箱中点
        bin_center = (k + current_bin_end) / 2 / 1000000  # 转换为百万单位
        # 统计当前分箱内的数据点数量
        count_in_bin = ((new_gff_2["ops"] > k) & (new_gff_2["ops"] <= current_bin_end)).sum()
        # 计算百分比
        if current_bin_end < m:
            percentage = (count_in_bin * 100) / m_2
        else:        # 处理最后一个不完整的分箱
            actual_bin_size = m - k
            scaling_factor = actual_bin_size / bin_size
            percentage = (count_in_bin * 100) / (m_2 * scaling_factor)
        new_gff_2_i.append(bin_center)
        new_gff_2_le.append(percentage)
    max_gff_1 = ((max(new_gff_1_le) + 9) // 10) * 10
    max_gff_2 = ((max(new_gff_2_le) + 9) // 10) * 10
    ax_top.plot(new_gff_1_i, new_gff_1_le, marker='.', linestyle='-', color="#FFCD16")
    ax_top.set_xlim(0, max(new_gff_1_i) * 1.03)
    ax_top.set_ylim(0, max_gff_1)
    ax_top.spines['top'].set_visible(False)    # 隐藏上边框和右边框
    ax_top.spines['right'].set_visible(False)    # 隐藏上边框和右边框
    ax_top.set_xticks([])  # 隐藏x轴刻度     # ax_top.set_yticks([])  # 隐藏y轴刻度
    for i in range(len(ture_gff_1)):
        ax_top.text(ture_gff_1[i], -0.05 * max_gff_1, f"{int(ture_gff_1[i])}", ha='center', va='center',
                    fontsize=8, weight='bold')
    for i in range(len(new_gff_1_i)):
        ax_top.text(new_gff_1_i[i], new_gff_1_le[i], f"{new_gff_1_le[i]:.1f}", ha='center', va='bottom', fontsize=9,
                    weight='bold')
    ax_top.text(max(new_gff_1_i) * 1.03, -0.05 * max_gff_1, "M", ha='center', va='center', fontsize=9,
                weight='bold')
    ax_top.text(0, 1.06 * max_gff_1, "%", ha='center', va='center', fontsize=9,
                weight='bold')
    ax_top.text(max(new_gff_1_i) * 0.5, 1.05 * max_gff_1, fr"$\boldsymbol{{{name_1}}}$ density", ha='center',
                va='center',
                fontsize=9, weight='bold')
    ax_right.plot(new_gff_2_le, new_gff_2_i, marker='.', linestyle='-', color="#FFCD16")
    ax_right.set(xlim=(0, max_gff_2),  # x轴范围保持0-100
                 ylim=(0, max(new_gff_2_i) * 1.03))  # y轴范围调整使柱子居中
    ax_right.xaxis.set_ticks_position('top')
    ax_right.tick_params(axis='x', labelrotation=270)
    ax_right.tick_params(axis='y', labelrotation=270)
    ax_right.invert_yaxis()  # 保持y轴反转
    ax_right.spines['bottom'].set_visible(False)    # 边框控制
    ax_right.spines['right'].set_visible(False)    # 边框控制
    # ax_right.set_xticks([])  # 隐藏x轴刻度
    ax_right.set_yticks([])  # 隐藏y轴刻度
    for i in range(len(ture_gff_2)):
        ax_right.text(-0.05 * max_gff_2, ture_gff_2[i], f"{int(ture_gff_2[i])}", ha='center', va='center',
                      fontsize=8, weight='bold', rotation=270)
    for i in range(len(new_gff_2_i)):
        ax_right.text(new_gff_2_le[i] * 1.03, new_gff_2_i[i], f"{new_gff_2_le[i]:.1f}", ha='center', va='center',
                      fontsize=8, weight='bold', rotation=270)
    ax_right.text(-0.05 * max_gff_2, max(new_gff_2_i) * 1.03, "M", ha='center', va='center', fontsize=9,
                  weight='bold', rotation=270)
    ax_right.text(1.06 * max_gff_2, 0, "%", ha='center', va='center', fontsize=9,
                  weight='bold', rotation=270)
    ax_right.text(1.05 * max_gff_2, max(new_gff_2_i) * 0.5, fr"$\boldsymbol{{{name_2}}}$ density", ha='center',
                  va='center',
                  fontsize=9, weight='bold', rotation=270)
def half_col_2_1(ax,ax_4,posd):
    if posd == "bot_left":
        ax.invert_yaxis()  # 反转y轴，y轴的方向将变为递减向上
        ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
        ax_4.invert_xaxis()  # 反转x轴，x轴的方向将变为递减向上
        ax_4.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    if posd == "bot_right":
        ax_4.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
        ax.invert_yaxis()  # 反转y轴，y轴的方向将变为递减向上
        ax.invert_xaxis()  # 反转x轴，x轴的方向将变为递减向上
        ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    ax.spines['top'].set_visible(False)  # 隐藏上边框
    ax.spines['right'].set_visible(False)  # 隐藏右边框
    ax.spines['left'].set_visible(False)  # 隐藏左边框
    ax.spines['bottom'].set_visible(False)  # 隐藏下边框
    ax_4.spines['top'].set_visible(False)  # 隐藏上边框
    ax_4.spines['right'].set_visible(False)  # 隐藏右边框
    ax_4.spines['left'].set_visible(False)  # 隐藏左边框
    ax_4.spines['bottom'].set_visible(False)  # 隐藏下边框
def half_col_2_2(ax,long1,long1_1,x_lens,y_lens,m,n,x_len,y_len,dot_len_1,dot_star_1,dot_end_1,name_1,name_2,block):
    # 使用plt.text()添加染色体标签
    if block == "False":
        for i, chr_num in enumerate(long1["chr"]):
            ax.text(x_lens[i], 1.02 * m, f" {str(chr_num)}", ha='center', va='center', fontsize=9,
                    weight='bold')  # ha='center'	文本的中心对齐坐标点 va='top'文本的顶部对齐坐标点
    for i, chr_num in enumerate(long1_1["chr"]):
        ax.text(-0.02 * n, y_lens[i], f" {str(chr_num)}", ha='center', va='center', fontsize=9,
                weight='bold')  # ha='right' 文本的右边界对齐坐标点  va='center'文本的中线对齐坐标点
    # plt.text(x坐标位置, y坐标位置, 染色体编号,水平对齐方式, 垂直对齐方式, 字体大小。)
    for x in x_len:
        ax.vlines(x, x, m, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    # ax.vlines(x_len, 0, m-x_len, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    for y in y_len:
        ax.hlines(y, 0, y, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画水平线
    # ax.hlines(y_len, 0,y_len, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画水平线
    if block == "True":
        for i in range(len(dot_len_1)):
            if dot_len_1[i] < 10:
                rect = Rectangle((dot_star_1[i], m * 1.01), abs(dot_end_1[i] - dot_star_1[i]), m * 0.02,
                                 color='blue', alpha=1, linewidth=0, clip_on=False)
                ax.add_patch(rect)
        for i in range(len(dot_len_1)):
            if 10 < dot_len_1[i] < 20:
                rect = Rectangle((dot_star_1[i], m * 1.01), abs(dot_end_1[i] - dot_star_1[i]), m * 0.02,
                                 color='#ADD8E6', alpha=1, linewidth=0, clip_on=False)
                ax.add_patch(rect)
        for i in range(len(dot_len_1)):
            if 20 < dot_len_1[i] < 30:
                rect = Rectangle((dot_star_1[i], m * 1.01), abs(dot_end_1[i] - dot_star_1[i]), m * 0.02,
                                 color='yellow', alpha=1, linewidth=0, clip_on=False)
                ax.add_patch(rect)
        for i in range(len(dot_len_1)):
            if 30 < dot_len_1[i] < 40:
                rect = Rectangle((dot_star_1[i], m * 1.01), abs(dot_end_1[i] - dot_star_1[i]), m * 0.02,
                                 color='orange', alpha=1, linewidth=0, clip_on=False)
                ax.add_patch(rect)
        for i in range(len(dot_len_1)):
            if 40 < dot_len_1[i]:
                rect = Rectangle((dot_star_1[i], m * 1.01), abs(dot_end_1[i] - dot_star_1[i]), m * 0.02,
                                 color='r', linewidth=0, alpha=1, clip_on=False)
                ax.add_patch(rect)
    ax.vlines(0, 0, m, lw=1.5, color='k', linestyle='-', alpha=1)  # 画外框架
    ax.hlines(m, 0, n, lw=1.5, color='k', linestyle='-', alpha=1)  # 画外框架
    ax.plot([0, m], [0, n], lw=1.5, color='k', linestyle='-', alpha=0.5)
    ax.text(n / 2, 1.1 * m, fr"$\boldsymbol{{{name_1}}}$", fontsize=14, ha='center', va='top')  # 添加坐标轴标签
    ax.text(-0.1 * n, m / 2, fr"$\boldsymbol{{{name_2}}}$", fontsize=14, ha='right', va='center',rotation=90)  # 添加坐标轴标签
    # ax.text(x坐标位置, y坐标位置, 使用LaTeX语法加粗显示物种名称,字体大小，水平对齐方式, 垂直对齐方式,旋转角度)
def half_dot_1(dot_red_1,dot_red_2):
    reds_2 = []
    reds_1 = []
    for i in range(len(dot_red_1)):
        if dot_red_2[i] >= dot_red_1[i]:
            reds_1.append(dot_red_1[i])
            reds_2.append(dot_red_2[i])
    return reds_1,reds_2
def half_f_1(ax):
    ax.invert_yaxis()  # 反转y轴，y轴的方向将变为递减向上
    ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    ax.spines['top'].set_visible(False)  # 隐藏上边框
    ax.spines['right'].set_visible(False)  # 隐藏右边框
    ax.spines['left'].set_visible(False)  # 隐藏左边框
    ax.spines['bottom'].set_visible(False)  # 隐藏下边框
    # 使用plt.text()添加染色体标签
def half_f_2(ax,long1,long1_1,x_lens,y_lens,x_len,y_len,m,n,name_1,name_2):
    for i, chr_num in enumerate(long1["chr"]):
        ax.text(x_lens[i], 1.02 * m, f" {str(chr_num)}", ha='center', va='center', fontsize=9,weight='bold')  # ha='center'	文本的中心对齐坐标点 va='top'文本的顶部对齐坐标点
    for i, chr_num in enumerate(long1_1["chr"]):
        ax.text(-0.02 * n, y_lens[i], f" {str(chr_num)}", ha='center', va='center', fontsize=9,weight='bold')  # ha='right' 文本的右边界对齐坐标点  va='center'文本的中线对齐坐标点
    # plt.text(x坐标位置, y坐标位置, 染色体编号,水平对齐方式, 垂直对齐方式, 字体大小。)
    ax.vlines(x_len, 0, m, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画水平线
    ax.hlines(y_len, 0, n, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    ax.vlines((0, x_len[-1]), 0, m, lw=1.5, color='k', linestyle='-', alpha=1)  # 画外框架
    ax.hlines((0, y_len[-1]), 0, n, lw=1.5, color='k', linestyle='-', alpha=1)  # 画外框架
    ax.text(n / 2, 1.1 * m, fr"$\boldsymbol{{{name_1}}}$", fontsize=12, ha='center', va='top')  # 添加坐标轴标签
    ax.text(-0.1 * n, m / 2, fr"$\boldsymbol{{{name_2}}}$", fontsize=12, ha='right', va='center',rotation=90)  # 添加坐标轴标签
    # ax.text(x坐标位置, y坐标位置, 使用LaTeX语法加粗显示物种名称,字体大小，水平对齐方式, 垂直对齐方式,旋转角度)
def half_t_1(pos,ax):
    if pos == "bot_left":
        ax.invert_yaxis()  # 反转y轴，y轴的方向将变为递减向上
        ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    if pos == "top_left":
        ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    if pos == "bot_right":
        ax.invert_yaxis()  # 反转y轴，y轴的方向将变为递减向上
        ax.invert_xaxis()  # 反转x轴，x轴的方向将变为递减向上
        ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    if pos == "top_right":
        ax.invert_xaxis()  # 反转x轴，x轴的方向将变为递减向上
        ax.xaxis.set_ticks_position('top')  # 将x轴的刻度线移动到坐标轴的上方
    ax.spines['top'].set_visible(False)  # 隐藏上边框
    ax.spines['right'].set_visible(False)  # 隐藏右边框
    ax.spines['left'].set_visible(False)  # 隐藏左边框
    ax.spines['bottom'].set_visible(False)  # 隐藏下边框
    # 使用plt.text()添加染色体标签
def half_t_2(ax,long1,long1_1,x_lens,y_lens,x_len,y_len,m,n,name_1,name_2):
    for i, chr_num in enumerate(long1["chr"]):
        ax.text(x_lens[i], 1.02 * m, f" {str(chr_num)}", ha='center', va='center', fontsize=9, weight='bold')  # ha='center'	文本的中心对齐坐标点 va='top'文本的顶部对齐坐标点
    for i, chr_num in enumerate(long1_1["chr"]):
        ax.text(-0.02 * n, y_lens[i], f" {str(chr_num)}", ha='center', va='center', fontsize=9, weight='bold')  # ha='right' 文本的右边界对齐坐标点  va='center'文本的中线对齐坐标点
    # plt.text(x坐标位置, y坐标位置, 染色体编号,水平对齐方式, 垂直对齐方式, 字体大小。)
    for x in x_len:
        ax.vlines(x, x, m, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    # ax.vlines(x_len, 0, m-x_len, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    for y in y_len:
        ax.hlines(y, 0, y, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画水平线
    # ax.hlines(y_len, 0,y_len, lw=0.5, color='k', linestyle='-', alpha=0.5)  # 画水平线
    ax.vlines(0, 0, m, lw=1.5, color='k', linestyle='-', alpha=1)  # 画外框架
    ax.hlines(m, 0, n, lw=1.5, color='k', linestyle='-', alpha=1)  # 画外框架
    ax.plot([0, m], [0, n], lw=1.5, color='k', linestyle='-', alpha=0.5)
    ax.text(n / 2, 1.1 * m, fr"$\boldsymbol{{{name_1}}}$", fontsize=12, ha='center', va='top')  # 添加坐标轴标签
    ax.text(-0.1 * n, m / 2, fr"$\boldsymbol{{{name_2}}}$", fontsize=12, ha='right', va='center',rotation=90)  # 添加坐标轴标签
    # ax.text(x坐标位置, y坐标位置, 使用LaTeX语法加粗显示物种名称,字体大小，水平对齐方式, 垂直对齐方式,旋转角度)
def get_tree_1(gene_show,gff_1,gff_chr_1,gff_torder_1,chr_torder_1,gff_dict_1,chr_dict_1,gff_2,gff_chr_2,gff_torder_2,chr_torder_2,gff_dict_2,chr_dict_2,mod):
    tree_1 = pd.read_csv(
        gene_show,
        sep='\t',
        header=None,
        usecols=[0, 1],
        names=["chr", "tnew_name"],
        dtype={
            "chr": int,
            "tnew_name": str,
            "tpos_start": int,
            "tpos_end": int,
            "tchain": str,
            "torder": int,
            "told_name": str
        }
    )
    tree_1_sorted = tree_1.sort_values(by=["chr", "torder"])  # 进行排序  先按染色体，再按顺序
    x_gene = []
    y_gene = []
    for i in tree_1_sorted["tnew_name"]:
        date_xy = []
        x_gene_1 = []
        y_gene_1 = []
        x_gene_1s = []
        y_gene_1s = []
        if i in gff_1["tnew_name"].values:  # y=max
            if mod == "torder":
                tree_chr_1 = gff_chr_1[i]
                tree_dict_1 = gff_torder_1[i]
                tree_dict_x = chr_torder_1[tree_chr_1] + tree_dict_1
            elif mod == "tpos":
                tree_chr_1 = gff_chr_1[i]
                tree_dict_1 = gff_dict_1[i]
                tree_dict_x = chr_dict_1[tree_chr_1] + tree_dict_1
            x_gene_1.append(i)
            x_gene_1s.append(tree_dict_x)
            x_genes = gene_id(i,
                              tree_dict_x,
                              x_gene_1,
                              x_gene_1s,
                              date_xy
                              )
            x_gene.append(x_genes)
        elif i in gff_2["tnew_name"].values:  # x=max
            if mod == "torder":
                tree_chr_2 = gff_chr_2[i]
                tree_dict_2 = gff_torder_2[i]
                tree_dict_y = chr_torder_2[tree_chr_2] + tree_dict_2
            elif mod == "tpos":
                tree_chr_2 = gff_chr_2[i]
                tree_dict_2 = gff_dict_2[i]
                tree_dict_y = chr_dict_2[tree_chr_2] + tree_dict_2
            y_gene_1.append(i)
            y_gene_1s.append(tree_dict_y)
            y_genes = gene_id(i,
                              tree_dict_y,
                              y_gene_1,
                              y_gene_1s,
                              date_xy
                              )
            y_gene.append(y_genes)
        else:
            print(f"{i}不存在于GFF文件中")
    return x_gene,y_gene
def get_tree_2(x_gene):
    i = 0
    while i < len(x_gene):
        date_xy = []
        current_x = x_gene[i]
        if i < len(x_gene) - 1:
            next_x = x_gene[i + 1]
            if abs(next_x.xy - current_x.xys[-1]) < 11:
                merged_genes = current_x.genes + [next_x.gene_name]
                merged_xys = current_x.xys + [next_x.xy]
                merged_x = gene_id(
                    current_x.gene_name,
                    (current_x.xys[0] + next_x.xy) / 2,
                    merged_genes,
                    merged_xys,
                    date_xy
                )
                x_gene[i] = merged_x
                kx = x_gene.pop(i + 1)
            else:
                date_xy = np.linspace(current_x.xy - (len(current_x.xys) - 1) * wide / 2,
                                      current_x.xy + (len(current_x.xys) - 1) * wide / 2,
                                      len(current_x.xys))
                merged_x = gene_id(
                    current_x.gene_name,
                    current_x.xy,
                    current_x.genes,
                    current_x.xys,
                    date_xy
                )
                x_gene[i] = merged_x
                i = i + 1
        else:
            date_xy = np.linspace(current_x.xy - (len(current_x.xys) - 1) * wide / 2,
                                  current_x.xy + (len(current_x.xys) - 1) * wide / 2,
                                  len(current_x.xys))
            merged_x = gene_id(
                current_x.gene_name,
                current_x.xy,
                current_x.genes,
                current_x.xys,
                date_xy
            )
            x_gene[i] = merged_x
            i = i + 1
    return x_gene
def tree_co_1(color_date,color_dates):
    if color_date == "False":
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        date_colors='NAN'
    if color_date == "True":
        date_color = pd.read_csv(
            color_dates,
            sep='\t',
            header=None,
            names=["tnew_name", "color"],
            dtype={
                "tnew_name": str,
                "color": str,
            }
        )
        date_colors = date_color.set_index("tnew_name")["color"].to_dict()
        colors='NAN'
    return colors,date_colors
def tree_y_gene_1(y_gene,mod,n_2,color_date,colors,ax,date_colors,n):
    for idx, pos in enumerate(y_gene):  # enumerate() 会返回一个迭代器，生成 (index, value) 元组，其中：idx 是当前元素的索引（从 0 开始）。pos 是 y_gene 中的当前元素
        if mod == "torder":
            has_overlap = False
            # 检查当前元素的 date_xy 是否与其他元素的 date_xy 相交
            for j, pos2 in enumerate(y_gene):
                if idx == j:
                    continue
                start1, end1 = min(pos.date_xy), max(pos.date_xy)
                start2, end2 = min(pos2.date_xy), max(pos2.date_xy)

                if not (end1 < start2 or end2 < start1):
                    has_overlap = True
                    break  # 只要有一个重叠，就可以停止检查

            # 根据是否有交集来统一计算 name_xs
            if has_overlap:
                name_xs = [1.07 * n_2 + i * longs for i in range(len(pos.genes))]
            else:
                name_xs = [1.07 * n_2 for _ in range(len(pos.genes))]
            if color_date == "False":
                line_color = colors[
                    idx % len(colors)]  # 取模运算，确保索引 idx 始终落在 colors 的有效范围内  当 idx >= len(colors) 时，循环回到列表开头
                if len(pos.xys) == 1:
                    # 单个名称直接连接
                    ax.plot([1.05 * n_2, name_xs[0]], [pos.xy, pos.xy], color=line_color, lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        # 从垂直线到名称的水平线
                        ax.plot([n_2, 1.05 * n_2], [pos.xys[i], pos.date_xy[i]], color=line_color, lw=1,
                                clip_on=False)
                        ax.plot([1.05 * n_2, name_xs[i]], [pos.date_xy[i], pos.date_xy[i]], color=line_color, lw=1,
                                clip_on=False)
                    ax.text(name_xs[i], pos.date_xy[i], f'{pos.genes[i]}',
                            va='center', ha='left', fontsize=8)
            if color_date == "True":
                if len(pos.xys) == 1:
                    # 单个名称直接连接
                    if pos.gene_name in date_colors.keys():
                        ax.plot([1.05 * n_2, name_xs[0]], [pos.xy, pos.xy], color=date_colors[pos.gene_name], lw=1,
                                clip_on=False)
                    else:
                        ax.plot([1.05 * n_2, name_xs[0]], [pos.xy, pos.xy], color="#808080", lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        # 从垂直线到名称的水平线
                        if pos.genes[i] in date_colors.keys():
                            ax.plot([n_2, 1.05 * n_2], [pos.xys[i], pos.date_xy[i]],
                                    color=date_colors[pos.genes[i]],
                                    lw=1, clip_on=False)
                            ax.plot([1.05 * n_2, name_xs[i]], [pos.date_xy[i], pos.date_xy[i]],
                                    color=date_colors[pos.genes[i]], lw=1, clip_on=False)
                        else:

                            ax.plot([n_2, 1.05 * n_2], [pos.xys[i], pos.date_xy[i]], color="#808080",
                                    lw=1,
                                    clip_on=False)
                            ax.plot([1.05 * n_2, name_xs[i]], [pos.date_xy[i], pos.date_xy[i]],
                                    color="#808080", lw=1,
                                    clip_on=False)
                    ax.text(name_xs[i], pos.date_xy[i], f'{pos.genes[i]}',
                            va='center', ha='left', fontsize=8)
        if mod == "tpos":
            has_overlap = False
            for j, pos2 in enumerate(y_gene):
                if idx == j:
                    continue
                start1, end1 = min(pos.date_xy), max(pos.date_xy)
                start2, end2 = min(pos2.date_xy), max(pos2.date_xy)

                # 判断是否有重叠
                if not (end1 < start2 or end2 < start1):
                    # 计算名称坐标（自动防重叠）
                    has_overlap = True
                    break  # 只要有一个重叠，就可以停止检查
            # 计算名称坐标（自动防重叠）
            if has_overlap:
                name_xs = [1.07 * n + i * longs for i in range(len(pos.genes))]
            else:
                name_xs = [1.07 * n for _ in range(len(pos.genes))]
            if color_date == "False":
                line_color = colors[
                    idx % len(colors)]  # 取模运算，确保索引 idx 始终落在 colors 的有效范围内  当 idx >= len(colors) 时，循环回到列表开头
                if len(pos.xys) == 1:
                    ax.plot([1.05 * n, name_xs[0]], [pos.xy, pos.xy], color=line_color, lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        ax.plot([n, 1.05 * n], [pos.xys[i], pos.date_xy[i]], color=line_color, lw=1,
                                clip_on=False)
                        ax.plot([1.05 * n, name_xs[i]], [pos.date_xy[i], pos.date_xy[i]], color=line_color, lw=1,
                                clip_on=False)
                    ax.text(name_xs[i], pos.date_xy[i], f'{pos.genes[i]}',
                            va='center', ha='left', fontsize=8)
            if color_date == "True":
                if len(pos.xys) == 1:
                    # 单个名称直接连接
                    if pos.gene_name in date_colors.keys():
                        ax.plot([1.05 * n, name_xs[0]], [pos.xy, pos.xy], color=date_colors[pos.gene_name], lw=1,
                                clip_on=False)
                    else:
                        ax.plot([1.05 * n, name_xs[0]], [pos.xy, pos.xy], color="#808080", lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        # 从垂直线到名称的水平线
                        if pos.genes[i] in date_colors.keys():
                            ax.plot([n, 1.05 * n], [pos.xys[i], pos.date_xy[i]], color=date_colors[pos.genes[i]],
                                    lw=1,
                                    clip_on=False)
                            ax.plot([1.05 * n, name_xs[i]], [pos.date_xy[i], pos.date_xy[i]],
                                    color=date_colors[pos.genes[i]], lw=1,
                                    clip_on=False)
                        else:
                            ax.plot([n, 1.05 * n], [pos.xys[i], pos.date_xy[i]], color="#808080",
                                    lw=1,
                                    clip_on=False)
                            ax.plot([1.05 * n, name_xs[i]], [pos.date_xy[i], pos.date_xy[i]],
                                    color="#808080", lw=1,
                                    clip_on=False)
                    ax.text(name_xs[i], pos.date_xy[i], f'{pos.genes[i]}',va='center', ha='left', fontsize=8)
def tree_x_gene_1(x_gene,mod,m_2,color_date,colors,ax,date_colors,m):
    for idx, pos in enumerate(x_gene):
        if mod == "torder":
            has_overlap = False
            for j, pos2 in enumerate(x_gene):
                if idx == j:
                    continue
                start1, end1 = min(pos.date_xy), max(pos.date_xy)
                start2, end2 = min(pos2.date_xy), max(pos2.date_xy)

                # 判断是否有重叠
                if not (end1 < start2 or end2 < start1):
                    # 计算名称坐标（自动防重叠）
                    has_overlap = True
                    break  # 只要有一个重叠，就可以停止检查
            # 计算名称坐标（自动防重叠）
            if has_overlap:
                name_ys = [-0.07 * m_2 - i * longs for i in range(len(pos.genes))]
            else:
                name_ys = [-0.07 * m_2 for i in range(len(pos.genes))]

            if color_date == "False":
                line_color = colors[
                    idx % len(colors)]  # 取模运算，确保索引 idx 始终落在 colors 的有效范围内  当 idx >= len(colors) 时，循环回到列表开头
                if len(pos.xys) == 1:
                    ax.plot([pos.xy, pos.xy], [-0.05 * m_2, name_ys[0]], color=line_color, lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        ax.plot([pos.xys[i], pos.date_xy[i]], [0, -0.05 * m_2], color=line_color, lw=1,
                                clip_on=False)
                        ax.plot([pos.date_xy[i], pos.date_xy[i]], [-0.05 * m_2, name_ys[i]], color=line_color, lw=1,
                                clip_on=False)
                    ax.text(pos.date_xy[i], name_ys[i], f'{pos.genes[i]}',
                            va='bottom', ha='center', fontsize=8, rotation=90)
            if color_date == "True":
                if len(pos.xys) == 1:
                    if pos.gene_name in date_colors.keys():
                        ax.plot([pos.xy, pos.xy], [-0.05 * m_2, name_ys[0]], color=date_colors[pos.gene_name], lw=1,
                                clip_on=False)
                    else:
                        ax.plot([pos.xy, pos.xy], [-0.05 * m_2, name_ys[0]], color="#808080", lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        if pos.genes[i] in date_colors.keys():
                            ax.plot([pos.xys[i], pos.date_xy[i]], [0, -0.05 * m_2], color=date_colors[pos.genes[i]],
                                    lw=1,
                                    clip_on=False)
                            ax.plot([pos.date_xy[i], pos.date_xy[i]], [-0.05 * m_2, name_ys[i]],
                                    color=date_colors[pos.genes[i]], lw=1,
                                    clip_on=False)
                        else:
                            ax.plot([pos.xys[i], pos.date_xy[i]], [0, -0.05 * m_2], color="#808080", lw=1,
                                    clip_on=False)
                            ax.plot([pos.date_xy[i], pos.date_xy[i]], [-0.05 * m_2, name_ys[i]], color="#808080",
                                    lw=1,
                                    clip_on=False)
                    ax.text(pos.date_xy[i], name_ys[i], f'{pos.genes[i]}',
                            va='bottom', ha='center', fontsize=8, rotation=90)

        if mod == "tpos":
            has_overlap = False
            for j, pos2 in enumerate(x_gene):
                if idx == j:
                    continue
                start1, end1 = min(pos.date_xy), max(pos.date_xy)
                start2, end2 = min(pos2.date_xy), max(pos2.date_xy)

                # 判断是否有重叠
                if not (end1 < start2 or end2 < start1):
                    # 计算名称坐标（自动防重叠）
                    has_overlap = True
                    break  # 只要有一个重叠，就可以停止检查
            # 计算名称坐标（自动防重叠）
            if has_overlap:
                name_ys = [-0.07 * m - i * longs for i in range(len(pos.genes))]
            else:
                name_ys = [-0.07 * m for i in range(len(pos.genes))]

            if color_date == "False":
                line_color = colors[
                    idx % len(colors)]  # 取模运算，确保索引 idx 始终落在 colors 的有效范围内  当 idx >= len(colors) 时，循环回到列表开头
                if len(pos.xys) == 1:
                    ax.plot([pos.xy, pos.xy], [-0.05 * m, name_ys[0]], color=line_color, lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        ax.plot([pos.xys[i], pos.date_xy[i]], [0, -0.05 * m], color=line_color, lw=1,
                                clip_on=False)
                        ax.plot([pos.date_xy[i], pos.date_xy[i]], [-0.05 * m, name_ys[i]], color=line_color, lw=1,
                                clip_on=False)
                    ax.text(pos.date_xy[i], name_ys[i], f'{pos.genes[i]}',
                            va='bottom', ha='center', fontsize=8, rotation=90)
            if color_date == "True":
                if len(pos.xys) == 1:
                    if pos.gene_name in date_colors.keys():
                        ax.plot([pos.xy, pos.xy], [-0.05 * m, name_ys[0]], color=date_colors[pos.gene_name], lw=1,
                                clip_on=False)
                    else:
                        ax.plot([pos.xy, pos.xy], [-0.05 * m, name_ys[0]], color="#808080", lw=1, clip_on=False)
                for i in range(len(pos.xys)):
                    if len(pos.xys) > 1:
                        if pos.genes[i] in date_colors.keys():
                            ax.plot([pos.xys[i], pos.date_xy[i]], [0, -0.05 * m], color=date_colors[pos.genes[i]],
                                    lw=1,
                                    clip_on=False)
                            ax.plot([pos.date_xy[i], pos.date_xy[i]], [-0.05 * m, name_ys[i]],
                                    color=date_colors[pos.genes[i]], lw=1,
                                    clip_on=False)
                        else:
                            ax.plot([pos.xys[i], pos.date_xy[i]], [0, -0.05 * m], color="#808080", lw=1,
                                    clip_on=False)
                            ax.plot([pos.date_xy[i], pos.date_xy[i]], [-0.05 * m, name_ys[i]], color="#808080",
                                    lw=1,
                                    clip_on=False)
                    ax.text(pos.date_xy[i], name_ys[i], f'{pos.genes[i]}',
                            va='bottom', ha='center', fontsize=8, rotation=90)
def ha_f_set(ax,m,n_2,n,m_2):
    ax2 = ax.twinx()  # 添加右侧y轴
    tick_interval = 50 * 1000000  # 50M的固定间隔
    max_tick = min(int(m // tick_interval) * tick_interval, m)  # 确保不超过m
    num_ticks = int(max_tick / tick_interval) + 1  # 计算刻度数量
    # 生成刻度位置（从0开始，50M间隔）
    tick_positions = np.linspace(0, max_tick, num_ticks)
    # 设置y轴范围（上限设为m，但刻度只到不超过m的最后一个50M点）
    ax2.set_ylim(0, m)
    # 生成刻度标签（0不加M，其他加M）
    tick_labels = [f"{int(pos / 1000000)}M" if pos else "0" for pos in tick_positions]
    ax2.invert_yaxis()
    ax2.xaxis.set_visible(False)  # 隐藏 x 轴
    # 隐藏 ax2 的所有轴线（上、下、左、右）
    ax2.spines['top'].set_visible(False)  # 隐藏上边框
    ax2.spines['right'].set_visible(False)  # 隐藏右边框
    ax2.spines['left'].set_visible(False)  # 隐藏左边框
    ax2.spines['bottom'].set_visible(False)  # 隐藏下边框
    ax2.yaxis.set_visible(True)  # 确保 y 轴刻度可见
    # ax2.yaxis.tick_left()  # 刻度移到左侧
    # ax2.yaxis.set_label_position("left")  # 标签移到左侧
    ax2.set_yticks(tick_positions)
    ax2.set_yticklabels(tick_labels)
    # 创建顶部 x 轴 (ax3)
    ax3 = ax.twiny()  # 使用 twiny() 而不是 twinx()
    # 保持相同的刻度计算逻辑
    tick_interval = 50 * 1000000
    max_tick = min(int(n // tick_interval) * tick_interval, n)
    num_ticks = int(max_tick / tick_interval) + 1
    tick_positions = np.linspace(0, max_tick, num_ticks)
    # 设置 x 轴范围
    ax3.set_xlim(0, n)
    # ax3.xaxis.set_ticks_position('bottom')  # 将刻度放在底部
    # ax3.xaxis.set_label_position('bottom')  # 将标签也放在底部
    # 生成刻度标签
    tick_labels = [f"{int(pos / 1e6)}M" if pos else "0" for pos in tick_positions]
    # 隐藏不需要的元素
    ax3.yaxis.set_visible(False)  # 隐藏 y 轴
    ax3.spines['top'].set_visible(False)  # 隐藏上边框
    ax3.spines['right'].set_visible(False)  # 隐藏右边框
    ax3.spines['left'].set_visible(False)  # 隐藏左边框
    ax3.spines['bottom'].set_visible(False)  # 隐藏下边框
    # 设置刻度
    ax3.set_xticks(tick_positions)
    ax3.set_xticklabels(tick_labels)
    # 如果需要反转刻度方向
    # ax3.invert_xaxis()
    if mod == "tpos":
        ax.set_xlim(0, n)  # 设置x轴的显示范围
        ax.set_ylim(m, 0)  # 因为y轴已经反转  设置y轴的显示范围
    elif mod == "torder":
        ax.set_xlim(0, n_2)  # 设置x轴的显示范围
        ax.set_ylim(m_2, 0)  # 因为y轴已经反转  设置y轴的显示范围
    # 移除默认的刻度标签
    ax.set_xticks([])
    ax.set_yticks([])
def rate_set(alignment_id1_date,ax_top,alignment_id2_date,ax_right,name_3,name_4):
    values = []
    for key, value in alignment_id1_date.items():
        values.extend([key] * value)  # 假设每个key出现value次（模拟频数）

    # 绘制直方图（使用 hist 而不是 bar）
    ax_top.hist(values,
                bins=[0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2],  # 5个柱子对应x=1,2,3,4,5
                range=(1, 6),  # 修正：range应包含右边界
                edgecolor='white',
                linewidth=0.5,
                alpha=1,
                color="#FFCD16")

    # 设置坐标轴范围和刻度
    ax_top.set(xlim=(0.5, 5.5),  # 调整为让柱子居中显示
               ylim=(0, 100))  # y轴保持0-100
    ax_top.set_xticks([1, 2, 3, 4, 5])  # 明确设置x轴刻度
    ax_top.set_yticks(np.arange(0, 101, 20))  # y轴刻度0,20,40,...,100

    # 隐藏上边框和右边框
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)

    # 设置标题
    ax_top.set_title(fr"$\boldsymbol{{{name_3}}}$", fontsize=14, pad=20)

    values = []
    for y_val, count in alignment_id2_date.items():
        values.extend([y_val] * count)  # 每个y值重复count次
    # 绘制水平直方图
    ax_right.hist(values,
                  bins=[0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2],  # 使柱子正好对应y=1,2,3,4,5
                  range=(0.5, 5.5),
                  orientation='horizontal',  # 水平方向
                  edgecolor='white',
                  linewidth=0.5,
                  color="#FFCD16",
                  alpha=1)
    # 坐标轴设置
    ax_right.set(xlim=(0, 100),  # x轴范围保持0-100
                 ylim=(0.5, 5.5))  # y轴范围调整使柱子居中
    ax_right.set_xticks(np.arange(0, 101, 20))  # x轴刻度
    ax_right.set_yticks([1, 2, 3, 4, 5])  # y轴刻度

    # 刻度位置和标签旋转
    ax_right.xaxis.set_ticks_position('top')
    ax_right.tick_params(axis='x', labelrotation=270)
    ax_right.tick_params(axis='y', labelrotation=270)
    ax_right.invert_yaxis()  # 保持y轴反转

    # 标题（使用text模拟原位置）
    ax_right.text(120, 3, fr"$\boldsymbol{{{name_4}}}$",
                  ha='center', va='center',
                  fontsize=9, rotation=270)

    # 边框控制
    ax_right.spines['bottom'].set_visible(False)
    ax_right.spines['right'].set_visible(False)

    for key, value in alignment_id2_date.items():
        ax_right.text(value + 4, key, str(f"{value}%"), ha='center', va='center', fontsize=9,
                      rotation=270)  # 显示数字
    for key, value in alignment_id1_date.items():
        ax_top.text(key, value, str(f"{value}%"), ha='center', va='bottom', fontsize=9)  # 显示数字
def ks_dot(dot_red_1,dot_red_2,dot_c):
    dot_cs = []
    reds_1 = []
    reds_2 = []
    for i in range(len(dot_red_1)):
        if dot_red_2[i] >= dot_red_1[i]:
            reds_1.append(dot_red_1[i])
            reds_2.append(dot_red_2[i])
            dot_cs.append(dot_c[i])
    return dot_cs,reds_1,reds_2
def set_ks_1(alignment_id,mod):
    alignment_ids = []
    for i in range(len(alignment_id)):
        if mod == "torder":
            if alignment_id[i].y_end1 >= alignment_id[i].x_start1:
                alignment_ids.append(alignment_id[i])

        if mod == "tpos":
            if alignment_id[i].y_end2 >= alignment_id[i].x_start2:
                alignment_ids.append(alignment_id[i])
    return alignment_ids
def set_ks_2(alignment_ids, ax, Ks, mod, ks_mod, NSD):
    i = 0
    if Ks == "True":
        for alignment in alignment_ids:
            if mod == "torder":
                if ks_mod == "Average":
                    if alignment.ns > NSD:
                        ax.annotate(f"{alignment.ks1:.2f}",xy=(alignment.x_start1, alignment.y_end1),
                                    xytext=(0, 0),  # 设置文本距离目标点的偏移量（单位是像素）
                                    textcoords='offset points',fontsize=5,  # 表示 xytext 是相对于点的偏移量
                                    )
                elif ks_mod == "Median":
                    if alignment.ns > NSD:
                        ax.annotate(f"{alignment.ks2:.2f}",xy=(alignment.x_start1, alignment.y_end1),
                                    xytext=(0, 0),  # 设置文本距离目标点的偏移量（单位是像素）
                                    textcoords='offset points',fontsize=5,  # 表示 xytext 是相对于点的偏移量
                                    )
            elif mod == "tpos":
                if ks_mod == "Average":
                    if alignment.ns > NSD:
                        ax.annotate(f"{alignment.ks1:.2f}",
                                    xy=(alignment.x_start2, alignment.y_end2),xytext=(0, 0),  # 设置文本距离目标点的偏移量（单位是像素）
                                    textcoords='offset points',fontsize=5,  # 表示 xytext 是相对于点的偏移量
                                    )
                elif ks_mod == "Median":
                    if alignment.ns > NSD:
                        ax.annotate(f"{alignment.ks2:.2f}",xy=(alignment.x_start2, alignment.y_end2),
                                    xytext=(0, 0),  # 设置文本距离目标点的偏移量（单位是像素）
                                    textcoords='offset points',fontsize=5,  # 表示 xytext 是相对于点的偏移量
                                    )
#前置处理部分
gff_1=read_gff(blast_gff_1)
gff_2=read_gff(blast_gff_2)
gff_1["tpos"] = (gff_1["tpos_start"] + gff_1["tpos_end"]) / 2
gff_2["tpos"] = (gff_2["tpos_start"] + gff_2["tpos_end"]) / 2
gff_chr_1 = gff_1.set_index("tnew_name")["chr"].to_dict()  # 得到基因对应的染色体
gff_chr_2 = gff_2.set_index("tnew_name")["chr"].to_dict()  # 得到基因对应的染色体
gff_dict_1 = gff_1.set_index("tnew_name")["tpos"].to_dict()  # 得到基因对应的染色体实际位置
gff_dict_2 = gff_2.set_index("tnew_name")["tpos"].to_dict()  # 得到基因对应的染色体实际位置
gff_torder_1 = gff_1.set_index("tnew_name")["torder"].to_dict()  # 得到基因对应的染色体相对位置
gff_torder_2 = gff_2.set_index("tnew_name")["torder"].to_dict()  # 得到基因对应的染色体相对位置
x_lens ,y_lens,x_len,y_len,x_lens_2,y_lens_2,x_len_2,y_len_2= [],[],[],[],[],[],[],[]
x_len.append(0)
y_len.append(0)
x_len_2.append(0)
y_len_2.append(0)
n, n_2,m,m_2= 0,0,0,0
chr_dict_1,chr_dict_2,chr_torder_1,chr_torder_2 = {},{},{},{}
long1=read_lens(blast_lens_1)
long1_1=read_lens(blast_lens_2)
chr_lens_1 = long1.set_index("chr")["tpos"].to_dict()
chr_lens_2 = long1_1.set_index("chr")["tpos"].to_dict()
lens_torder_1 = long1.set_index("chr")["torder"].to_dict()
lens_torder_2 = long1_1.set_index("chr")["torder"].to_dict()# 构建字典   染色体 : 染色体大小    得到2种形式    tpos  torder
avg_xdate,avg_ydate,avg_xl,avg_yl = [],[],[],[]
for x in long1["chr"]:# 构建字典   x轴 染色体 : 染色体画线位置     染色体 : 染色体标签位置
    chr_dict_1[x] = n
    n = n + chr_lens_1[x]
    z = chr_dict_1[x] + 0.5 * chr_lens_1[x]
    x_lens.append(z)  # 染色体标签位置    tpos
    x_len.append(n)  # 画线位置         tpos
    chr_torder_1[x] = n_2
    n_2 = n_2 + lens_torder_1[x]
    z_2 = chr_torder_1[x] + 0.5 * lens_torder_1[x]
    x_lens_2.append(z_2)  # 染色体标签位置    torder
    x_len_2.append(n_2)  # 画线位置         torder
for x in long1_1["chr"]:# 构建字典   y轴 染色体 : 染色体画线位置     染色体 : 染色体标签位置
    chr_dict_2[x] = m
    m = m + chr_lens_2[x]
    z = chr_dict_2[x] + 0.5 * chr_lens_2[x]
    y_lens.append(z)  # 染色体标签位置    tpos
    y_len.append(m)  # 画线位置         tpos
    chr_torder_2[x] = m_2
    m_2 = m_2 + lens_torder_2[x]
    z_2 = chr_torder_2[x] + 0.5 * lens_torder_2[x]
    y_lens_2.append(z_2)  # 染色体标签位置    torder
    y_len_2.append(m_2)  # 画线位置         torder
if half == "self_col_2" or half == "self_blco_2":
    gff_1_4=read_gff(blast_gff_1_2)
    gff_2_4=read_gff(blast_gff_2_2)
    gff_1_4["tpos"] = (gff_1_4["tpos_start"] + gff_1_4["tpos_end"]) / 2
    gff_2_4["tpos"] = (gff_2_4["tpos_start"] + gff_2_4["tpos_end"]) / 2
    gff_chr_1_4 = gff_1_4.set_index("tnew_name")["chr"].to_dict()  # 得到基因对应的染色体
    gff_chr_2_4 = gff_2_4.set_index("tnew_name")["chr"].to_dict()  # 得到基因对应的染色体
    gff_dict_1_4 = gff_1_4.set_index("tnew_name")["tpos"].to_dict()  # 得到基因对应的染色体实际位置
    gff_dict_2_4 = gff_2_4.set_index("tnew_name")["tpos"].to_dict()  # 得到基因对应的染色体实际位置
    gff_torder_1_4 = gff_1_4.set_index("tnew_name")["torder"].to_dict()  # 得到基因对应的染色体相对位置
    gff_torder_2_4 = gff_2_4.set_index("tnew_name")["torder"].to_dict()  # 得到基因对应的染色体相对位置
    x_lens_4,y_lens_4,x_len_4,y_len_4,x_lens_2_4,y_lens_2_4,x_len_2_4,y_len_2_4 = [],[],[],[],[],[],[],[]
    x_len_4.append(0)
    y_len_4.append(0)
    x_len_2_4.append(0)
    y_len_2_4.append(0)
    n_4,n_2_4,m_4,m_2_4 = 0,0,0,0
    chr_dict_1_4,chr_dict_2_4,chr_torder_1_4,chr_torder_2_4 = {},{},{},{}
    long1_4=read_lens(blast_lens_1_2)
    long1_1_4=read_lens(blast_lens_2_2)
    chr_lens_1_4 = long1_4.set_index("chr")["tpos"].to_dict()
    chr_lens_2_4 = long1_1_4.set_index("chr")["tpos"].to_dict()
    lens_torder_1_4 = long1_4.set_index("chr")["torder"].to_dict()
    lens_torder_2_4 = long1_1_4.set_index("chr")["torder"].to_dict()# 构建字典   染色体 : 染色体大小    得到2种形式    tpos  torder
    avg_xdate_4,avg_ydate_4,avg_xl_4,avg_yl_4 = [],[],[],[]
    for x in long1_4["chr"]:# 构建字典   y轴 染色体 : 染色体画线位置     染色体 : 染色体标签位置
        chr_dict_1_4[x] = n_4
        n_4 = n_4 + chr_lens_1_4[x]
        z_4 = chr_dict_1_4[x] + 0.5 * chr_lens_1_4[x]
        x_lens_4.append(z_4)  # 染色体标签位置    tpos
        x_len_4.append(n_4)  # 画线位置         tpos
        chr_torder_1_4[x] = n_2_4
        n_2_4 = n_2_4 + lens_torder_1_4[x]
        z_2_4 = chr_torder_1_4[x] + 0.5 * lens_torder_1_4[x]
        x_lens_2_4.append(z_2_4)  # 染色体标签位置    torder
        x_len_2_4.append(n_2_4)  # 画线位置         torder
    for x in long1_1_4["chr"]:# 构建字典   y轴 染色体 : 染色体画线位置     染色体 : 染色体标签位置
        chr_dict_2_4[x] = m_4
        m_4 = m_4 + chr_lens_2_4[x]
        z_4 = chr_dict_2_4[x] + 0.5 * chr_lens_2_4[x]
        y_lens_4.append(z_4)  # 染色体标签位置    tpos
        y_len_4.append(m_4)  # 画线位置         tpos
        chr_torder_2_4[x] = m_2_4
        m_2_4 = m_2_4 + lens_torder_2_4[x]
        z_2_4 = chr_torder_2_4[x] + 0.5 * lens_torder_2_4[x]
        y_lens_2_4.append(z_2_4)  # 染色体标签位置    torder
        y_len_2_4.append(m_2_4)  # 画线位置         torder
#主要部分
if mods == "blast":
    new_blasts_1=read_blast(blast_name,gff_1,gff_2,gff_chr_1,chr_lens_1,gff_chr_2,chr_lens_2)
    new_blasts_1['count'] = new_blasts_1.groupby("query_id").cumcount() + 1
    # 分配匹配优先级 - 保持原变量名但改变筛选逻辑
    dot_2 = new_blasts_1[new_blasts_1['count'] == 1]  # 第1次匹配 -> 红色点 (原dot_2)
    dot_3 = new_blasts_1[new_blasts_1['count'] == 2]  # 第2次匹配 -> 蓝色点 (原dot_3)
    dot_7 = new_blasts_1[(new_blasts_1['count'] >= 3) & (new_blasts_1['count'] <= 5)]  # 第3次及以上 -> 灰色点 (相当于原dot_4+dot_5+dot_6合并)
    if mod == "torder":  # 得到  torder   情况下的所有点的  xy  坐标
        dot_red_1, dot_red_2 = bla_dot(dot_2, gff_chr_1, gff_torder_1, chr_torder_1, gff_chr_2, gff_torder_2,chr_torder_2)
        dot_blue_1, dot_blue_2 = bla_dot(dot_3, gff_chr_1, gff_torder_1, chr_torder_1, gff_chr_2, gff_torder_2,chr_torder_2)
        dot_grey_1, dot_grey_2 = bla_dot(dot_7, gff_chr_1, gff_torder_1, chr_torder_1, gff_chr_2, gff_torder_2,chr_torder_2)
    elif mod == "tpos":  # 得到  tpos   情况下的所有点的  xy  坐标
        dot_red_1, dot_red_2 = bla_dot(dot_2, gff_chr_1, gff_dict_1, chr_dict_1, gff_chr_2, gff_dict_2,chr_dict_2)
        dot_blue_1, dot_blue_2 = bla_dot(dot_3, gff_chr_1, gff_dict_1, chr_dict_1, gff_chr_2, gff_dict_2,chr_dict_2)
        dot_grey_1, dot_grey_2 = bla_dot(dot_7, gff_chr_1, gff_dict_1, chr_dict_1, gff_chr_2, gff_dict_2,chr_dict_2)
    if half == "self_blco_2":
        blast_l_4 = read_collinearity(Collinearity_2)
        if collinearity_ks_mod == "False":
            if mod == "torder":  # 得到点的坐标
                dot_red_1_4, dot_red_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = order_col_dot(blast_l_4, gff_chr_1_4,gff_chr_2_4,gff_torder_1_4,gff_torder_2_4,chr_torder_1_4,chr_torder_2_4, name_id)
            elif mod == "tpos":
                dot_red_1_4, dot_red_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = tpos_col_dot(blast_l_4, gff_chr_1_4,gff_chr_2_4, gff_torder_1_4,gff_torder_2_4,gff_dict_1_4, chr_dict_1_4,gff_dict_2_4, chr_dict_2_4,name_id)
        if collinearity_ks_mod == "True":  # 读取 Collinearity_ks 文件情况下
            final_result_4 = read_ks(Collinearity_ks_2, blast_l_4)
            if mod == "torder":  # 得到点的坐标
                if top_block == "False":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = order_ks_dot(final_result_4, gff_chr_1_4, gff_chr_2_4, gff_torder_1_4, gff_torder_2_4, chr_torder_1_4,chr_torder_2_4, top_block, ks_mod, ks_dem_2, ks_dem_2_2, name_id)
                if top_block == "True":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_block_1_4, dot_block_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = order_ks_dot(final_result_4, gff_chr_1_4, gff_chr_2_4, gff_torder_1_4, gff_torder_2_4, chr_torder_1_4,chr_torder_2_4, top_block, ks_mod, ks_dem_2, ks_dem_2_2, name_id)
            elif mod == "tpos":
                if top_block == "False":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = tpos_ks_dot(final_result_4, gff_chr_1_4, gff_chr_2_4, gff_torder_1_4, gff_torder_2_4, gff_dict_1_4,chr_dict_1_4, gff_dict_2_4, chr_dict_2_4, top_block, ks_mod, ks_dem_2, ks_dem_2_2, name_id)
                if top_block == "True":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_block_1_4, dot_block_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = tpos_ks_dot(final_result_4, gff_chr_1_4, gff_chr_2_4, gff_torder_1_4, gff_torder_2_4, gff_dict_1_4,chr_dict_1_4, gff_dict_2_4, chr_dict_2_4, top_block, ks_mod, ks_dem_2, ks_dem_2_2, name_id)
            alignment_id_4 = get_block(final_result_4, gff_chr_1_4, gff_chr_2_4, gff_torder_1_4, gff_torder_2_4,chr_torder_1_4, gff_dict_1_4, chr_dict_1_4, chr_torder_2_4, gff_dict_2_4,chr_dict_2_4, name_id, mod, Distance, N)
    if half == "self_col_2":
        new_blasts_1_4=read_blast(blast_name_2,gff_1_4,gff_2_4,gff_chr_1_4,chr_lens_1_4,gff_chr_2_4,chr_lens_2_4)
        new_blasts_1_4['count'] = new_blasts_1_4.groupby("query_id").cumcount() + 1
        # 分配匹配优先级 - 保持原变量名但改变筛选逻辑
        dot_2_4 = new_blasts_1_4[new_blasts_1_4['count'] == 1]  # 第1次匹配 -> 红色点 (原dot_2)
        dot_3_4 = new_blasts_1_4[new_blasts_1_4['count'] == 2]  # 第2次匹配 -> 蓝色点 (原dot_3)
        dot_7_4 = new_blasts_1_4[(new_blasts_1_4['count'] >= 3) & (new_blasts_1_4['count'] <= 5)]  # 第3次及以上 -> 灰色点 (相当于原dot_4+dot_5+dot_6合并)
        if mod == "torder":  # 得到  torder   情况下的所有点的  xy  坐标
            dot_red_1_4, dot_red_2_4 = bla_dot(dot_2_4, gff_chr_1_4, gff_torder_1_4, chr_torder_1_4, gff_chr_2_4, gff_torder_2_4,chr_torder_2_4)
            dot_blue_1_4, dot_blue_2_4 = bla_dot(dot_3_4, gff_chr_1_4, gff_torder_1_4, chr_torder_1_4, gff_chr_2_4, gff_torder_2_4,chr_torder_2_4)
            dot_grey_1_4, dot_grey_2_4 = bla_dot(dot_7_4, gff_chr_1_4, gff_torder_1_4, chr_torder_1_4, gff_chr_2_4, gff_torder_2_4,chr_torder_2_4)
        elif mod == "tpos":  # 得到  tpos   情况下的所有点的  xy  坐标
            dot_red_1_4, dot_red_2_4 = bla_dot(dot_2_4, gff_chr_1_4, gff_dict_1_4, chr_dict_1_4, gff_chr_2_4, gff_dict_2_4,chr_dict_2_4)
            dot_blue_1_4, dot_blue_2_4 = bla_dot(dot_3_4, gff_chr_1_4, gff_dict_1_4, chr_dict_1_4, gff_chr_2_4, gff_dict_2_4,chr_dict_2_4)
            dot_grey_1_4, dot_grey_2_4 = bla_dot(dot_7_4, gff_chr_1_4, gff_dict_1_4, chr_dict_1_4, gff_chr_2_4, gff_dict_2_4,chr_dict_2_4)
if mods == "collinearity":  # 只读取collinearity文件情况下
    blast_l = read_collinearity(Collinearity)
    if half == "self_col_2":
        blast_l_4 = read_collinearity(Collinearity_2)
    if collinearity_ks_mod == "False":
        if mod == "torder":  # 得到点的坐标
            dot_red_1, dot_red_2, dot_star_1, dot_end_1, dot_len_1=order_col_dot(blast_l,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,chr_torder_1,chr_torder_2,name_id)
        elif mod == "tpos":
            dot_red_1, dot_red_2, dot_star_1, dot_end_1, dot_len_1 = tpos_col_dot(blast_l,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,gff_dict_1,chr_dict_1,gff_dict_2,chr_dict_2,name_id)
        if half == "self_col_2":
            if mod == "torder":  # 得到点的坐标
                dot_red_1_4, dot_red_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = order_col_dot(blast_l_4, gff_chr_1_4, gff_chr_2_4,gff_torder_1_4, gff_torder_2_4,chr_torder_1_4, chr_torder_2_4,name_id)
            elif mod == "tpos":
                dot_red_1_4, dot_red_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = tpos_col_dot(blast_l_4, gff_chr_1_4, gff_chr_2_4,gff_torder_1_4, gff_torder_2_4,gff_dict_1_4, chr_dict_1_4,gff_dict_2_4, chr_dict_2_4, name_id)
        if mod_show == "rate":
            alignment_id1_date, alignment_id2_date=rate_date(blast_l,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,name_id)
    if collinearity_ks_mod == "True":  # 读取 Collinearity_ks 文件情况下
        final_result=read_ks(Collinearity_ks,blast_l)
        if mod == "torder":  # 得到点的坐标
            if top_block == "False":
                dot_red_1, dot_red_2, dot_c, dot_star_1, dot_end_1, dot_len_1=order_ks_dot(final_result,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,chr_torder_1,chr_torder_2,top_block,ks_mod,ks_dem_1,ks_dem_1_1,name_id)
            if top_block == "True":
                dot_red_1, dot_red_2, dot_c, dot_block_1, dot_block_2, dot_star_1, dot_end_1, dot_len_1=order_ks_dot(final_result,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,chr_torder_1,chr_torder_2,top_block,ks_mod,ks_dem_1,ks_dem_1_1,name_id)
        elif mod == "tpos":
            if top_block == "False":
                dot_red_1, dot_red_2, dot_c, dot_star_1, dot_end_1, dot_len_1=tpos_ks_dot(final_result,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,gff_dict_1,chr_dict_1,gff_dict_2,chr_dict_2,top_block,ks_mod,ks_dem_1,ks_dem_1_1,name_id)
            if top_block == "True":
                dot_red_1, dot_red_2, dot_c, dot_block_1, dot_block_2, dot_star_1, dot_end_1, dot_len_1=tpos_ks_dot(final_result,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,gff_dict_1,chr_dict_1,gff_dict_2,chr_dict_2,top_block,ks_mod,ks_dem_1,ks_dem_1_1,name_id)
        if half == "self_col_2":
            final_result_4=read_ks(Collinearity_ks_2,blast_l_4)
            if mod == "torder":  # 得到点的坐标
                if top_block == "False":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = order_ks_dot(final_result_4,gff_chr_1_4, gff_chr_2_4,gff_torder_1_4,gff_torder_2_4,chr_torder_1_4,chr_torder_2_4,top_block, ks_mod,ks_dem_2, ks_dem_2_2,name_id)
                if top_block == "True":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_block_1_4, dot_block_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = order_ks_dot(final_result_4,gff_chr_1_4, gff_chr_2_4,gff_torder_1_4,gff_torder_2_4,chr_torder_1_4,chr_torder_2_4,top_block, ks_mod,ks_dem_2, ks_dem_2_2,name_id)
            elif mod == "tpos":
                if top_block == "False":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = tpos_ks_dot(final_result_4, gff_chr_1_4,gff_chr_2_4, gff_torder_1_4,gff_torder_2_4,gff_dict_1_4, chr_dict_1_4,gff_dict_2_4, chr_dict_2_4,top_block, ks_mod,ks_dem_2, ks_dem_2_2,name_id)
                if top_block == "True":
                    dot_red_1_4, dot_red_2_4, dot_c_4, dot_block_1_4, dot_block_2_4, dot_star_1_4, dot_end_1_4, dot_len_1_4 = tpos_ks_dot(final_result_4, gff_chr_1_4,gff_chr_2_4, gff_torder_1_4,gff_torder_2_4,gff_dict_1_4, chr_dict_1_4,gff_dict_2_4, chr_dict_2_4,top_block, ks_mod,ks_dem_2, ks_dem_2_2,name_id)
        if mod_show == "rate":  # 显示 基因对比 的比例情况下
            alignment_id1_date, alignment_id2_date = rate_date(final_result, gff_chr_1, gff_chr_2, gff_torder_1,gff_torder_2,name_id)
        alignment_id=get_block(final_result,gff_chr_1,gff_chr_2,gff_torder_1,gff_torder_2,chr_torder_1,gff_dict_1,chr_dict_1,chr_torder_2,gff_dict_2,chr_dict_2,name_id,mod,Distance,N)
        if half == "self_col_2":
            alignment_id_4 = get_block(final_result_4,gff_chr_1_4,gff_chr_2_4,gff_torder_1_4,gff_torder_2_4,chr_torder_1_4,gff_dict_1_4,chr_dict_1_4,chr_torder_2_4,gff_dict_2_4,chr_dict_2_4,name_id,mod,Distance,N)
cdict = {
    'red': [
        (0.0, 1.0, 1.0),  # 黄色 (R=1.0)
        (0.5, 1.0, 1.0),  # 橙色 (R=1.0)
        (1.0, 1.0, 1.0)  # 红色 (R=1.0)
    ],
    'green': [
        (0.0, 1.0, 1.0),  # 黄色 (G=1.0)
        (0.5, 0.65, 0.65),  # 橙色 (G=0.65)
        (1.0, 0.0, 0.0)  # 红色 (G=0.0)
    ],
    'blue': [
        (0.0, 0.0, 0.0),  # 黄色 (B=0.0)
        (0.5, 0.0, 0.0),  # 橙色 (B=0.0)
        (1.0, 0.0, 0.0)  # 红色 (B=0.0)
    ]
}  # 定义了从 黄色（Yellow）→ 橙色（Orange）→ 红色（Red） 的颜色渐变规则
# 定义 RGB 通道在关键位置的值，控制颜色渐变
norm = plt.Normalize(vmin=0, vmax=3)  # 控制颜色映射范围 数据归一化（Normalize 对象，控制颜色映射范围）
cmap = LinearSegmentedColormap("YellowOrangeRed", cdict)  # 将 cdict 转换为可用的色图对象 颜色映射 （colormap）
if mods == "blast":
    if mod_show == "density":
        fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
        # 使用gridspec创建3个区域
        gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1],#   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
                               height_ratios=[1, 4], wspace=0.2, hspace=0.2)#  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
        ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
        ax_top = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
        ax_right = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
        get_density(gff_1, long1, gff_2, long1_1, chr_dict_1, gff_dict_1, chr_dict_2, gff_dict_2, bin_sizes, n, m, n_2,m_2, ax_top, ax_right, name_1, name_2)
    if mod_show == "nan":
        if half == "self_col_2" or half == "self_blco_2":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1],height_ratios=[1, 1], wspace=-0.2, hspace=-0.2)
            if posd == "bot_left":
                ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
                ax_4 = fig.add_subplot(gs[0, 1])  # 第 1 行，第 2 列
            if posd == "bot_right":
                ax_4 = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
                ax = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
        else:
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))  # 创建图形和坐标轴
        fig.canvas.draw()
    if mod_show == "tree":
        fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
        # 使用gridspec创建3个区域
        gs = gridspec.GridSpec(2, 2, width_ratios=[10, 1],  #   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
                               height_ratios=[1, 10], wspace=0.2, hspace=0.2)#  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
        ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
    if half == "self_blco_2":
        if collinearity_ks_mod == "False":
            half_col_2_1(ax, ax_4, posd)
            if mod == "tpos":
                half_t_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_4, y_lens_4, m_4, n_4, x_len_4, y_len_4, dot_len_1_4,dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            elif mod == "torder":
                half_t_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_2_4, y_lens_2_4, m_2_4, n_2_4, x_len_2_4, y_len_2_4,dot_len_1_4, dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            greys_1,greys_2=half_dot_1(dot_grey_1,dot_grey_2)
            blues_1,blues_2=half_dot_1(dot_blue_1,dot_blue_2)
            reds_1,reds_2=half_dot_1(dot_red_1,dot_red_2)
            ax.scatter(greys_1, greys_2, **dot_mod_1)
            ax.scatter(blues_1, blues_2, **dot_mod_2)
            ax.scatter(reds_1, reds_2, **dot_mod_3)
            ax.set_xticks([])  # 隐藏x轴刻度
            ax.set_yticks([])  # 隐藏y轴刻度
            reds_1_4, reds_2_4 = half_dot_1(dot_red_1_4, dot_red_2_4)
            ax_4.scatter(reds_1_4, reds_2_4, c='blue', s=0.7, alpha=1, marker='o', edgecolors=None,linewidths=0)  # x轴数据  y轴数据
            # s=1,     点的大小（size）   alpha=1, marker='o',   透明度    点的形状   edgecolors=None, linewidths=0  点边缘颜色（None 表示无边框）  点边缘线宽
            ax_4.set_xticks([])  # 隐藏x轴刻度
            ax_4.set_yticks([])  # 隐藏y轴刻度
        if collinearity_ks_mod == "True":
            half_col_2_1(ax, ax_4, posd)
            if mod == "tpos":
                half_t_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_4, y_lens_4, m_4, n_4, x_len_4, y_len_4, dot_len_1_4,dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            elif mod == "torder":
                half_t_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_2_4, y_lens_2_4, m_2_4, n_2_4, x_len_2_4, y_len_2_4,dot_len_1_4, dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            dot_cs_4, reds_1_4, reds_2_4 = ks_dot(dot_red_1_4, dot_red_2_4, dot_c_4)
            greys_1, greys_2 = half_dot_1(dot_grey_1, dot_grey_2)
            blues_1, blues_2 = half_dot_1(dot_blue_1, dot_blue_2)
            reds_1, reds_2 = half_dot_1(dot_red_1, dot_red_2)
            ax.scatter(greys_1, greys_2, **dot_mod_1)
            ax.scatter(blues_1, blues_2, **dot_mod_2)
            ax.scatter(reds_1, reds_2, **dot_mod_3)
            ax.set_xticks([])  # 隐藏x轴刻度
            ax.set_yticks([])  # 隐藏y轴刻度
            if top_block == "False":
                scatter = ax_4.scatter(reds_1_4, reds_2_4,  # x轴数据  y轴数据
                                       c=dot_cs_4, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
                                       norm=norm, s=1.2,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
                                       alpha=1, marker='o',  # 透明度    点的形状
                                       edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
                                       )  # 画出点
            if top_block == "True":
                scatter = ax_4.scatter(reds_1_4, reds_2_4,  # x轴数据  y轴数据
                                       c='grey', s=1.2, alpha=1, marker='o', edgecolors=None, linewidths=0
                                       # 点边缘颜色（None 表示无边框）     点边缘线宽
                                       )  # 画出点
                dot_blocks_1_4 = []
                dot_blocks_2_4 = []
                for i in range(len(dot_block_1_4)):
                    if dot_block_2_4[i] >= dot_block_1_4[i]:
                        dot_blocks_1_4.append(dot_block_1_4[i])
                        dot_blocks_2_4.append(dot_block_2_4[i])
                scatter = ax_4.scatter(dot_blocks_1_4, dot_blocks_2_4,  # x轴数据  y轴数据
                                       c='r', s=1.2, alpha=1, marker='o', edgecolors=None, linewidths=0
                                       # 点边缘颜色（None 表示无边框）     点边缘线宽
                                       )  # 画出点
            alignment_ids_4=set_ks_1(alignment_id_4,mod)
            if Ks == "True":
                set_ks_2(alignment_ids_4, ax_4, Ks, mod, ks_mod, NSD)
            if top_block == "False":
                if posd == "bot_left":
                    cbar = plt.colorbar(scatter, ax=ax_4, fraction=0.03, pad=0.1)
                    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）
                    # 画出颜色条
                    cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
                    cbar.ax.set_position([0.855, 0.6, 0.03, 0.2])
                if posd == "bot_right":
                    cbar = plt.colorbar(scatter, ax=ax_4, fraction=0.03, pad=0.1, location='left')
                    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）
                    # 画出颜色条
                    cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
                    cbar.ax.set_position([0.12, 0.62, 0.03, 0.2])
                    cbar.ax.yaxis.set_ticks_position('left')
            ax_4.set_xticks([])  # 隐藏x轴刻度
            ax_4.set_yticks([])  # 隐藏y轴刻度
    if half == "self_col_2":
        half_col_2_1(ax, ax_4, posd)
        if mod == "tpos":
            half_t_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
            half_t_2(ax_4, long1_4, long1_1_4, x_lens_4, y_lens_4, x_len_4, y_len_4, m_4, n_4, name_1_4, name_2_4)
        elif mod == "torder":
            half_t_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
            half_t_2(ax_4, long1_4, long1_1_4, x_lens_2_4, y_lens_2_4, x_len_2_4, y_len_2_4, m_2_4, n_2_4, name_1_4, name_2_4)
        greys_1,greys_2=half_dot_1(dot_grey_1,dot_grey_2)
        blues_1,blues_2=half_dot_1(dot_blue_1,dot_blue_2)
        reds_1,reds_2=half_dot_1(dot_red_1,dot_red_2)
        ax.scatter(greys_1, greys_2, **dot_mod_1)
        ax.scatter(blues_1, blues_2, **dot_mod_2)
        ax.scatter(reds_1, reds_2, **dot_mod_3)
        ax.set_xticks([])  # 隐藏x轴刻度
        ax.set_yticks([])  # 隐藏y轴刻度
        greys_1_4, greys_2_4 = half_dot_1(dot_grey_1_4, dot_grey_2_4)
        blues_1_4, blues_2_4 = half_dot_1(dot_blue_1_4, dot_blue_2_4)
        reds_1_4, reds_2_4 = half_dot_1(dot_red_1_4, dot_red_2_4)
        ax_4.scatter(greys_1_4, greys_2_4, **dot_mod_1)
        ax_4.scatter(blues_1_4, blues_2_4, **dot_mod_2)
        ax_4.scatter(reds_1_4, reds_2_4, **dot_mod_3)
        ax_4.set_xticks([])  # 隐藏x轴刻度
        ax_4.set_yticks([])  # 隐藏y轴刻度
    if half == "False":
        half_f_1(ax)
        if mod == "tpos":
            half_f_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
        elif mod == "torder":
            half_f_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
        ax.scatter(dot_grey_1, dot_grey_2, **dot_mod_1)
        ax.scatter(dot_blue_1, dot_blue_2, **dot_mod_2)
        ax.scatter(dot_red_1, dot_red_2, **dot_mod_3)
        # "blast"情况下打点
    if half == "True":
        half_t_1(pos, ax)
        if mod == "tpos":
            half_t_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
        elif mod == "torder":
            half_t_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
        greys_1, greys_2 = half_dot_1(dot_grey_1, dot_grey_2)
        blues_1, blues_2 = half_dot_1(dot_blue_1, dot_blue_2)
        reds_1, reds_2 = half_dot_1(dot_red_1, dot_red_2)
        ax.scatter(greys_1, greys_2, **dot_mod_1)
        ax.scatter(blues_1, blues_2, **dot_mod_2)
        ax.scatter(reds_1, reds_2, **dot_mod_3)
        ax.set_xticks([])  # 隐藏x轴刻度
        ax.set_yticks([])  # 隐藏y轴刻度
    if mod_show == "tree":
        x_gene, y_gene=get_tree_1(gene_show,gff_1,gff_chr_1,gff_torder_1,chr_torder_1,gff_dict_1,chr_dict_1,gff_2,gff_chr_2,gff_torder_2,chr_torder_2,gff_dict_2,chr_dict_2,mod)
        x_gene=get_tree_2(x_gene)
        y_gene = get_tree_2(y_gene)
        colors, date_colors=tree_co_1(color_date,color_dates)
        tree_y_gene_1(y_gene, mod, n_2, color_date, colors, ax, date_colors, n)
        tree_x_gene_1(x_gene, mod, m_2, color_date, colors, ax, date_colors, m)
    if half == "False":
        ha_f_set(ax, m, n_2, n, m_2)
if mods == "collinearity":
    if collinearity_ks_mod == "False":
        if mod_show == "nan":
            if half == "self_col_2":
                fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
                # 使用gridspec创建3个区域
                gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1],
                                       height_ratios=[1, 1], wspace=-0.2, hspace=-0.2)
                if posd == "bot_left":
                    ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
                    ax_4 = fig.add_subplot(gs[0, 1])  # 第 1 行，第 2 列# 放大ax子图
                if posd == "bot_right":
                    ax_4 = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
                    ax = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
            else:
                fig, ax = plt.subplots(1, 1, figsize=(8, 8))  # 创建图形和坐标轴
            fig.canvas.draw()
        if mod_show == "density":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1],height_ratios=[1, 4], wspace=0.2, hspace=0.2)
            #   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
            #  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
            ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
            ax_top = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
            ax_right = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
            get_density(gff_1, long1, gff_2, long1_1, chr_dict_1, gff_dict_1, chr_dict_2, gff_dict_2, bin_sizes, n, m,n_2, m_2, ax_top, ax_right, name_1, name_2)
        if mod_show == "rate":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1],height_ratios=[1, 4], wspace=0.2, hspace=0.2)
            #   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
            #  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
            ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
            ax_top = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
            ax_right = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
            rate_set(alignment_id1_date, ax_top, alignment_id2_date, ax_right,name_3,name_4)
        if mod_show == "tree":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[10, 1],height_ratios=[1, 10], wspace=0.2, hspace=0.2)
            #   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
            #  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
            ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
        if half == "False":
            half_f_1(ax)
            if mod == "tpos":
                half_f_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
            elif mod == "torder":
                half_f_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
            ax.scatter(dot_red_1, dot_red_2, c='blue', s=0.7, alpha=1, marker='o', edgecolors=None,linewidths=0)  # x轴数据  y轴数据
            # s=1, 点的大小（size）   alpha=1, marker='o',   透明度    点的形状   edgecolors=None, linewidths=0  点边缘颜色（None 表示无边框）  点边缘线宽
        if half == "self_col_2":
            half_col_2_1(ax, ax_4, posd)
            if mod == "tpos":
                half_col_2_2(ax, long1, long1_1, x_lens, y_lens, m, n, x_len, y_len, dot_len_1, dot_star_1, dot_end_1,name_1, name_2, block)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_4, y_lens_4, m_4, n_4, x_len_4, y_len_4, dot_len_1_4,dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            elif mod == "torder":
                half_col_2_2(ax, long1, long1_1, x_lens_2, y_lens_2, m_2, n_2, x_len_2, y_len_2, dot_len_1, dot_star_1,dot_end_1, name_1, name_2, block)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_2_4, y_lens_2_4, m_2_4, n_2_4, x_len_2_4, y_len_2_4,dot_len_1_4, dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            reds_1,reds_2=half_dot_1(dot_red_1,dot_red_2)
            ax.scatter(reds_1, reds_2, c='blue', s=0.7, alpha=1, marker='o', edgecolors=None,linewidths=0)  # x轴数据  y轴数据
            # s=1,     点的大小（size）   alpha=1, marker='o',   透明度    点的形状   edgecolors=None, linewidths=0  点边缘颜色（None 表示无边框）  点边缘线宽
            ax.set_xticks([])  # 隐藏x轴刻度
            ax.set_yticks([])  # 隐藏y轴刻度
            reds_1_4, reds_2_4 = half_dot_1(dot_red_1_4, dot_red_2_4)
            ax_4.scatter(reds_1_4, reds_2_4, c='blue', s=0.7, alpha=1, marker='o', edgecolors=None,linewidths=0)  # x轴数据  y轴数据
            # s=1,     点的大小（size）   alpha=1, marker='o',   透明度    点的形状   edgecolors=None, linewidths=0  点边缘颜色（None 表示无边框）  点边缘线宽
            ax_4.set_xticks([])  # 隐藏x轴刻度
            ax_4.set_yticks([])  # 隐藏y轴刻度
        if half == "True":
            half_t_1(pos, ax)
            if mod == "tpos":
                half_t_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
            elif mod == "torder":
                half_t_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
            reds_1, reds_2 = half_dot_1(dot_red_1, dot_red_2)
            ax.scatter(reds_1, reds_2, c='blue', s=0.7, alpha=1, marker='o', edgecolors=None,linewidths=0)  # x轴数据  y轴数据
            # s=1,     点的大小（size）   alpha=1, marker='o',   透明度    点的形状   edgecolors=None, linewidths=0  点边缘颜色（None 表示无边框）  点边缘线宽
            ax.set_xticks([])  # 隐藏x轴刻度
            ax.set_yticks([])  # 隐藏y轴刻度
        if mod_show == "tree":
            x_gene, y_gene = get_tree_1(gene_show, gff_1, gff_chr_1, gff_torder_1, chr_torder_1, gff_dict_1, chr_dict_1,gff_2, gff_chr_2, gff_torder_2, chr_torder_2, gff_dict_2, chr_dict_2, mod)
            x_gene = get_tree_2(x_gene)
            y_gene = get_tree_2(y_gene)
            colors, date_colors = tree_co_1(color_date, color_dates)
            tree_y_gene_1(y_gene, mod, n_2, color_date, colors, ax, date_colors, n)
            tree_x_gene_1(x_gene, mod, m_2, color_date, colors, ax, date_colors, m)
        if half == "False":
            ha_f_set(ax, m, n_2, n, m_2)
    if collinearity_ks_mod == "True":
        if mod_show == "nan":
            if half == "self_col_2":
                fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
                # 使用gridspec创建3个区域
                gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1], wspace=-0.2, hspace=-0.2)
                if posd == "bot_left":
                    ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
                    ax_4 = fig.add_subplot(gs[0, 1])  # 第 1 行，第 2 列
                if posd == "bot_right":
                    ax_4 = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
                    ax = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
            else:
                fig, ax = plt.subplots(1, 1, figsize=(8, 8))  # 创建图形和坐标轴
            fig.canvas.draw()
        if mod_show == "density":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1],
                                   height_ratios=[1, 4], wspace=0.2, hspace=0.2)
            #   2, 2   画布划分为2行×2列的网格。
            #  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
            #  height_ratios=[1, 4] 第1行:第2行高度 =1:4
            # wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
            ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
            ax_top = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
            ax_right = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
            get_density(gff_1, long1, gff_2, long1_1, chr_dict_1, gff_dict_1, chr_dict_2, gff_dict_2, bin_sizes, n, m,n_2, m_2, ax_top, ax_right, name_1, name_2)
        if mod_show == "rate":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1],height_ratios=[1, 4], wspace=0.1, hspace=0.2)
            #   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
            #  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
            ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
            ax_top = fig.add_subplot(gs[0, 0])  # 第 1 行，第 1 列
            ax_right = fig.add_subplot(gs[1, 1])  # 第 2 行，第 2 列
            rate_set(alignment_id1_date, ax_top, alignment_id2_date, ax_right, name_3, name_4)
        if mod_show == "tree":
            fig = plt.figure(figsize=(10, 10))  # 整体图形大小调整
            # 使用gridspec创建3个区域
            gs = gridspec.GridSpec(2, 2, width_ratios=[10, 1],height_ratios=[1, 10], wspace=0.2, hspace=0.2)
            #   2, 2   画布划分为2行×2列的网格。#  width_ratios=[4, 1]  第1列:第2列宽度 =4:1
            #  height_ratios=[1, 4] 第1行:第2行高度 =1:4# wspace=0.3, hspace=0.3 子图水平间距垂直间距为 0.3
            ax = fig.add_subplot(gs[1, 0])  # 第 2 行，第 1 列
        if half == "False":
            half_f_1(ax)
            if mod == "tpos":
                half_f_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
            elif mod == "torder":
                half_f_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
            scatter = ax.scatter(dot_red_1, dot_red_2,  # x轴数据  y轴数据
                                 c=dot_c, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
                                 norm=norm, s=1.2,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
                                 alpha=1, marker='o',  # 透明度    点的形状
                                 edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
                                 )  # 画出点
            set_ks_2(alignment_id, ax, Ks, mod, ks_mod, NSD)
            cbar = plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.1)
            # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）
            # 画出颜色条
            cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
            # (0.5, 1.05)：0.5 表示水平居中，1.05 表示稍微高于颜色条。# rotation = 0：标签水平排列（0度旋转）。# transform=cbar.ax.transAxes：使用相对坐标（0-1 范围）。
            # va = 'bottom'（或verticalalignment = 'bottom'）：垂直对齐方式为底部。# ha = 'center'（或horizontalalignment = 'center'）：水平居中对齐。
        if half == "self_col_2":
            half_col_2_1(ax, ax_4, posd)
            if mod == "tpos":
                half_col_2_2(ax, long1, long1_1, x_lens, y_lens, m, n, x_len, y_len, dot_len_1, dot_star_1, dot_end_1,name_1, name_2, block)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_4, y_lens_4, m_4, n_4, x_len_4, y_len_4, dot_len_1_4,dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            elif mod == "torder":
                half_col_2_2(ax, long1, long1_1, x_lens_2, y_lens_2, m_2, n_2, x_len_2, y_len_2, dot_len_1, dot_star_1,dot_end_1, name_1, name_2, block)
                half_col_2_2(ax_4, long1_4, long1_1_4, x_lens_2_4, y_lens_2_4, m_2_4, n_2_4, x_len_2_4, y_len_2_4,dot_len_1_4, dot_star_1_4, dot_end_1_4, name_1_4, name_2_4, block)
            dot_cs, reds_1, reds_2=ks_dot(dot_red_1,dot_red_2,dot_c)
            dot_cs_4, reds_1_4, reds_2_4 = ks_dot(dot_red_1_4, dot_red_2_4, dot_c_4)
            if top_block == "False":
                scatter = ax.scatter(reds_1, reds_2,  # x轴数据  y轴数据
                                     c=dot_cs, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
                                     norm=norm, s=1.2,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
                                     alpha=1, marker='o',  # 透明度    点的形状
                                     edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
                                     )  # 画出点
            if top_block == "True":
                scatter = ax.scatter(reds_1, reds_2,  # x轴数据  y轴数据
                                     c='grey', s=1.2, alpha=1, marker='o', edgecolors=None, linewidths=0
                                     # 点边缘颜色（None 表示无边框）     点边缘线宽
                                     )  # 画出点
                dot_blocks_1 = []
                dot_blocks_2 = []
                for i in range(len(dot_block_1)):
                    if dot_block_2[i] >= dot_block_1[i]:
                        dot_blocks_1.append(dot_block_1[i])
                        dot_blocks_2.append(dot_block_2[i])
                scatter = ax.scatter(dot_blocks_1, dot_blocks_2,  # x轴数据  y轴数据
                                     c='r', s=1.2, alpha=1, marker='o', edgecolors=None, linewidths=0
                                     # 点边缘颜色（None 表示无边框）     点边缘线宽
                                     )  # 画出点
            alignment_ids=set_ks_1(alignment_id,mod)
            set_ks_2(alignment_ids, ax, Ks, mod, ks_mod, NSD)
            if top_block == "False":
                if posd == "bot_left":
                    cbar = plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.1)
                    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）# 画出颜色条
                    cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
                    cbar.ax.set_position([0.145, 0.2, 0.03, 0.2])
                    cbar.ax.yaxis.set_ticks_position('left')
                if posd == "bot_right":
                    cbar = plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.1)
                    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）# 画出颜色条
                    cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
                    cbar.ax.set_position([0.86, 0.2, 0.03, 0.2])
            ax.set_xticks([])  # 隐藏x轴刻度
            ax.set_yticks([])  # 隐藏y轴刻度
            if top_block == "False":
                scatter = ax_4.scatter(reds_1_4, reds_2_4,  # x轴数据  y轴数据
                                       c=dot_cs_4, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
                                       norm=norm, s=1.2,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
                                       alpha=1, marker='o',  # 透明度    点的形状
                                       edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
                                       )  # 画出点
            if top_block == "True":
                scatter = ax_4.scatter(reds_1_4, reds_2_4,  # x轴数据  y轴数据
                                       c='grey', s=1.2, alpha=1, marker='o', edgecolors=None, linewidths=0
                                       # 点边缘颜色（None 表示无边框）     点边缘线宽
                                       )  # 画出点
                dot_blocks_1_4 = []
                dot_blocks_2_4 = []
                for i in range(len(dot_block_1_4)):
                    if dot_block_2_4[i] >= dot_block_1_4[i]:
                        dot_blocks_1_4.append(dot_block_1_4[i])
                        dot_blocks_2_4.append(dot_block_2_4[i])
                scatter = ax_4.scatter(dot_blocks_1_4, dot_blocks_2_4,  # x轴数据  y轴数据
                                       c='r', s=1.2, alpha=1, marker='o', edgecolors=None, linewidths=0
                                       # 点边缘颜色（None 表示无边框）     点边缘线宽
                                       )  # 画出点
            alignment_ids_4=set_ks_1(alignment_id_4,mod)
            set_ks_2(alignment_ids_4, ax_4, Ks, mod, ks_mod, NSD)
            if top_block == "False":
                if posd == "bot_left":
                    cbar = plt.colorbar(scatter, ax=ax_4, fraction=0.03, pad=0.1)
                    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）
                    # 画出颜色条
                    cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
                    cbar.ax.set_position([0.855, 0.6, 0.03, 0.2])
                if posd == "bot_right":
                    cbar = plt.colorbar(scatter, ax=ax_4, fraction=0.03, pad=0.1, location='left')
                    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）
                    # 画出颜色条
                    cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes,ha='center', va='bottom', rotation=0)
                    cbar.ax.set_position([0.12, 0.62, 0.03, 0.2])
                    cbar.ax.yaxis.set_ticks_position('left')
            ax_4.set_xticks([])  # 隐藏x轴刻度
            ax_4.set_yticks([])  # 隐藏y轴刻度
        if half == "True":
            half_t_1(pos, ax)
            if mod == "tpos":
                half_t_2(ax, long1, long1_1, x_lens, y_lens, x_len, y_len, m, n, name_1, name_2)
            elif mod == "torder":
                half_t_2(ax, long1, long1_1, x_lens_2, y_lens_2, x_len_2, y_len_2, m_2, n_2, name_1, name_2)
            dot_cs, reds_1, reds_2 = ks_dot(dot_red_1, dot_red_2, dot_c)
            scatter = ax.scatter(reds_1, reds_2,  # x轴数据  y轴数据
                                 c=dot_cs, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
                                 norm=norm, s=1.2,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
                                 alpha=1, marker='o',  # 透明度    点的形状
                                 edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
                                 )  # 画出点
            alignment_ids=set_ks_1(alignment_id,mod)
            set_ks_2(alignment_ids, ax, Ks, mod, ks_mod, NSD)
            cbar = plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.1)
            cbar.ax.text(0.5, 1.05, 'Ks range', transform=cbar.ax.transAxes, ha='center', va='bottom', rotation=0)
            ax.set_xticks([])  # 隐藏x轴刻度
            ax.set_yticks([])  # 隐藏y轴刻度
        if mod_show == "tree":
            cbar.ax.set_position([0.055, 0.5, 0.03, 0.23])
        if mod_show == "rate":
            cbar.ax.set_position([0.065, 0.47, 0.03, 0.19])
        if mod_show == "tree":
            x_gene, y_gene = get_tree_1(gene_show, gff_1, gff_chr_1, gff_torder_1, chr_torder_1, gff_dict_1, chr_dict_1,gff_2, gff_chr_2, gff_torder_2, chr_torder_2, gff_dict_2, chr_dict_2, mod)
            x_gene = get_tree_2(x_gene)
            y_gene = get_tree_2(y_gene)
            colors, date_colors = tree_co_1(color_date, color_dates)
            tree_y_gene_1(y_gene, mod, n_2, color_date, colors, ax, date_colors, n)
            tree_x_gene_1(x_gene, mod, m_2, color_date, colors, ax, date_colors, m)
        if half == "False":
            ha_f_set(ax, m, n_2, n, m_2)
if half == "self_col_2" or half == "self_blco_2":
    for a in [ax, ax_4]:
        a.set_facecolor('none')
        a.spines[:].set_visible(False)
        a.set_xticks([])
        a.set_yticks([])
    if block == "True":
        colorsd = ['blue', '#ADD8E6', 'yellow', 'orange', 'r']
        cmapd = LinearSegmentedColormap.from_list('custom_gradient', colorsd, N=256)
        boundsd = [0, 10, 20, 30, 40]
        normd = Normalize(vmin=0, vmax=40)  # 线性映射到0-40范围
        cb_ax = inset_axes(
            ax_4,
            width="20%",  # 颜色条宽度（相对主图）
            height="7%",  # 颜色条高度
            loc='upper center',  # 定位在主图下方居中
            bbox_to_anchor=(-0.21, 0.9, 1, 0.1),  # 微调位置（x0, y0, width, height）
            bbox_transform=ax_4.transAxes,
            borderpad=0
        )
        cb = fig.colorbar(
            plt.cm.ScalarMappable(norm=normd, cmap=cmapd),cax=cb_ax,
            orientation='horizontal',ticks=boundsd
        )
        cb.set_ticklabels(['0', '10', '20', '30', '40'])
        cb.ax.text(-0.35, 0.5, 'Block Length', transform=cb.ax.transAxes,ha='center', va='center', fontweight='bold', rotation=0)
    if posd == "bot_left":
        if block == "False":
            if half == "self_blco_2":
                pos_ax = ax.get_position()  # 获取当前的位置盒子
                ax.set_position([pos_ax.x0 + 0.056, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                ax_4.set_position([pos_ax4.x0 - 0.286, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
            else:
                pos_ax = ax.get_position()  # 获取当前的位置盒子
                ax.set_position([pos_ax.x0 + 0.056, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                ax_4.set_position([pos_ax4.x0 - 0.286, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
        if block == "True":
            if collinearity_ks_mod == "True":
                if half == "self_blco_2":
                    if top_block == "False":
                        pos_ax = ax.get_position()  # 获取当前的位置盒子
                        ax.set_position([pos_ax.x0+ 0.053 , pos_ax.y0, pos_ax.width * 1.52, pos_ax.height * 1.72])
                        pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                        ax_4.set_position([pos_ax4.x0 - 0.289, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
                    if top_block == "True":
                        pos_ax = ax.get_position()  # 获取当前的位置盒子
                        ax.set_position([pos_ax.x0 + 0.053, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.73])
                        pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                        ax_4.set_position(
                            [pos_ax4.x0 - 0.289, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
                else:
                    pos_ax = ax.get_position()  # 获取当前的位置盒子
                    ax.set_position([pos_ax.x0 + 0.053, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                    pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                    ax_4.set_position([pos_ax4.x0 - 0.289, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
            if collinearity_ks_mod == "False":
                if half == "self_blco_2":
                    pos_ax = ax.get_position()  # 获取当前的位置盒子
                    ax.set_position([pos_ax.x0 + 0.014, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                    pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                    ax_4.set_position(
                        [pos_ax4.x0 - 0.32, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
                else:
                    pos_ax = ax.get_position()  # 获取当前的位置盒子
                    ax.set_position([pos_ax.x0 + 0.014, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                    pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                    ax_4.set_position([pos_ax4.x0 - 0.32, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
    if posd == "bot_right":
        if block == "False":
            if half == "self_blco_2":
                pos_ax = ax.get_position()  # 获取当前的位置盒子
                ax.set_position([pos_ax.x0 - 0.283, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
            else:
                pos_ax = ax.get_position()  # 获取当前的位置盒子
                ax.set_position([pos_ax.x0 - 0.283, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])

        if block == "True":
            if collinearity_ks_mod == "True":
                if half == "self_blco_2":
                    if top_block == "False":
                        pos_ax = ax.get_position()  # 获取当前的位置盒子
                        ax.set_position([pos_ax.x0 - 0.283, pos_ax.y0, pos_ax.width * 1.52, pos_ax.height * 1.72])
                        pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                        ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
                    if top_block == "True":
                        pos_ax = ax.get_position()  # 获取当前的位置盒子
                        ax.set_position([pos_ax.x0 - 0.33, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.74])
                        pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                        ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
                else:
                    pos_ax = ax.get_position()  # 获取当前的位置盒子
                    ax.set_position([pos_ax.x0 - 0.283, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                    pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                    ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.3, pos_ax4.width * 1.75, pos_ax4.height * 1.75])

            if collinearity_ks_mod == "False":
                if half == "self_blco_2":
                    pos_ax = ax.get_position()  # 获取当前的位置盒子
                    ax.set_position([pos_ax.x0 - 0.34, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                    pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                    ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.29, pos_ax4.width * 1.75, pos_ax4.height * 1.75])
                else:
                    pos_ax = ax.get_position()  # 获取当前的位置盒子
                    ax.set_position([pos_ax.x0 - 0.34, pos_ax.y0, pos_ax.width * 1.75, pos_ax.height * 1.75])
                    pos_ax4 = ax_4.get_position()  # 获取当前的位置盒子
                    ax_4.set_position([pos_ax4.x0, pos_ax4.y0 - 0.29, pos_ax4.width * 1.75, pos_ax4.height * 1.75])

plt.savefig(savefile, dpi=1000, format='png')


