#蛋白质序列  BLOSUM62 15-25个残基 19
#核酸序列   DNA+5/-4  25-50个碱基   29
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import math
import toml
date_norm="False"     #  "False"   "True"  是否归一化   "False" 不进行归一化 "True" 进行归一化
Max_date= "True"     #  "MAX"   "True"  最大值设置   "MAX"是理论完美匹配的得分 "True" 实际最大匹配的得分
A_dna = 25     # 20-30  基准窗口常数
J_dna = 0.22     # 0.15-0.35 谨慎程度
B_dna =-8        # -3 至 -10  补偿强度
C_dna = 200      #150-300    补偿饱和点

A_protein = 15     # 12-18 基准窗口常数
J_protein = 0.28     # 0.15-0.35 谨慎程度
B_protein =-5        # -3 至 -10  补偿强度
C_protein = 100      #50-150 补偿饱和点
figs="False"    # "False"  "True" 图像显示时是否使用真实比例
min_prop="True"     # "False"  "True" 图像显示时是否限制下限
min_props = 0     # min_prop="True" 时，下限为多少
Matrix_dna="BLASTN"     #选择dna矩阵 "BLASTN","TRANS","NUC.4.2"
Matrix_protein="PAM250" #选择蛋白质矩阵 "BLOSUM62" ,"BLOSUM80","BLOSUM30","PAM250","PAM120","PAM30"
show_leng_dna=500       #设置标签的分割
show_leng_protein=100   #设置标签的分割

date_1 = "1.pep"
date_2 = "2.pep"
out_part_date_1 = "1.2pep_seq"
out_part_date_2 = "1.2pep_sco"
out_part_date_3 = "1.2pep_pro"
savefile = "output_pep.png"

current_vars_1=['date_norm','Max_date','A_dna','J_dna','B_dna','C_dna','A_protein','J_protein','B_protein',
              'C_protein','figs','min_prop','min_props','Matrix_dna','Matrix_protein','show_leng_dna','show_leng_protein',]
current_vars_2=['date_1','date_2','out_part_date_1','out_part_date_2','out_part_date_3','savefile']

try:
    # 配置文件路径 - 根据实际情况调整
    config_file_path = '../config.toml'  # 配置文件在上级目录
    config_dir = os.path.dirname(os.path.abspath(config_file_path))

    # 加载TOML配置文件
    with open(config_file_path, "r", encoding="utf-8") as f:
        config = toml.load(f)

    if "Fragment" in config:
        intragenic_config = config["Fragment"]

        # 先处理非路径变量
        for var_name in current_vars_1:
            if var_name in intragenic_config:
                globals()[var_name] = intragenic_config[var_name]

        # 再处理路径变量
        for var_name in current_vars_2:
            if var_name in intragenic_config:
                value = intragenic_config[var_name]  # 直接从配置文件获取
                if isinstance(value, str):
                    if not os.path.isabs(value):  # 相对路径
                        # 转换为相对于配置文件目录的绝对路径
                        value = os.path.join(config_dir, value)
                globals()[var_name] = value

        print("✅ 配置加载成功！")

except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")




dna_tag="bp"
protein_tag="aa"
dna = set('ATCGN')
protein = set('ACDEFGHIKLMNPQRSTVWY')
maxs=0

"""
自适应滑动窗口大小计算器
公式：W = floor( A·[1 + α·ln(max/min)] + B·[min/(C + min)] )
"""

#计算窗口大小
def calc_window(L1, L2, A, α, B, C):
    """
    计算窗口大小
    输入：L1,L2=序列长度, A,α,B,C=公式参数
    输出：窗口大小
    """
    # 计算长序列和短序列
    long = max(L1, L2)
    short = min(L1, L2)
    # 计算长度比值
    ratio = long / short
    # 长度差异调整项
    adjust = A * (1 + α * math.log(ratio)) if ratio > 1.001 else A
    # 短序列补偿项
    compensate = B * (short / (C + short))
    # 合并并取整
    window = int(round(adjust + compensate))
    # 确保在合理范围内（5-100）
    return max(5, min(window, 100))

def no_r_max(Matrix):
    matrix = substitution_matrices.load(Matrix)
    # 将矩阵转换为numpy数组
    size = len(matrix.alphabet)
    matrix_array = np.zeros((size, size))

    for i, a in enumerate(matrix.alphabet):
        for j, b in enumerate(matrix.alphabet):
            matrix_array[i, j] = matrix[a, b]

    return matrix_array.max()

record_1 = next(SeqIO.parse(date_1, "fasta"))
#读取序列
record_1.seq = Seq(str(record_1.seq).replace('-', '').replace('X', '').replace('*', ''))
leng_x=len(record_1.seq)
name_1=record_1.id
#清理序列seq中的特殊字符
record_2 = next(SeqIO.parse(date_2, "fasta"))
#读取序列
record_2.seq = Seq(str(record_2.seq).replace('-', '').replace('X', '').replace('*', ''))
leng_y=len(record_2.seq)
name_2=record_2.id
#清理序列seq中的特殊字符
sequence_1 = set(str(record_1.seq).upper())
sequence_2 = set(str(record_2.seq).upper())
"""
自动判断是 dna protein
"""
if sequence_1.issubset(dna):
    seq_1 = "dna"
elif sequence_1.issubset(protein):
    seq_1 = "protein"
else:
    seq_1 = "unknown"
    print("序列1你输入的文件错误")
if sequence_2.issubset(dna):
    seq_2 = "dna"
elif sequence_2.issubset(protein):
    seq_2 = "protein"
else:
    seq_2 = "unknown"
    print("序列2你输入的文件错误")
if seq_1 == seq_2:
    seqs = seq_1
else:
    seqs = "error"
    print("error,你输入的文件错误")
class SEQ1 :
    def __init__(self,start,end,seq,point):
        self.start = int(start)
        self.end = int(end)
        self.seq = str(seq)
        self.point = int(point)
class SEQ2 :
    def __init__(self,point1,point2,score,start1,start2):
        self.point1 = int(point1)
        self.point2 = int(point2)
        self.score = float(score)
        self.start1 = int(start1)
        self.start2 = int(start2)
alignment_results = []
dots_x=[]
dots_y = []
scores=[]
id_1=[]
id_2=[]
n_seq_1=[]
n_seq_2=[]
if seqs=="dna":
    window_size_dna = calc_window(leng_x,leng_y,A_dna,J_dna,B_dna,C_dna)
    seqs_1 = []
    sequence = record_1.seq
    for start in range(0, len(sequence) - window_size_dna +1):
        end = start + window_size_dna
        point = (end + start + 1)/2
        record_1_seq = sequence[start:end]
        dot=SEQ1(
            start+1,
            end,
            str(record_1_seq),
            point
        )
        id_1.append(int(start+1))
        n_seq_1.append(str(record_1_seq))
        seqs_1.append(dot)
    seqs_2 = []
    sequence = record_2.seq
    for start in range(0, len(sequence) - window_size_dna + 1):
        end = start + window_size_dna
        point = (end + start + 1) / 2
        record_2_seq = sequence[start:end]
        dot = SEQ1(
            start+1,
            end,
            str(record_2_seq),
            point
        )
        id_2.append(int(start + 1))
        n_seq_2.append(str(record_2_seq))
        seqs_2.append(dot)
    matrix = substitution_matrices.load(Matrix_dna)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.mode = 'global'     #设置全局比对
    aligner.open_gap_score = -1000  # 极大的空位开启罚分
    aligner.extend_gap_score = -1000  # 极大的空位延伸罚分
    date_scos = []
    date_pros = []
    maxsd=no_r_max(Matrix_dna) * window_size_dna
    for i, seq1_obj in enumerate(seqs_1):
        date_sco = []
        date_pro = []
        for j, seq2_obj in enumerate(seqs_2):
            score = aligner.score(seq1_obj.seq, seq2_obj.seq)
            alignment_result = SEQ2(
                point1=seq1_obj.point,
                point2=seq2_obj.point,
                score=score,
                start1=seq1_obj.start,
                start2=seq2_obj.start,
            )
            date_sco.append(score)
            date_pro.append(score)
            if maxs < score:
                maxs = score
            alignment_results.append(alignment_result)
            if score >= 0:
                dots_x.append(seq1_obj.point)
                dots_y.append(seq2_obj.point)
                scores.append(score)
        date_scos.append(date_sco)
        date_pros.append(date_pro)
    if Max_date == "True":
        vmax = maxsd
    elif Max_date == "MAX":
        vmax = maxs

    scores_np = np.asarray(scores, dtype=np.float32)
    dots_x_np = np.asarray(dots_x, dtype=np.int32)
    dots_y_np = np.asarray(dots_y, dtype=np.int32)
    if date_norm=="True":
        range_val = vmax
        scores_np = scores_np / range_val
        # 原地归一化（最快）
        if min_prop == "False":
            # === 第2步：快速计算分位数和阈值 ===
            # 方法A：使用quantile（最快最准）
            Q1, Q3 = np.quantile(scores_np, [0.25, 0.75], method='linear')
            lower_bound = Q1 - 1.5 * (Q3 - Q1)
        if min_prop == "True":
            lower_bound = min_props
        # === 第3步：过滤和排序 ===
        # 创建布尔掩码
        mask = scores_np > lower_bound
        # 直接获取过滤后的数据（避免中间变量）
        filtered_scores = scores_np[mask]

        # 使用np.argsort的优化参数
        # kind='quicksort'：默认，对小数据最快
        # kind='stable'：对有重复值的数据稳定
        sorted_idx = np.argsort(filtered_scores, kind='stable')  # 升序

        # 直接构建最终结果（避免多次tolist转换）
        scores = filtered_scores[sorted_idx].tolist()
        dots_x = dots_x_np[mask][sorted_idx].tolist()
        dots_y = dots_y_np[mask][sorted_idx].tolist()
    if date_norm == "False":
        if min_prop == "False":
            lower_bound = 0
        if min_prop == "True":
            lower_bound = min_props
        # === 第3步：过滤和排序 ===
        # 创建布尔掩码
        mask = scores_np > lower_bound
        valid_count = np.sum(mask)
        # 直接获取过滤后的数据（避免中间变量）
        filtered_scores = scores_np[mask]

        # 使用np.argsort的优化参数
        # kind='quicksort'：默认，对小数据最快
        # kind='stable'：对有重复值的数据稳定
        sorted_idx = np.argsort(filtered_scores, kind='stable')  # 升序

        # 直接构建最终结果（避免多次tolist转换）
        scores = filtered_scores[sorted_idx].tolist()
        dots_x = dots_x_np[mask][sorted_idx].tolist()
        dots_y = dots_y_np[mask][sorted_idx].tolist()

if seqs=="protein":
    window_size_protein = calc_window(leng_x, leng_y, A_protein, J_protein, B_protein, C_protein)
    seqs_1 = []
    sequence = record_1.seq
    for start in range(0, len(sequence) - window_size_protein +1):
        end = start + window_size_protein
        point = (end + start + 1)/2
        record_1_seq = sequence[start:end]
        dot=SEQ1(
            start+1,
            end,
            str(record_1_seq),
            point
        )
        id_1.append(int(start + 1))
        n_seq_1.append(str(record_1_seq))
        seqs_1.append(dot)
    seqs_2 = []
    sequence = record_2.seq
    for start in range(0, len(sequence) - window_size_protein + 1):
        end = start + window_size_protein
        point = (end + start + 1) / 2
        record_2_seq = sequence[start:end]
        dot = SEQ1(
            start+1,
            end,
            str(record_2_seq),
            point
        )
        id_2.append(int(start + 1))
        n_seq_2.append(str(record_2_seq))
        seqs_2.append(dot)
    matrix = substitution_matrices.load(Matrix_protein)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.mode = 'global'     #设置全局比对
    aligner.open_gap_score = -1000  # 极大的空位开启罚分
    aligner.extend_gap_score = -1000  # 极大的空位延伸罚分
    date_scos=[]
    date_pros = []
    maxsd=no_r_max(Matrix_protein) * window_size_protein
    for i, seq1_obj in enumerate(seqs_1):
        date_sco = []
        date_pro=[]
        for j, seq2_obj in enumerate(seqs_2):
            score = aligner.score(seq1_obj.seq, seq2_obj.seq)
            alignment_result = SEQ2(
                point1=seq1_obj.point,
                point2=seq2_obj.point,
                score=score,
                start1=seq1_obj.start,
                start2=seq2_obj.start,
            )
            date_sco.append(score)
            date_pro.append(score)
            if maxs < score:
                maxs = score
            alignment_results.append(alignment_result)
            if score >= 0:
                dots_x.append(seq1_obj.point)
                dots_y.append(seq2_obj.point)
                scores.append(score)
        date_scos.append(date_sco)
        date_pros.append(date_pro)

    if Max_date == "True":
        vmax = maxsd
    elif Max_date == "MAX":
        vmax = maxs
    scores_np = np.asarray(scores, dtype=np.float32)
    dots_x_np = np.asarray(dots_x, dtype=np.int32)
    dots_y_np = np.asarray(dots_y, dtype=np.int32)
    if date_norm == "True":
        if Max_date == "True":
            vmax = maxsd
        elif Max_date == "MAX":
            vmax = maxs
        scores_np = np.asarray(scores, dtype=np.float32)
        dots_x_np = np.asarray(dots_x, dtype=np.int32)
        dots_y_np = np.asarray(dots_y, dtype=np.int32)
        # 原地归一化（最快）
        range_val = vmax
        scores_np = scores_np / range_val

        if min_prop == "False":
            # === 第2步：快速计算分位数和阈值 ===
            # 方法A：使用quantile（最快最准）
            Q1, Q3 = np.quantile(scores_np, [0.25, 0.75], method='linear')
            lower_bound = Q1 - 1.5 * (Q3 - Q1)

        if min_prop == "True":
            lower_bound = min_props
        # === 第3步：过滤和排序 ===
        # 创建布尔掩码
        mask = scores_np > lower_bound
        # 直接获取过滤后的数据（避免中间变量）
        filtered_scores = scores_np[mask]

        # 使用np.argsort的优化参数
        # kind='quicksort'：默认，对小数据最快
        # kind='stable'：对有重复值的数据稳定
        sorted_idx = np.argsort(filtered_scores, kind='stable')  # 升序

        # 直接构建最终结果（避免多次tolist转换）
        scores = filtered_scores[sorted_idx].tolist()
        dots_x = dots_x_np[mask][sorted_idx].tolist()
        dots_y = dots_y_np[mask][sorted_idx].tolist()
    if date_norm == "False":
        if min_prop == "False":
            lower_bound = 0
        if min_prop == "True":
            lower_bound = min_props
        # === 第3步：过滤和排序 ===
        # 创建布尔掩码
        mask = scores_np > lower_bound
        valid_count = np.sum(mask)
        # 直接获取过滤后的数据（避免中间变量）
        filtered_scores = scores_np[mask]

        # 使用np.argsort的优化参数
        # kind='quicksort'：默认，对小数据最快
        # kind='stable'：对有重复值的数据稳定
        sorted_idx = np.argsort(filtered_scores, kind='stable')  # 升序

        # 直接构建最终结果（避免多次tolist转换）
        scores = filtered_scores[sorted_idx].tolist()
        dots_x = dots_x_np[mask][sorted_idx].tolist()
        dots_y = dots_y_np[mask][sorted_idx].tolist()

if seqs=="protein" or seqs=="dna":
    if seqs=="protein":
        show_leng=show_leng_protein
        tag=protein_tag
    else:
        show_leng = show_leng_dna
        tag = dna_tag

    maxs_length = max(len(id_1), len(n_seq_1), len(id_2), len(n_seq_2))

    # 填充到相同长度（用 None 或空字符串填充）
    name_1_filled = id_1 + [None] * (maxs_length - len(id_1))
    n_seq_1_filled = n_seq_1 + [None] * (maxs_length - len(n_seq_1))
    name_2_filled = id_2 + [None] * (maxs_length - len(id_2))
    n_seq_2_filled = n_seq_2 + [None] * (maxs_length - len(n_seq_2))

    df_1 = pd.DataFrame({
        'id_1': name_1_filled,
        'seq_1': n_seq_1_filled,
        'id_2': name_2_filled,
        'seq_2': n_seq_2_filled
    })
    df_2 = pd.DataFrame(date_scos).T  # 转置
    df_2.columns = id_1
    df_2.index = id_2
    if date_norm == "True":
        df_3 = pd.DataFrame(date_pros).T  # 转置
        df_3.columns = id_1
        df_3.index = id_2

    if figs == "True":
        base_size = 8
        fig_width = base_size
        fig_height = base_size * (leng_y / leng_x)
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))  # 创建图形和坐标轴
    ax.spines['top'].set_visible(False)  # 隐藏上边框
    ax.spines['right'].set_visible(False)  # 隐藏右边框
    ax.spines['left'].set_visible(False)  # 隐藏左边框
    ax.spines['bottom'].set_visible(False)  # 隐藏下边框
    ax.set_xticks([])  # 隐藏x轴刻度
    ax.set_yticks([])  # 隐藏y轴刻度

    ax.vlines(leng_x, 0, leng_y, lw=1, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    ax.vlines(0, 0, leng_y, lw=1, color='k', linestyle='-', alpha=0.5)  # 画垂直线
    ax.hlines(0, 0, leng_x, lw=1, color='k', linestyle='-', alpha=0.5)  # 画水平线
    ax.hlines(leng_y, 0, leng_x, lw=1, color='k', linestyle='-', alpha=0.5)  # 画水平线
    ax.text(leng_x/2,1.08*leng_y,fr"$\boldsymbol{{{name_1}}}$ ({leng_x}{tag})", fontsize=12, ha='center', va='top')
    ax.text(-0.15 * leng_x, leng_y / 2, fr"$\boldsymbol{{{name_2}}}$ ({leng_y}{tag})", fontsize=12, ha='right', va='center', rotation=90)
    colors = ['#ffffff', '#ffcccc','#ff9999',  '#ff0000']
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
    if date_norm == "True":
        norm = plt.Normalize(vmin=0, vmax=1)
    else:
        norm = plt.Normalize(vmin=0, vmax=vmax)
    ax.invert_yaxis()
    scatter = ax.scatter(dots_x, dots_y,  # x轴数据  y轴数据
                         c=scores, cmap=cmap,  # 颜色数据（可以是数值或颜色列表）  颜色映射（colormap）
                         norm=norm, s=1.2,  # 数据归一化（Normalize 对象，控制颜色映射范围）   点的大小（size）
                         alpha=1, marker='o',  # 透明度    点的形状
                         edgecolors=None, linewidths=0  # 点边缘颜色（None 表示无边框）     点边缘线宽
                         )  # 画出点
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.1, aspect=15)
    # fraction：颜色条相对于绘图区域的高度/宽度比例（默认0.15）# pad：颜色条与主图的间距（单位是图形宽度/高度的比例）
    # 画出颜色条
    if date_norm == "True":
        cbar.ax.text(0.5, 1.05, 'PCT', transform=cbar.ax.transAxes, ha='center', va='bottom', rotation=0)
    else:
        cbar.ax.text(0.5, 1.05, 'score', transform=cbar.ax.transAxes, ha='center', va='bottom', rotation=0)
    cbar.ax.set_position([0.855, 0.65, 0.1, 0.15])
    if min_prop == "True":
        target_value = min_props  # 要标记的数值
        y_pos = (target_value - norm.vmin) / (norm.vmax - norm.vmin)

        # 添加箭头和文本
        cbar.ax.annotate('',
                         xy=(0.2, y_pos),  # 箭头终点（颜色条中心）
                         xytext=(-2.5, y_pos),  # 箭头起点（颜色条右侧）
                         xycoords='axes fraction',
                         arrowprops=dict(arrowstyle='->', color='black', lw=0.8),
                         annotation_clip=False)

        # 添加数值标签
        cbar.ax.text(-3, y_pos, f'{target_value}',
                     transform=cbar.ax.transAxes,
                     ha='right', va='center',
                     color='red', fontsize=10,
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    cbar.ax.set_aspect('auto')
    ax3 = ax.twiny()
    ax3.set_xlim(ax.get_xlim())  # 确保与ax的x轴范围一致
    ax3.spines['top'].set_position(('data', 0))
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    # 设置刻度位置在顶部
    ax3.xaxis.set_ticks_position('top')
    ax3.xaxis.set_label_position('top')
    ax3_ticks = np.arange(0, leng_x, show_leng)
    ax3.set_xticks(ax3_ticks)
    ax3.set_xticklabels(ax3_ticks)

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.spines['left'].set_position(('data', 0))
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    # 设置刻度位置在顶部
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2_ticks = np.arange(0, leng_y, show_leng)
    ax2.set_yticks(ax2_ticks)
    ax2.set_yticklabels(ax2_ticks)
    #ax2.invert_yaxis()     # 反转y轴方向
    df_1.to_csv(out_part_date_1, sep='\t', index=False, encoding='utf-8')
    df_2.to_csv(out_part_date_2, sep='\t', index=True, encoding='utf-8')
    if date_norm == "True":
        df_3.to_csv(out_part_date_3, sep='\t', index=True, encoding='utf-8')
    plt.savefig(savefile, dpi=1000, format='png', bbox_inches='tight')
    plt.close()
    # plt.show()