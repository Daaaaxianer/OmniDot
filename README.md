# OmniDot
A multi-scale dot plot visualization tool for comparative genomics. Analyze alignments from genomic fragments and chromosomes down to intragenic micro-structures.

**OmniDot**
 is a comprehensive dot plot visualization tool designed for multi-scale comparative genomics.

Traditional tools often struggle to bridge the gap between whole-genome synteny and sequence-level details. OmniDot is built to provide a unified visualization experience across all genomic dimensions:

* **Chromosomal Level**
: Visualize macro-synteny and large-scale structural variations between genomes.
* **Fragment Level**: efficiently handle and visualize unplaced **genomic fragments (scaffolds/contigs)**
 to assist in assembly validation.
* **Intragenic Level**: Zoom into the micro-scale to reveal **intragenic**
 alignments, exon/intron boundaries, and local repeats.

From broad evolutionary patterns to precise gene-level architecture, OmniDot helps you see the whole picture.



本工具由3个模块组成，基于**Linux 操作系统**编写的工具

## 系统要求

### 软件环境

- **操作系统**：Linux（推荐 Ubuntu 20.04+）、macOS 或 Windows（需 WSL2）
- **Python**：3.8 或更高版本
- **必要工具**：minimap2（已包含在 `./../bin/` 目录中）

## 核心模块介绍

### 1. Chromosomal Level

**功能**：生成高质量的染色体比对可视化图，支持多种分析模式。核心原理是 **基因组比对数据的坐标映射和可视化**。

- **BLAST模式**：显示基因间的一对一、一对多比对关系

- **共线性模式**：显示基因组间的共线性区块
- **组合模式**：可同时显示BLAST和共线性结果

#### 主要工作流程

1. **数据读取与处理**
   - 读取GFF、染色体长度、比对结果等文件
   - 构建基因位置映射字典
   - 过滤和预处理比对数据
2. **坐标计算**
   - 根据选择的坐标系统（torder或tpos）计算基因位置
   - 处理染色体边界和标签位置
3. **图形生成**
   - 创建主比对图
   - 添加染色体边界、标签
   - 绘制比对点或共线性区块
   - 添加颜色条、图例等辅助元素
4. **布局调整**
   - 调整子图位置和大小
   - 处理重叠和间距
   - 保存最终图像

#### 原理

#### 1. **数据模型与坐标系统**

脚本支持两种基因定位方式：

##### a. **torder（基因顺序）坐标系统**

- **原理**：按基因在染色体上的顺序编号进行定位
- **计算方法**：基因在总坐标系中的位置 = 前序染色体上的基因总数 + 该基因在当前染色体上的顺序号
- **优点**：均匀分布，不受染色体大小影响
- **适用**：基因密度比较、共线性模式分析

##### b. **tpos（物理位置）坐标系统**

- **原理**：按基因在染色体上的实际物理位置（bp）进行定位
- **计算方法**：基因在总坐标系中的位置 = 前序染色体的总长度 + 该基因在当前染色体上的中点位置

- **优点**：反映真实物理距离
- **适用**：物理距离相关的分析、密度计算

#### 2. **染色体坐标系构建**

把所有染色体"首尾相接"地排列在一条直线上，构建了一个连续的坐标空间。

1号染色体的结尾下一位是2号染色体的开头。

```
n = 0  # 累积位置
for 每个染色体 in 物种1:
    染色体起始位置[染色体] = n
    n = n + 染色体长度[染色体]
    染色体标签位置[染色体] = 起始位置 + 染色体长度/2  # 用于显示染色体编号
```

#### 3. **数据映射原理**

```
query基因位置 = query基因所在的染色体起始位置 + query基因在染色体上的位置
subject基因位置 = subject基因所在的染色体起始位置 + subject基因在染色体上的位置
```

每个比对关系在图中表现为一个点：(query位置, subject位置)

#### 4.共线性区块的识别

1. **区块检测**：基于连续的基因对构建共线性区块
2. **区块合并**：根据距离阈值合并相邻的共线性区块
3. **Ks值整合**：计算区块内基因对的Ks统计值（平均值/中位数）

#### 5. **关键算法原理**

##### a. **共线性区块合并算法**

过滤小区块

合并两个距离相近的区块

##### b**对称性与半数显示原理**

当分析同一物种的自身比对时（如全基因组复制分析），比对图具有对称性。利用这一特性：在图上显示2个物种自身比对比对图

### 2. Fragment Level

**功能**：生成序列间的相似性点图，用于直观显示序列间的相似区域。

**点阵图（Dot Plot）**生成工具。点阵图是生物信息学中用于可视化两个序列之间相似性的经典方法。

#### 主要功能

#### **序列比对分析**

#### **自适应窗口大小**

$$
W=⌊A⋅[1+α⋅ln( 
max(L 
1
​
 ,L 
2
​
 )/
 min(L 
1
​
 ,L 
2
​
 )
​
 )]+B⋅( 
min(L 
1
​
 ,L 
2
​
 )/
 (C+min(L 
1
​
 ,L 
2
​
 ))
​
 )⌋
$$

##### **输入变量**

- `L1`, `L2`: 两个输入序列的长度
- `max`: 较长序列的长度 `max(L1, L2)`
- `min`: 较短序列的长度 `min(L1, L2)`
- `ratio = max/min`: 序列长度比值

##### **参数说明**

| 参数      | 名称         | 作用                 | 典型范围                    |
| :-------- | :----------- | :------------------- | :-------------------------- |
| **A**     | 基准窗口常数 | 确定基本窗口大小     | DNA: 20-30 蛋白质: 12-18    |
| **α (J)** | 谨慎程度     | 控制长度差异的敏感性 | 0.15-0.35                   |
| **B**     | 补偿强度     | 调整短序列的补偿量   | -3 至 -10                   |
| **C**     | 补偿饱和点   | 控制补偿的饱和程度   | DNA: 150-300 蛋白质: 50-150 |

##### **第一部分：长度差异调整项**

```
A·[1 + α·ln(ratio)]
```

- **功能**：根据序列长度差异调整窗口大小
- **工作原理**：
  - 当`ratio=1`（序列等长）：`ln(1)=0`，该项简化为`A`
  - 当`ratio>1`（长度不同）：`ln(ratio)>0`，窗口增大
  - `α`控制对长度差异的敏感度

##### **第二部分：短序列补偿项**

```
B·[min/(C + min)]
```

- **功能**：对短序列进行补偿调整
- **工作原理**：
  - 当`min`很小：`min/(C+min) ≈ min/C`，补偿较小
  - 当`min`很大：`min/(C+min) ≈ 1`，达到饱和补偿
  - `B`为负值，实际是减少窗口大小
  - `C`决定补偿的饱和点

####  **公式的生物学意义**

| 序列情况       | 窗口大小趋势 | 原因                             |
| :------------- | :----------- | :------------------------------- |
| **序列等长**   | 接近基准值A  | 不需要特殊调整                   |
| **长度差异大** | 适当增大     | 适应不同尺度的比较               |
| **序列很短**   | 相对较大     | 避免窗口太小失去统计意义         |
| **序列很长**   | 相对较小     | 提高分辨率，避免过大窗口模糊细节 |

#### **数据处理选项**

- **归一化**：可选择是否对得分进行归一化
- **最大值设置**：可使用理论最大值或实际最大值
- **过滤阈值**：可设置显示的下限值
- **值过滤**：基于0自动过滤低信息值

####  **可视化输出**

- 生成高质量点阵图

#### **数据输出**

- 生成三个数据文件：
  - 序列窗口信息（CSV格式）
  - 原始得分矩阵（CSV格式）
  - 归一化后得分矩阵（CSV格式，可选）

#### 核心算法流程

1. **序列预处理**：清理特殊字符，识别序列类型
2. **窗口计算**：根据序列长度自适应计算窗口大小
3. **序列分割**：将序列分割为滑动窗口
4. **比对打分**：使用指定矩阵进行全局比对
5. **数据过滤**：根据阈值过滤低分值点
6. **可视化**：生成点阵图并保存结果

### Intragenic Level

本项目提供了一个完整的基因组片段自比对分析流程，通过4个独立的Python脚本实现了从基因组序列处理到可视化展示的全过程。

#### 文件结构说明

#### 主要脚本文件

###### 1.`make_date_1.py` - 基因组片段切割模块

- **功能**：将大基因组文件按指定参数切割成小片段
  - 输入：完整的染色体FASTA文件
  - 输出：按chunk_size（默认2000bp）切割的片段FASTA文件
  - 关键参数：
    - `new_chr`：目标染色体编号
    - `start/end`：切割的起始和结束位置（可设置为"start"/"end"使用整个染色体）
    - `chunk_size`：每个片段的长度
  - 输出序列ID格式：`chr1:20000000:20002000`（包含位置信息）

######  2.`make_minimap_2.py` - 序列比对模块

**功能**：使用minimap2进行自比对

- 调用minimap2工具将切割后的片段进行自比对
- 参数设置：使用`-ax ava-ont`模式，适合长序列比对
- 输出：SAM格式的比对结果文件

###### 3.`make_date_bed_3.py` - 比对结果解析模块

**功能**：解析SAM文件，提取比对统计信息

- 使用pysam库读取SAM文件
- 解析CS标签获取详细的比对统计：
  - 匹配数、错配数、插入/缺失事件
  - 计算三种同一性百分比（基于匹配、基于事件、总体）
- 对称化处理：使比对结果双向对称
- 输出：TSV格式的比对统计表格

###### 4.`make_fignew_4.py` - 可视化模块

**功能**：生成点图可视化比对结果

- 读取比对统计表格

- 创建带颜色映射的散点图：
  - x轴：比对片段的中点位置
  - y轴：比对片段的距离差
  - 颜色：基于同一性百分比（byevents）

- 添加同一性分布直方图

- 支持四种布局（top/bot/left/right）

- 输出：高分辨率PNG图片

#### 整体流程

完整基因组FASTA
    ↓ (make_date_1.py)
切割成小片段FASTA
    ↓ (make_minimap_2.py)
minimap2自比对 → SAM文件
    ↓ (make_date_bed_3.py)
解析统计信息 → TSV表格
    ↓ (make_fignew_4.py)
生成可视化点图 → PNG图像

## 使用

### 要求

- **操作系统**：Linux（推荐 Ubuntu 20.04+）、macOS 或 Windows（需 WSL2）
- **Python**：3.8 或更高版本
- **必要工具**：minimap2（已包含在 `./../bin/` 目录中）

使用Chromosomal Level

需要安装matplotlib，pandas，numpy，toml，statistics，math，os库

使用 Fragment Level

需要安装Bio，pandas，numpy，toml，math，matplotlib，os库

使用Intragenic Level

需要安装os，pandas，toml，Bio，sys，re，pysam，subprocess，pathlib，matplotlib，numpy库

全部使用需要安装os，pandas，toml，Bio，sys，re，pysam，subprocess，pathlib，matplotlib，numpy，statistics，math库

### 输入文件格式

#### Chromosomal 

##### GFF格式

| 列   | 信息      | 解释                        |
| ---- | --------- | --------------------------- |
| 1    | Chr       | 染色体                      |
| 2    | ID        | 基因名称                    |
| 3    | Start     | 基因的起始位置              |
| 4    | End       | 基因的末端位置              |
| 5    | Direction | 基因序列的方向              |
| 6    | Order     | 每条染色体的顺序，从 1 开始 |
| 7    | told_name | 原始 ID                     |

##### lens

| 列   | 信息   | 解释           |
| ---- | ------ | -------------- |
| 1    | Chr    | 染色体         |
| 2    | Length | 染色体序列长度 |
| 3    | Number | 染色体基因数量 |

##### Collinearity

```python
# Alignment 1: score=1451 pvalue=0.2214 N=49 1&1 plus
vivi01g00514 514 vivi01g00515 515 1
```

##### KS

```python
id1 id2    ka_NG86    ks_NG86    ka_YN00    ks_YN00
vivi01g00514    vivi01g00515   0.0686 0.1182 0.0678 0.1258
vivi01g00517    vivi01g00516   0.0796 0.137  0.0835 0.1197
vivi01g00518    vivi01g00517   0.0278 0.0306 0.0299 0.0248
```

##### color_dates

| 列   | 信息  | 解释           |
| ---- | ----- | -------------- |
| 1    | ID    | 基因名称       |
| 2    | color | 基因显示的颜色 |

##### gene_show

| 列   | 信息 | 解释     |
| ---- | ---- | -------- |
| 1    | Chr  | 染色体   |
| 2    | ID   | 基因名称 |

#### Fragment 

2个标准的只有一个序列的fast文件

```
>schi1g0000008
MSKTGYYLDLCVFLVTLILIITTRTAAKDPSTILSRFQQYLQINTAQPHPNYYEAAEFII
SQAKLLSLESQTLEFVKGKPLILLKWPGKDPTLPSILLNSHTDVVPSEHHKWTHPPFSAH
LDSTTGNIFARGSQDMKCVGLQYLEAIRKLKSYGFRPLRTLYLSFLPDEEIGGNDGARKF
VDSDVFAKMNVGIVLDEGLASPTDNYRAFYGERSPWWLVVKAVGAPGHGAKLYDNTAMEN
LLKSIEIIRRFRAAQFDLVKAGQKAEGEVISVNMVFLKAGTPSPSGFVMNLQPSEAQAGF
DIRVPPTADQASLERLIADEWAPASRNMTFEFKQKVSVNDKLGRPAVTAVDSSNIWWALF
EEAIIKANARLGKPEIFPASTDARYFRERGLPAIGFSPMANTPILLHDHNEFLNKDEYLK
GIDVYESIIKTYASYIQYRRDDASREELKVSLSVNQYFMYDYATHTLVSVVDVSIPYCSN
NQTQFPDEENAKAISLFEDVVHLNDAIEDEAKRMENLVKGIFAGNLFDLGSAQLGNVSKM
THI*
```

```
>schi10g0000796
MESSGELVPFPLLTTPIESNYRACTIPYRFPSDNPKKPTPTELSWIDLFMNSIPSFRKRA
ESDDSVPDAPIRAEKFAQRYSAILEDMKKDPESHGGPPDCILLCRLREQVLREVGFRDIF
KKVKDEENAKAISLFKDVVSLNDAIEDEAKRVENLVRGIFAGNIFDLGSAKLAELFSEDG
ISFLASFQNLVPRPWVIDDLDVFITKWSKKTWKKAVIFVDNSGADVILGILPFARELLRH
GAQVVLAANDLPSINDVTYPELVEIISKLKDEHGKLIGVDTSNLLVANSGNDLPVIDLTT
VSQELAYLASDADLVIVEGMGRGIETNLYARFKCDSLKIGMVKHPEVAQFLGGRLYDCVF
KFNEASS*
```

#### Intragenic 

1个标准的只有一个序列的cds fast文件

```
>chr1
CAGTAAAGTTTGCAAAGAGATTCTGGCAAAGTTGTAATCATCGAATGTGAAATGTCAAGT
AAGTGTAAATGGATTAGTTTGCCAATGGAGGTCGGCAACTCTTTCACATCTGCTCCAAAC
AATTTTAGGACATGCAAGTATTTGAACTTTGATAACATATCATCAGCTATGCCACCCTCC
AGAAATAATGTGCGAAGTGATGCTGATTTTTTTTCATTGATGATTTCCACCATTCTTTCG
GATGAGTATACTGCAAGGTAGCGGTCCTCATTGCTGCTGCTACGATTCAAAATTGATTTT
GCAAAATCGTGCACAAGGTCATGCATTTTATACCTATTTTCTTTTACTTCTTCCAATAAG
GAAGTTTGCAGCAAAATTCTCAAATACTCACATCCTATTTTCTCCATCGTTCTTTCATTT
TGGGAATCCGGTTGAAGAAAGCCTTCAGCCATCCAAAGCTCAACTAGTAGATCTTGTTCC
AACTCAGCATCTTGATCAAAAATTGAGCAATATGCAAAACATTTCTTAACCGGTGCAAGT
GACAGATGATCAAAACTCACCTTAATTATTTGCTCGATCCGACCTTGGTCTCCATTCAAC
AGACTCTCCTCCAAAATAGATTGCCACTCCTCTTTTCTCTTTTTATATAACAAACCTCCA
ATTAACTTTGCTGCCAGAGGTAGACCGTCACATCTTCTTAAAACTCGCTCCCTTATGTCT
TCCAATTCTTTTGGTACTTCTTCCCCTACAGTTGCCCATTTATTAATGATGGACCAGCAA
TCATCATTGCATAGCTTTCCTAGCCCATGGCGAGTAAAATTGATCTGCGGAAGACCGGAC
AGAGTAGTTTCCACTTCTTGCAGACGAGTGGTAACAAGACACCAGCTCCCTTTCTTCGCT
TCGAGTGCCTTCAATGTGGTGAAAAAGTCATCCAATTGTTGATGATTCCACAAATCGTCA
AGAACAAGCAAACATCTTTTTCCCTTGATTTCATTTTGAATTCCTCGAACTATTACGTCC
CTAACATCCGCGTCAGGCTTGTTTCCTGTTGACGATTCTAGAATCATTTTGAAGAGATCC
ATGGTTTCAACTTCTTTAGCCACACAAACCCAAATTTTTTTGTCAAAATGATTATCATTA
AACTGTGGATTGTTGTAAACGGCTTTAG
```

开头必须为     >chr

如  >chr1  >chr2  >chr3  >chr4

### 使用方法

配置文件

#### 第一步

本工具运行需要运行all_run.py
首先配置一个与all_run.py在同一目录下的config.toml文件

config.toml

```
[Config]
config_1 = "Fragment"    
config_2 = "blast"          
```

config_1 = "Fragment" 

3种模式输入 "Chromosomal" "Fragment" "Intragenic"以调用使用的模块

config_2 = "blast" 

当config_1 ="Chromosomal" 时需要写，代表使用的是Chromosomal情况下的哪一个小板块，如"blast","Collinearity","blast_blast","Collinearity_Collinearity","blast_Collinearity"

#### 第二步

在 config_1 的对应模块文件中
"Chromosomal" : Chromosomal Level文件夹中配置一个config.toml文件

"Fragment "  :  Fragment Level文件夹中配置一个config.toml文件

"Intragenic "  : Intragenic Level文件夹中的project文件夹中的src文件夹中配置一个config.toml文件

##### Chromosomal 配置文件格式

###### config.toml文件

使用表分割不同的功能

表是键值对的集合，用方括号表示

配置文件中所有的文件地址都用单引号 ' ' 包括

```
[blast]
half = "False"
mods = "blast"
blast_name = 'D:\Python_cod\vivi_coca.1e-5.blast'
...
[Collinearity]
half = "False"
mods = "collinearity"
Collinearity = 'D:\Python_cod\vivi_coca.collinearity_new.txt'
...
```

###### [blast]

画blast时的标准格式

```toml
[blast]
half = "False"
mods = "blast"
mod = "torder"  
blast_name = 'D:\Python_cod\vivi_coca.1e-5.blast'  
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens'  
name_1 = "V.vinifera" 
name_2 = "C.canephora" 
mod_show = "nan"    
bin_sizes = 20 
wide = 400  
longs = 350  
color_date = "False"  
gene_show = 'D:\Python_cod\vivi_coca_gene.gff'  
color_dates = 'D:\Python_cod\vivi_coca_c.gff'  
savefile= 'graph39_11.png'   
```

解释：

```toml
half = "False"
```

默认参数，不变或不写

```toml
mods = "blast"
```

默认参数，不变

```toml
mod = "torder"
```

"torder"或"tpos"        

根据gff文件的 tpos 或 torder 画图 "torder" 根据基因在染色体总的基因上的相对顺序 "tpos" 根据基因在染色体上的真实位置

```toml
blast_name = 'D:\Python_cod\vivi_coca.1e-5.blast'  
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens
```

写入blast文件的路径

写入blast第1列的gff文件的路径

写入blast第2列的gff文件的路径

写入blast第1列的lens文件的路径

写入blast第2列的lens文件的路径

```toml
name_1 = "V.vinifera" 
name_2 = "C.canephora" 
```

写入blast第1列物种名称

写入blast第2列物种名称

```toml
mod_show = "nan"    
```

"tree" ,"density" ,"nan"    

 "nan" 都不显示    "tree" 显示特定基因的名称  "density"显示基因的密度

```toml
bin_sizes = 20  
```

mod_show ="density"下定义分箱大小  默认单位为百万 值越小细节越丰富

bin_sizes = 20  分箱大小为 20M

```toml
wide = 400  
longs = 350  
```

mod_show ="tree"情况下控制显示的名称上下的距离差

mod_show ="tree"情况下控制显示的名称左右的距离差 

```toml
color_date = "False"  
```

"False" "True"

mod_show ="tree"情况下是否有显示基因名称的连接线的颜色文件 

```toml
gene_show = 'D:\Python_cod\vivi_coca_gene.gff'  
color_dates = 'D:\Python_cod\vivi_coca_c.gff'  
```

写入mod_show ="tree"情况下显示的基因名称的文件路径

写入mod_show ="tree"且color_date = "True"情况下显示的基因名称对应颜色的文件路径  

```toml
savefile= 'graph39_11.png'    
```

写入最后产生的图片路径

###### [Collinearity]

画Collinearity时的标准格式

```toml
[Collinearity]
half = "False"
mods = "collinearity"
mod = "torder"  
Collinearity = 'D:\Python_cod\vivi_coca.collinearity_new.txt'   
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens' 
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens'  
name_1 = "V.vinifera"  
name_2 = "C.canephora"  
mod_show = "nan"  
bin_sizes = 20  
name_3 = "one Vivi : Coca" 
name_4 = "one Coca : Vivi"
gene_show = 'D:\Python_cod\vivi_coca_gene.gff'  
wide = 400  
longs = 350  
color_date = "False"  
color_dates = 'D:\Python_cod\vivi_coca_c.gff'  
collinearity_ks_mod = "True"  
Collinearity_ks = 'D:\Python_cod\vivi_coca.collinearity.ks.txt' 
Ks = "True"  
ks_mod = "Median"  
name_id = "False"  
N = 5  
NSD = 10  
Distance = 250  
savefile= 'graph39_11.png'   
```

解释：

```toml
half = "False"
```

默认参数，不变

```toml
mods = "collinearity"
```

默认参数，不变

```toml
mod = "torder"
```

"torder"或"tpos"        

根据gff文件的 tpos 或 torder 画图 "torder" 根据基因在染色体总的基因上的相对顺序 "tpos" 根据基因在染色体上的真实位置

```toml
Collinearity = 'D:\Python_cod\vivi_coca.collinearity_new.txt'
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens
```

写入collinearity文件路径

写入collinearity第1列的gff文件的路径

写入collinearity第2列的gff文件的路径

写入collinearity第1列的lens文件的路径

写入collinearity第2列的lens文件的路径

```toml
name_1 = "V.vinifera" 
name_2 = "C.canephora" 
```

写入blast第1列物种名称

写入blast第2列物种名称

```
mod_show = "nan" 
```

"rate","tree","density","nan"

"rate"显示基因的比例 "tree"显示特定基因的名称 "density"显示基因的密度 "nan"都不显示

```toml
bin_sizes = 20  
```

mod_show ="density"下定义分箱大小  默认单位为百万 值越小细节越丰富

bin_sizes = 20  分箱大小为 20M

```toml
wide = 400  
longs = 350  
```

mod_show ="tree"情况下控制显示的名称上下的距离差

mod_show ="tree"情况下控制显示的名称左右的距离差 

```toml
color_date = "False"  
```

"False" "True"

mod_show ="tree"情况下是否有显示基因名称的连接线的颜色文件 

```toml
gene_show = 'D:\Python_cod\vivi_coca_gene.gff'  
color_dates = 'D:\Python_cod\vivi_coca_c.gff'  
```

写入mod_show ="tree"情况下显示的基因名称的文件路径

写入mod_show ="tree"且color_date = "True"情况下显示的基因名称对应颜色的文件路径  

```toml
name_3 = "one Vivi : Coca" 
name_4 = "one Coca : Vivi" 
```

在mod_show="rate"下需要填写显示基因的比例的名称 第1列：第2列

在mod_show="rate"下需要填写显示基因的比例的名称 第2列：第1列

```toml
collinearity_ks_mod = "True"  
```

"False" "True"   

是否输入ks文件 

```toml
Collinearity_ks = 'D:\Python_cod\vivi_coca.collinearity.ks.txt'
```

collinearity_ks_mod="True"情况下 写入ks文件路径

```toml
Ks = "True"  
```

"False" "True"

collinearity_ks_mod="True"情况下是否显示ks的值

```toml
ks_mod = "Median"
```

 "Average" 或 "Median"  

collinearity_ks_mod="True"情况下 "Average"显示ks的平均数"Median"显示ks的中位数

```toml
name_id = "False"  
```

"True", "False" 

是否去除物种自身比较时的对角线 (当并非是物种自身比较时 应为 "False" )

```toml
N = 5 
NSD = 10  
Distance = 250 
```

最小对齐数阈值 ≥N 个对齐的区块才合并

显示阈值 ≥NSD 个对齐的区块才显示

区块合并距离 相对顺序差＜250 的区块合并

```toml
savefile= 'graph39_11.png'    
```

写入最后产生的图片路径

###### [blast_blast]

将2个自身与自身blast比对图以三角形拼接到一起

画blast_blast时的标准格式

```toml
[blast_blast]
half = "self_col_2"
mods = "blast"
mod = "torder" 
blast_name = 'D:\Python_cod\vivi_coca.1e-5.blast'  
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens'  
blast_name_2= 'D:\Python_cod\coca_coca.1e-5.blast'  
blast_gff_1_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1_2 = 'D:\Python_cod\blast文件\coca.lens' 
blast_lens_2_2 = 'D:\Python_cod\blast文件\coca.lens' 
name_1 = "V.vinifera"  
name_2 = "C.canephora" 
name_1_4 = "B.hispida" 
name_2_4 = "B.hispida"  
posd = "bot_left"  
savefile= 'graph39_11.png'   
```

解释：

```toml
half = "self_col_2"
```

默认参数，不变

```toml
mods = "blast"
```

默认参数，不变

```toml
mod = "torder"
```

"torder"或"tpos"        

根据gff文件的 tpos 或 torder 画图 "torder" 根据基因在染色体总的基因上的相对顺序 "tpos" 根据基因在染色体上的真实位置

```toml
blast_name = 'D:\Python_cod\vivi_coca.1e-5.blast'  
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens
blast_name_2= 'D:\Python_cod\coca_coca.1e-5.blast'  
blast_gff_1_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1_2 = 'D:\Python_cod\blast文件\coca.lens' 
blast_lens_2_2 = 'D:\Python_cod\blast文件\coca.lens' 
```

写入blast文件的路径

写入blast第1列的gff文件的路径

写入blast第2列的gff文件的路径

写入blast第1列的lens文件的路径

写入blast第2列的lens文件的路径

写入第2个blast文件的路径

写入第2个blast第1列的gff文件的路径

写入第2个blast第2列的gff文件的路径

写入第2个blast第1列的lens文件的路径

写入第2个blast第2列的lens文件的路径

```toml
name_1 = "V.vinifera" 
name_2 = "C.canephora" 
name_1_4 = "B.hispida"  
name_2_4 = "B.hispida" 
```

写入blast第1列物种名称

写入blast第2列物种名称

写入第2个blast第1列物种名称

写入第2个blast第2列物种名称

```toml
posd = "bot_left" 
```

"bot_left"   "bot_right"

图像的格式，第一幅图在左上方还是左下方

```toml
savefile= 'graph39_11.png'    
```

写入最后产生的图片路径

###### [Collinearity_Collinearity]

将2个自身与自身Collinearity比对图以三角形拼接到一起

画Collinearity_Collinearity时的标准格式

```toml
[Collinearity_Collinearity]
half = "self_col_2"
mods = "collinearity"
mod = "torder"
Collinearity = 'D:\Python_cod\vivi_coca.collinearity_new.txt'   
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens'  
Collinearity_2='D:\Python_cod\coca_coca.collinearity.txt'  
blast_gff_1_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1_2 = 'D:\Python_cod\blast文件\coca.lens' 
blast_lens_2_2 = 'D:\Python_cod\blast文件\coca.lens' 
name_1 = "V.vinifera"  
name_2 = "C.canephora"  
name_1_4 = "B.hispida" 
name_2_4 = "B.hispida"  
posd = "bot_left"  
collinearity_ks_mod = "True"  
Collinearity_ks = 'D:\Python_cod\vivi_coca.collinearity.ks.txt'  
Collinearity_ks_2= 'D:\Python_cod\coca_coca.collinearity.ks.txt'  
Ks = "False" 
ks_mod = "Median"  
name_id = "False"  
N = 5  
NSD = 10  
Distance = 250  
block = "True" 
top_block = "True"  
ks_dem_1 = 1.5  
ks_dem_1_1=0.8  
ks_dem_2 = 1.8  
ks_dem_2_2= 1.2 
savefile= 'graph39_11.png'  
```

解释：

```toml
half = "self_col_2"
```

默认参数，不变

```toml
mods = "collinearity"
```

默认参数，不变

```toml
mod = "torder"
```

"torder"或"tpos"        

根据gff文件的 tpos 或 torder 画图 "torder" 根据基因在染色体总的基因上的相对顺序 "tpos" 根据基因在染色体上的真实位置

```toml
Collinearity = 'D:\Python_cod\vivi_coca.collinearity_new.txt'  
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens
Collinearity_2='D:\Python_cod\coca_coca.collinearity.txt' 
blast_gff_1_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1_2 = 'D:\Python_cod\blast文件\coca.lens' 
blast_lens_2_2 = 'D:\Python_cod\blast文件\coca.lens' 
```

写入Collinearity文件的路径

写入Collinearity第1列的gff文件的路径

写入Collinearity第2列的gff文件的路径

写入blast第1列的lens文件的路径

写入Collinearity第2列的lens文件的路径

写入第2个Collinearity文件的路径

写入第2个Collinearity第1列的gff文件的路径

写入第2个Collinearity第2列的gff文件的路径

写入第2个Collinearity第1列的lens文件的路径

写入第2个Collinearity第2列的lens文件的路径

```toml
name_1 = "V.vinifera" 
name_2 = "C.canephora" 
name_1_4 = "B.hispida"  
name_2_4 = "B.hispida" 
```

写入Collinearity第1列物种名称

写入Collinearity第2列物种名称

写入第2个Collinearity第1列物种名称

写入第2个Collinearity第2列物种名称

```toml
posd = "bot_left" 
```

"bot_left"   "bot_right"

图像的格式，第一幅图在左上方还是左下方

```toml
collinearity_ks_mod = "True"  
```

"False" "True"   

是否输入ks文件 

```toml
Collinearity_ks = 'D:\Python_cod\vivi_coca.collinearity.ks.txt'
Collinearity_ks_2= 'D:\Python_cod\coca_coca.collinearity.ks.txt'
```

collinearity_ks_mod="True"情况下 写入第1个ks文件路径

collinearity_ks_mod="True"情况下 写入第2个ks文件路径

```toml
Ks = "True"  
```

"False" "True"

collinearity_ks_mod="True"情况下是否显示ks的值

```toml
ks_mod = "Median"
```

 "Average" 或 "Median"  

collinearity_ks_mod="True"情况下 "Average"显示ks的平均数"Median"显示ks的中位数

```toml
name_id = "False"  
```

"True", "False" 

是否去除物种自身比较时的对角线 (当并非是物种自身比较时 应为 "False" )

```toml
N = 5 
NSD = 10  
Distance = 250 
```

最小对齐数阈值 ≥N 个对齐的区块才合并

显示阈值 ≥NSD 个对齐的区块才显示

区块合并距离 相对顺序差＜250 的区块合并

```toml
block = "True"  
```

"False" "True"  

是否显示 block 的大小

```toml
top_block = "True"  
```

"False" "True"

在collinearity_ks_mod="True" 且Ks="False" 情况下  是否突出显示最近的一次加倍 需要ks文件

```
ks_dem_1 = 1.5  
ks_dem_1_1=0.8  
ks_dem_2 = 1.8 
ks_dem_2_2= 1.2 
```

top_block= "True"情况下 分割第一幅图最近的一次加倍的ks的上限

top_block= "True"情况下 分割第一幅图最近的一次加倍的ks的上限

top_block ="True"情况下 分割第二幅图最近的一次加倍的ks的上限

top_block ="True"情况下 分割第二幅图最近的一次加倍的ks的上限



```toml
savefile= 'graph39_11.png'    
```

写入最后产生的图片路径

###### [blast_Collinearity]

将1个自身与自身blast比对与将1个自身与自身Collinearity比对图以三角形拼接到一起

画blast_Collinearity时的标准格式

```toml
[blast_Collinearity]
half = "self_blco_2"
mods = "blast"
mod = "torder"  
blast_name = 'D:\Python_cod\coca_coca.1e-5.blast'  
blast_gff_1 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\coca.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens'  
Collinearity_2='D:\Python_cod\coca_coca.collinearity.txt' 
blast_gff_1_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1_2 = 'D:\Python_cod\blast文件\coca.lens' 
blast_lens_2_2 = 'D:\Python_cod\blast文件\coca.lens' 
name_1 = "C.canephora"  
name_2 = "C.canephora"  
name_1_4 = "C.canephora"  
name_2_4 = "C.canephora"  
posd = "bot_left"  
collinearity_ks_mod = "False"  
Collinearity_ks_2= 'D:\Python_cod\coca_coca.collinearity.ks.txt'  
Ks = "False"  
ks_mod = "Median"  
name_id = "True"  
N = 5  
NSD = 10  
Distance = 250  
block = "True"  
top_block = "False"  
ks_dem_1 = 1.5  
ks_dem_1_1=0.8  
ks_dem_2 = 1.8  
ks_dem_2_2= 1.2 
savefile= 'graph39_11.png'   
```

解释：

```toml
half = "self_col_2"
```

默认参数，不变

```toml
mods = "blast"
```

默认参数，不变

```toml
mod = "torder"
```

"torder"或"tpos"        

根据gff文件的 tpos 或 torder 画图 "torder" 根据基因在染色体总的基因上的相对顺序 "tpos" 根据基因在染色体上的真实位置

```toml
blast_name = 'D:\Python_cod\vivi_coca.1e-5.blast'  
blast_gff_1 = 'D:\Python_cod\vivi.new.gff'  
blast_gff_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1 = 'D:\Python_cod\blast文件\vivi.lens'  
blast_lens_2 = 'D:\Python_cod\blast文件\coca.lens
Collinearity_2= 'D:\Python_cod\coca_coca.1e-5.blast' 
blast_gff_1_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_gff_2_2 = 'D:\Python_cod\blast例子\Data\coca.new.gff'  
blast_lens_1_2 = 'D:\Python_cod\blast文件\coca.lens' 
blast_lens_2_2 = 'D:\Python_cod\blast文件\coca.lens' 
```

写入blast文件的路径

写入blast第1列的gff文件的路径

写入blast第2列的gff文件的路径

写入blast第1列的lens文件的路径

写入blast第2列的lens文件的路径

写入第2个Collinearity文件的路径

写入第2个Collinearity第1列的gff文件的路径

写入第2个Collinearity第2列的gff文件的路径

写入第2个Collinearity第1列的lens文件的路径

写入第2个Collinearity第2列的lens文件的路径

```toml
name_1 = "V.vinifera" 
name_2 = "C.canephora" 
name_1_4 = "B.hispida"  
name_2_4 = "B.hispida" 
```

写入blast第1列物种名称

写入blast第2列物种名称

写入第2个Collinearity第1列物种名称

写入第2个Collinearity第2列物种名称

```toml
posd = "bot_left" 
```

"bot_left"   "bot_right"

图像的格式，第一幅图在左上方还是左下方

```toml
collinearity_ks_mod = "True"  
```

"False" "True"   

是否输入ks文件 

```toml
Collinearity_ks_2= 'D:\Python_cod\coca_coca.collinearity.ks.txt'
```

collinearity_ks_mod="True"情况下 写入第2个ks文件路径

```toml
Ks = "True"  
```

"False" "True"

collinearity_ks_mod="True"情况下是否显示ks的值

```toml
ks_mod = "Median"
```

 "Average" 或 "Median"  

collinearity_ks_mod="True"情况下 "Average"显示ks的平均数"Median"显示ks的中位数

```toml
name_id = "False"  
```

"True", "False" 

是否去除物种自身比较时的对角线 (当并非是物种自身比较时 应为 "False" )

```toml
N = 5 
NSD = 10  
Distance = 250 
```

最小对齐数阈值 ≥N 个对齐的区块才合并

显示阈值 ≥NSD 个对齐的区块才显示

区块合并距离 相对顺序差＜250 的区块合并

```toml
block = "True"  
```

"False" "True"  

是否显示 block 的大小

```toml
top_block = "True"  
```

"False" "True"

在collinearity_ks_mod="True" 且Ks="False" 情况下  是否突出显示最近的一次加倍 需要ks文件

```
ks_dem_1 = 1.5  
ks_dem_1_1=0.8  
ks_dem_2 = 1.8 
ks_dem_2_2= 1.2 
```

top_block= "True"情况下 分割第一幅图最近的一次加倍的ks的上限

top_block= "True"情况下 分割第一幅图最近的一次加倍的ks的上限

top_block ="True"情况下 分割第二幅图最近的一次加倍的ks的上限

top_block ="True"情况下 分割第二幅图最近的一次加倍的ks的上限

```toml
savefile= 'graph39_11.png'    
```

写入最后产生的图片路径

##### Fragment 配置文件格式

###### config.toml文件

全部参数

```
# config.toml（标准格式）
# 归一化设置
date_norm = "False"     # "False" 不进行归一化，"True" 进行归一化
Max_date = "True"       # "MAX" 理论完美匹配的得分，"True" 实际最大匹配的得分
# DNA 相关参数
A_dna = 25              # 20-30，基准窗口常数
J_dna = 0.22            # 0.15-0.35，谨慎程度
B_dna = -8              # -3 至 -10，补偿强度
C_dna = 200             # 150-300，补偿饱和点
# 蛋白质相关参数
A_protein = 15          # 12-18，基准窗口常数
J_protein = 0.28        # 0.15-0.35，谨慎程度
B_protein = -5          # -3 至 -10，补偿强度
C_protein = 100         # 50-150，补偿饱和点
# 图像显示设置
figs = "False"          # "False" 或 "True"，图像显示时是否使用真实比例
min_prop = "True"       # "False" 或 "True"，图像显示时是否限制下限
min_props = 0           # min_prop="True" 时，下限为多少
# 矩阵选择
Matrix_dna = "BLASTN"           # 选择DNA矩阵："BLASTN"、"TRANS"、"NUC.4.2"
Matrix_protein = "PAM250"       # 选择蛋白质矩阵："BLOSUM62"、"BLOSUM80"、"BLOSUM30"、"PAM250"、"PAM120"、"PAM30"
# 标签分割设置
show_leng_dna = 500             # 设置DNA标签的分割长度
show_leng_protein = 100         # 设置蛋白质标签的分割长度
# 输入输出文件设置
date_1 = "1.pep"
date_2 = "2.pep"
out_part_date_1 = "1.2pep_seq"
out_part_date_2 = "1.2pep_sco"
out_part_date_3 = "1.2pep_pro"
savefile = "output_pep.png"
```

解释

###### **归一化控制**

```
date_norm = "False"  # "False" 或 "True"
```

- **False**：使用原始比对得分
- **True**：将得分归一化到0-1范围
- **作用**：方便不同序列间的比较

###### **最大值设置**

```
Max_date = "True"  # "MAX" 或 "True"
```

- **"MAX"**：使用实际计算得到的最大得分
- **"True"**：使用理论最大可能得分 = 矩阵最大值 × 窗口大小
- **影响**：决定颜色映射的范围上限

###### **图像比例**

```
figs = "False"  # "False" 或 "True"
```

- **False**：使用固定8×8尺寸
- **True**：根据序列长度比例调整图像尺寸
  - 宽度固定为8
  - 高度 = 8 × (序列Y长度 / 序列X长度)

###### **下限控制**

```
min_prop = "True"  # "False" 或 "True"
min_props = 0      # 下限值
```

- **False**：不设下限（默认为0）
- **True**：使用指定的下限值过滤数据点

###### DNA序列参数

###### **窗口计算参数**

```
A_dna = 25      # 20-30  基准窗口常数
J_dna = 0.22    # 0.15-0.35 谨慎程度
B_dna = -8      # -3 至 -10  补偿强度
C_dna = 200     # 150-300    补偿饱和点
```

###### **参数详细解释**

| 参数      | 作用               | 范围      | 设置建议                   |
| :-------- | :----------------- | :-------- | :------------------------- |
| **A_dna** | 基础窗口大小       | 20-30     | 序列复杂时选较大值         |
| **J_dna** | 对长度差异的敏感度 | 0.15-0.35 | 0.22为中等敏感度           |
| **B_dna** | 短序列补偿强度     | -3至-10   | 负值越大，对短序列补偿越强 |
| **C_dna** | 补偿饱和点         | 150-300   | 决定补偿效应的衰减速度     |

###### **DNA比对矩阵选择**

```
Matrix_dna = "BLASTN"  # "BLASTN", "TRANS", "NUC.4.2"
```



- **BLASTN**：标准的DNA比对矩阵
- **TRANS**：转录相关矩阵
- **NUC.4.2**：另一种核苷酸比对矩阵

###### 三、蛋白质序列参数

###### **窗口计算参数**

```
A_protein = 15     # 12-18 基准窗口常数
J_protein = 0.28   # 0.15-0.35 谨慎程度
B_protein = -5     # -3 至 -10  补偿强度
C_protein = 100    # 50-150 补偿饱和点
```

###### **蛋白质比对矩阵选择**

```
Matrix_protein = "PAM250"  # "BLOSUM62", "BLOSUM80", "BLOSUM30", "PAM250", "PAM120", "PAM30"
```

###### **常用矩阵特性**

| 矩阵         | 进化距离 | 适用场景               |
| :----------- | :------- | :--------------------- |
| **BLOSUM62** | 中等     | 通用蛋白质比对（默认） |
| **BLOSUM80** | 较近     | 近缘物种比对           |
| **BLOSUM30** | 较远     | 远缘物种比对           |
| **PAM250**   | 较远     | 25%相同性的序列        |
| **PAM120**   | 中等     | 40%相同性的序列        |
| **PAM30**    | 较近     | 近缘序列比对           |

###### 四、显示参数

###### **标签分割设置**

```
show_leng_dna = 500       # DNA序列刻度标签间隔
show_leng_protein = 100   # 蛋白质序列刻度标签间隔
```

###### **输入文件路径**

```
# 示例路径（需要根据实际情况修改）
date_1 = "1.pep"
date_2 = "2.pep"
```

###### **输出文件路径**

```
out_part_date_1 = "1.2pep_seq"
out_part_date_2 = "1.2pep_sco"
out_part_date_3 = "1.2pep_pro"
savefile = "output_pep.png"
```

## 

##### Intragenic 配置文件格式

config.toml文件

```toml
# config.toml（标准格式）
new_chr = 1
chunk_size = 2000   # 每个片段的长度
start = 20000000   #具体的数字或"start"
end = 30000000      #具体的数字或"end"
date_1 = "chr.fa"
out_part_date_1 = "chr.part.fa"
out_part_date_2 = "alignment_part.sam"
out_part_date_3 = "output_part.csv"
savefile = "output_part.csv.png"
pos = "bot"
```

`new_chr=1`	指定目标染色体的编号,从多染色体的FASTA文件中选择特定染色体进行处理

`chunk_size = 2000 `   每个片段的长度（单位：碱基对 bp,将长染色体序列切割成固定长度的小片段

`start = 20000000`  	提取序列的起始位置（单位：bp）`"start"`：表示从染色体开头（位置0）开始,具体数字：如20000000表示从第2000万个碱基开始

`end = 30000000` 	提取序列的结束位置（单位：bp）`"end"`：表示到染色体末尾结束, 具体数字：如30000000表示到第3000万个碱基结束

out_part_date_1 = "chr.part.fa"           # 步骤1输出：切割后的片段
out_part_date_2 = "alignment_part.sam"    # 步骤2输出：比对结果
out_part_date_3 = "output_part.csv"       # 步骤3输出：统计表格
savefile = "output_part.csv.png"          # 步骤4输出：可视化图像

pos = "bot"	控制**点图可视化时布局的方向**

#### 运行

上面2部完成配置文件后运行 all_run.py 即可



具体配置文件的配置与例图可在对应配置文件目录下的  介绍与使用文档  中查看
