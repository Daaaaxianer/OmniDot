import glob
import os
import sys
import re
import pysam
import pandas as pd
out_part_date_2="alignment_part.sam"
out_part_date_3='output_part.csv'
current_vars=['out_part_date_2','out_part_date_3']
try:
    # 加载TOML配置文件
    with open('config.toml', "r", encoding="utf-8") as f:
        config = toml.load(f)
    if "Intragenic" in config:
        intragenic_config = config["Intragenic"]
        for var_name in current_vars:# 遍历所有需要检查的变量
            if var_name in intragenic_config:
                globals()[var_name] = intragenic_config[var_name]
except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")

samfile = pysam.AlignmentFile(out_part_date_2, threads=4)
data_list = []
for read in samfile:
    if read.flag == 4:
        continue
    # 初始化变量
    matches = mismatch = ins = dele = insEvent = delEvent = 0
    query_name = read.query_name
    query_alignment_start = read.query_alignment_start
    query_alignment_end = read.query_alignment_end
    query_length = read.infer_query_length()
    reference_name = read.reference_name
    reference_start = read.reference_start
    reference_end = read.reference_end

    cs = read.get_tag("cs")
    pattern = ":[0-9]+"
    match = re.findall(pattern, cs)
    if match is not None:
        match = [int(x[1:]) for x in match]
        matches = sum(match)

    pattern = "-[a-z]+"
    match = re.findall(pattern, cs)
    if match is not None:
        match = [len((x[1:])) for x in match]
        delEvent = len(match)
        dele = sum(match)

    pattern = "\+[a-z]+"
    match = re.findall(pattern, cs)
    if match is not None:
        match = [len((x[1:])) for x in match]
        insEvent = len(match)
        ins = sum(match)

    mismatch = cs.count("*")
    strand = "+"
    if read.is_reverse:
        strand = "-"
    if (matches + mismatch) == 0:
        continue
    bymatches = (100.0 * matches) / (matches + mismatch)
    byevents = (100.0 * matches) / (matches + mismatch + insEvent + delEvent)
    byall = (100.0 * matches) / (matches + mismatch + ins + dele)
    data_list.append({
        'query_name': query_name,
        'query_alignment_start': query_alignment_start,
        'query_alignment_end': query_alignment_end,
        'query_length': query_length,
        'reference_name': reference_name,
        'reference_alignment_start': reference_start,
        'reference_alignment_end': reference_end,
        'bymatches': bymatches,     #基于匹配的同一性%
        'byevents': byevents,       #基于事件的同一性%
        'byall': byall,             #总体同一性%
        'matches': matches,         #匹配碱基数
        'mismatch': mismatch,       #错配碱基数
        'ins': ins,                 #插入碱基数
        'dele': dele,               #缺失碱基数
        'insEvent': insEvent,       #插入事件数
        'delEvent': delEvent,       #缺失事件数
        'strand': strand            #链方向
    })
    # 创建dataframe
df = pd.DataFrame(data_list)
df.sort_values(by=["query_name", "reference_name", "matches"], ascending=False)
df.drop_duplicates(subset=["query_name", "reference_name"], inplace=True) #去重 只考虑 query_name 和 reference_name 组合的唯一性
df['query_start'] = df['query_name'].str.split(':').str[1].str.replace('seq', '').astype(int)
df['query_end'] = df['query_name'].str.split(':').str[2].astype(int)
df['reference_start']=df['reference_name'].str.split(':').str[1].str.replace('seq', '').astype(int)
df['reference_end']=df['reference_name'].str.split(':').str[2].astype(int)
out = df
#对称化处理
out2 = out.copy()
out2.reference_name = out.query_name
out2.reference_start = out.query_start
out2.reference_end = out.query_end
out2.query_name = out.reference_name
out2.query_start = out.reference_start
out2.query_end = out.reference_end
out = (
            pd.concat([out, out2])
            .sort_values(
                by=[
                    "query_name",
                    "query_start",
                    "reference_name",
                    "reference_start",
                    "byevents",
                ],
                ascending=False,
            )
        )
columns_to_keep = [
    'query_name', 'query_start', 'query_end',
    'reference_name', 'reference_start', 'reference_end',
    'byevents', 'strand'
]
out = out.loc[:, columns_to_keep]
out.to_csv('output_part.csv', sep='\t', index=False,header=True)