import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
new_chr=1
chunk_size = 2000   # 每个片段的长度
start = 20000000       #具体的数字或"start"
end = 30000000   #具体的数字或"end"
date=r"D:\Python_cod\blast例子\Data\pseudomolecules.fa"
out_part_date=fr"D:\Python_cod\blast例子\Data\chr{new_chr}.part.fa"
new_date=[]
for gff3_lod in SeqIO.parse(date, "fasta"):
    if gff3_lod.id ==f"chr{new_chr}" :
        if start=="start":
            start=0
        if end=="end":
            end=len(gff3_lod.seq)
        seq_length = len(gff3_lod.seq)      #得到片段总长
        actual_end = min(end, seq_length)   #得到最终片段位置
        main_seq=gff3_lod.seq[start:actual_end] #第一次切片 得到总的所需片段
        for i, start_part in enumerate(range(0, len(main_seq), chunk_size)):  # 生成起始位置序列 同时获取索引i和起始位置start_part
            end_part = min(start_part + chunk_size, len(main_seq))  # 计算结束位置
            chunk = main_seq[start_part:end_part]  # 直接切片SeqRecord对象
            gff3_new_record = SeqRecord(
                seq=chunk,  # sep=   序列对象
                id=f"{gff3_lod.id}:{start+start_part}:{start+end_part}",  # id=    序列唯一标识符 包含位置信息
                description=''  # 设置description为空字符串
            )  # SeqRecord 建立一个SeqRecord对象
            new_date.append(gff3_new_record)
SeqIO.write(new_date, out_part_date, "fasta")