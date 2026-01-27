# run_all_4.py
import os

print("开始运行基因组分析全流程...")
print("=" * 60)

# 1. 切割基因组
print("\n[1/4] 运行 make_date_1.py...")
os.system("python make_date_1.py")
print("✓ 步骤1完成")

# 2. minimap2比对
print("\n[2/4] 运行 make_minimap_2.py...")
os.system("python make_minimap_2.py")
print("✓ 步骤2完成")

# 3. 解析SAM文件
print("\n[3/4] 运行 make_date_bed_3.py...")
os.system("python make_date_bed_3.py")
print("✓ 步骤3完成")

# 4. 生成可视化
print("\n[4/4] 运行 make_fignew_4.py...")
os.system("python make_fignew_4.py")
print("✓ 步骤4完成")

print("\n" + "=" * 60)
print("所有4个步骤都运行完成！")