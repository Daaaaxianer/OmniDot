import os
import toml
config_1 = "Chromosomal"    #3种模式 "Chromosomal" "Fragment" "Intragenic"
config_2 = "blast"          # 在"Chromosomal" 下的哪一个板块 "blast","Collinearity","blast_blast","Collinearity_Collinearity","blast_Collinearity"

current_vars=['config_1','config_2']
try:
    # 加载TOML配置文件
    with open('config.toml', "r", encoding="utf-8") as f:
        config = toml.load(f)
    if "Config" in config:
        intragenic_config = config["Config"]
        for var_name in current_vars:# 遍历所有需要检查的变量
            if var_name in intragenic_config:
                globals()[var_name] = intragenic_config[var_name]
except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")

if config_1== "Chromosomal":
    os.system('python "Chromosomal Level/graph39.py"')
if config_1== "Fragment":
    os.system('python "Fragment Level/dotter5.py"')
if config_1== "Intragenic":
    os.system('python "Intragenic Level/project/src/run.py"')