import subprocess
from pathlib import Path
import pandas as pd

out_part_date_1=fr"chr.part.fa"
out_part_date_2="alignment_part.sam"

current_vars=['out_part_date_1','out_part_date_2']
try:
    config_file_path = '../../../config.toml'
    config_dir = os.path.dirname(os.path.abspath(config_file_path))

    # 加载配置文件（在上级的上级目录）
    with open(config_file_path, "r", encoding="utf-8") as f:
        config = toml.load(f)

    if "Intragenic" in config:
        intragenic_config = config["Intragenic"]

        for var_name in current_vars:
            if var_name in intragenic_config:
                value = intragenic_config[var_name]  # 直接从配置文件获取
                if isinstance(value, str):
                    if not os.path.isabs(value):  # 相对路径
                        value = os.path.join(config_dir, value)
                globals()[var_name] = value  # 总是赋值

except FileNotFoundError:
    print("⚠ 警告: 未找到 config.toml 文件，使用所有默认配置")
except Exception as e:
    print(f"❌ 错误: 配置文件加载失败 - {str(e)}")
    print("⚠ 将使用默认配置继续运行")



def run_minimap2_simple(minimap2_path, reference, query, output):
    """
    简化的版本，使用pandas写出文件
    """
    # 检查文件是否存在
    for file_path in [minimap2_path, reference, query]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}")
            return None

    cmd = [
        minimap2_path,
        "-ax", "ava-ont",
        "-y", "--cs",
        "--dual=yes", "--eqx",  # 只添加这两个关键参数
        reference,
        query
    ]

    try:
        # 运行minimap2
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # 使用pandas将结果写入文件（最简单的方式）
        # 创建一个单列的DataFrame，每行包含一个SAM记录
        sam_lines = result.stdout.strip().split('\n')
        df = pd.DataFrame({'line': sam_lines})
        df.to_csv(output, index=False, header=False)

        print(f"Success! Output saved to: {output}")
        return output

    except subprocess.CalledProcessError as e:
        print(f"minimap2 failed: {e.stderr}")
        return None


# 使用示例
if __name__ == "__main__":
    result = run_minimap2_simple(
        "./../bin/minimap2-2.30_x64-linux/minimap2",
        out_part_date_1,
        out_part_date_1,
        out_part_date_2
    )
    if result:
        print("Alignment completed successfully!")





