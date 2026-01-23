import subprocess
from pathlib import Path
import pandas as pd


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
        "chr1.part1.fa",
        "chr1.part1.fa",
        "alignment_part.sam"
    )
    if result:
        print("Alignment completed successfully!")





