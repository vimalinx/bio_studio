#!/bin/bash
set -e  # 遇到错误立即停止

# 定义环境路径
PYTHON_EXEC="/home/vimalinx/miniforge3/envs/bio/bin/python"
PIP_EXEC="/home/vimalinx/miniforge3/envs/bio/bin/pip"

echo "🧬 Bio Studio 环境补全脚本"
echo "=============================="
echo "目标环境: $PYTHON_EXEC"
echo ""

# 检查pip是否存在
if [ ! -f "$PIP_EXEC" ]; then
    echo "❌ 错误: 找不到 pip ($PIP_EXEC)"
    echo "请检查 conda 环境 'bio' 是否正确创建。"
    exit 1
fi

# 定义要安装的包列表
PACKAGES=(
    "pandas"
    "scipy"
    "matplotlib"
    "seaborn"
    "scikit-learn"
    "jupyter"
    "notebook"
    "pysam"
    "scanpy"  # 推荐：单细胞分析标准
    "openpyxl" # 推荐：读写Excel
)

echo "📦 准备安装以下软件包:"
for pkg in "${PACKAGES[@]}"; do
    echo "  - $pkg"
done
echo ""

# 执行安装
echo "🚀 开始安装..."
"$PIP_EXEC" install "${PACKAGES[@]}"

echo ""
echo "✨ 所有依赖安装完成！"
echo "📊 验证安装版本:"
"$PIP_EXEC" list | grep -E "pandas|scipy|matplotlib|seaborn|scikit-learn|jupyter|pysam|scanpy"

echo ""
echo "✅ 环境补全成功。你可以使用 './start.sh' 启动工作区了。"
