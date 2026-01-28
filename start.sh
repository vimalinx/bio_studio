#!/bin/bash
# Bio Studio 快速启动脚本

echo "🧬 Bio Studio - AI协同生物工作区"
echo "================================"
echo ""

# 激活conda环境（如果存在）
if [ -f "$HOME/miniforge3/bin/activate" ]; then
    source "$HOME/miniforge3/bin/activate" bio 2>/dev/null || \
    source "$HOME/miniforge3/bin/activate" 2>/dev/null
    echo "✅ Conda环境已激活 (Miniforge)"
elif [ -f "$HOME/miniconda3/bin/activate" ]; then
    source "$HOME/miniconda3/bin/activate" bio 2>/dev/null || \
    source "$HOME/miniconda3/bin/activate" 2>/dev/null
    echo "✅ Conda环境已激活 (Miniconda)"
fi

# 显示当前环境
echo ""
echo "📍 当前位置: $(pwd)"
echo "🐍 Python: $(which python)"
echo "📦 环境: ${CONDA_DEFAULT_ENV:-系统Python}"
echo ""

# 检查核心依赖
echo "🔍 检查依赖..."
python -c "import Bio" 2>/dev/null && echo "  ✅ BioPython" || echo "  ❌ BioPython 未安装"
python -c "import pandas" 2>/dev/null && echo "  ✅ Pandas" || echo "  ❌ Pandas 未安装"
python -c "import torch" 2>/dev/null && echo "  ✅ PyTorch" || echo "  ⚠️  PyTorch 未安装 (AI功能需要)"
python -c "import esm" 2>/dev/null && echo "  ✅ ESM (蛋白质模型)" || echo "  ⚠️  ESM 未安装 (蛋白质预测需要)"

echo ""
echo "💡 使用提示:"
echo "  1. 将数据放入 data/raw/"
echo "  2. 在 projects/ 创建项目"
echo "  3. 告诉AI你要做什么"
echo ""
echo "🚀 准备就绪！开始你的研究吧！"
echo ""
