#!/bin/bash
# Biomni 快速启动脚本

set -e

echo "========================================"
echo "Biomni API 配置管理"
echo "========================================"
echo

# 检查conda环境
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "⚠️  警告: Conda环境未激活"
    echo "请先运行: conda activate biomni_e1"
    echo
    read -p "是否继续? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 运行配置切换工具
python switch_profile.py

echo
echo "========================================"
echo "配置完成！"
echo "========================================"
echo
echo "下一步："
echo "1. 如未配置API密钥，请编辑 .env 文件"
echo "2. 运行 Biomni:"
echo
echo "   Python代码:"
echo "   from biomni.agent import A1"
echo "   agent = A1(path='./data', llm='claude-sonnet-4-5')"
echo "   agent.go('你的生物医学任务')"
echo
echo "   或启动Gradio界面:"
echo "   from biomni.agent import A1"
echo "   agent = A1(path='./data', llm='claude-sonnet-4-5')"
echo "   agent.launch_gradio_demo()"
echo
