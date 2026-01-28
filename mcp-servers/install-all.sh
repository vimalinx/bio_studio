#!/bin/bash
# 安装所有Bio Studio MCP服务器

set -e

BIO_STUDIO="/media/vimalinx/Data/bio_studio"
MCP_DIR="$BIO_STUDIO/mcp-servers"

echo "🧬 Bio Studio MCP 服务器安装"
echo "================================"
echo ""

# 检查Python版本
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "📍 Python版本: $PYTHON_VERSION"
echo ""

# 安装MCP核心库
echo "📦 安装MCP核心库..."
pip install mcp>=0.9.0

# 安装BioPython
echo "📦 安装BioPython..."
pip install biopython>=1.83

# 安装numpy
echo "📦 安装numpy..."
pip install numpy>=1.24.0

# 询问是否安装ESM（结构预测）
echo ""
read -p "是否安装ESM-Fold（蛋白质结构预测）? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "📦 安装ESM和PyTorch（可能需要几分钟）..."
    pip install fair-esm>=2.0.0 torch>=2.0.0
    echo "✅ ESM-Fold已安装"
fi

echo ""
echo "🔧 配置MCP服务器..."
echo ""
echo "请将以下配置添加到Claude Code配置文件中:"
echo ""
echo "Linux: ~/.config/claude-code/config.json"
echo "macOS: ~/Library/Application Support/Claude Code/config.json"
echo ""
cat "$MCP_DIR/claude-config.json"
echo ""

# 创建测试脚本
cat > "$MCP_DIR/test-mcp.sh" << 'EOF'
#!/bin/bash
# 测试MCP服务器

echo "测试 bio-sequence-mcp..."
python3 /media/vimalinx/Data/bio_studio/mcp-servers/bio-sequence-mcp/sequence_server.py &
PID1=$!

echo "测试 bio-structure-mcp..."
python3 /media/vimalinx/Data/bio_studio/mcp-servers/bio-structure-mcp/structure_server.py &
PID2=$!

echo "测试 bio-database-mcp..."
python3 /media/vimalinx/Data/bio_studio/mcp-servers/bio-database-mcp/database_server.py &
PID3=$!

echo ""
echo "所有MCP服务器已启动（进程ID: $PID1, $PID2, $PID3）"
echo "按Ctrl+C停止所有服务器"

trap "kill $PID1 $PID2 $PID3 2>/dev/null" EXIT

wait
EOF

chmod +x "$MCP_DIR/test-mcp.sh"

echo ""
echo "✅ 安装完成！"
echo ""
echo "📝 下一步:"
echo "  1. 配置Claude Code（见上方配置）"
echo "  2. 重启Claude Code"
echo "  3. 开始使用MCP工具"
echo ""
echo "🧪 测试服务器:"
echo "  bash $MCP_DIR/test-mcp.sh"
echo ""
