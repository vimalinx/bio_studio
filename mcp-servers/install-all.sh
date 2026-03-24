#!/bin/bash
# 安装所有Bio Studio MCP服务器

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIO_STUDIO="$(cd "$SCRIPT_DIR/.." && pwd)"
MCP_DIR="$SCRIPT_DIR"

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

# 安装requests
echo "📦 安装requests..."
pip install requests>=2.31.0

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
python3 "$MCP_DIR/render_claude_config.py" --write "$MCP_DIR/claude-config.json"
echo "请将以下配置添加到Claude Code配置文件中:"
echo ""
echo "Linux: ~/.config/claude-code/config.json"
echo "macOS: ~/Library/Application Support/Claude Code/config.json"
echo ""
cat "$MCP_DIR/claude-config.json"
echo ""
echo "数据库查询前请设置 NCBI 身份："
echo "  export BIO_STUDIO_ENTREZ_EMAIL='you@example.com'"
echo "可选 API key："
echo "  export BIO_STUDIO_NCBI_API_KEY='your_ncbi_api_key'"
echo ""

# 创建测试脚本
cat > "$MCP_DIR/test-mcp.sh" << 'EOF'
#!/bin/bash
# 测试MCP服务器

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="$(python3 "$SCRIPT_DIR/render_claude_config.py" --python-path)"

echo "测试 bio-sequence-mcp..."
"$PYTHON_BIN" "$SCRIPT_DIR/bio-sequence-mcp/sequence_server.py" &
PID1=$!

echo "测试 bio-design-mcp..."
"$PYTHON_BIN" "$SCRIPT_DIR/bio-design-mcp/design_server.py" &
PID4=$!

echo "测试 bio-lab-mcp..."
"$PYTHON_BIN" "$SCRIPT_DIR/bio-lab-mcp/lab_server.py" &
PID5=$!

echo "测试 bio-structure-mcp..."
"$PYTHON_BIN" "$SCRIPT_DIR/bio-structure-mcp/structure_server.py" &
PID2=$!

echo "测试 bio-database-mcp..."
"$PYTHON_BIN" "$SCRIPT_DIR/bio-database-mcp/database_server.py" &
PID3=$!

echo ""
echo "所有MCP服务器已启动（进程ID: $PID1, $PID4, $PID5, $PID2, $PID3）"
echo "按Ctrl+C停止所有服务器"

trap "kill $PID1 $PID4 $PID5 $PID2 $PID3 2>/dev/null" EXIT

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
