#!/bin/bash
# Bio Studio - 快速添加工具
# 用法: ./add-tool.sh <GitHub仓库URL> [工具名称]

set -e

BIO_STUDIO="$HOME/bio_studio"
REPOS_DIR="$BIO_STUDIO/repositories/active"

# 颜色输出
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 Bio Studio - 添加工具${NC}"
echo "================================"

# 检查参数
if [ -z "$1" ]; then
    echo "用法: $0 <GitHub仓库URL> [工具名称]"
    echo ""
    echo "示例:"
    echo "  $0 https://github.com/facebookresearch/esm"
    echo "  $0 https://github.com/facebookresearch/esm esm-model"
    echo ""
    exit 1
fi

REPO_URL="$1"
TOOL_NAME="${2:-}"

# 提取仓库名称
if [ -z "$TOOL_NAME" ]; then
    TOOL_NAME=$(basename "$REPO_URL" .git)
fi

echo -e "${YELLOW}📥 克隆仓库: $TOOL_NAME${NC}"
echo "URL: $REPO_URL"
echo ""

# 克隆仓库
git clone "$REPO_URL" "$REPOS_DIR/$TOOL_NAME"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ 仓库已克隆到: $REPOS_DIR/$TOOL_NAME${NC}"
    echo ""
    echo "🔍 分析工具并生成Skill..."

    # 调用Python脚本生成Skill
    cd "$BIO_STUDIO"
    python3 tools/scripts/deploy_tool.py --repo "$REPOS_DIR/$TOOL_NAME"

    echo ""
    echo -e "${GREEN}🎉 完成！${NC}"
    echo ""
    echo "📝 下一步:"
    echo "  1. 编辑生成的Skill文档: $BIO_STUDIO/.claude/skills/$TOOL_NAME/SKILL.md"
    echo "  2. 添加详细的使用说明和示例"
    echo "  3. 测试工具是否正常工作"
    echo ""
else
    echo -e "${YELLOW}❌ 克隆失败${NC}"
    exit 1
fi
