# Bio Studio MCP 服务器生态系统

为Bio Studio提供的完整MCP（Model Context Protocol）服务器套件。

## 🎯 概述

这套MCP服务器让Claude Code能够直接访问专业的生物信息学工具和数据库，无需手动运行脚本或切换工具。

### 核心优势

✨ **无缝集成** - AI直接调用生物信息学工具
🔬 **专业功能** - 序列分析、结构预测、数据库查询
🚀 **自动执行** - 用自然语言描述，AI自动选择工具
📊 **结果可视化** - 结构化的JSON输出，易于解析

---

## 📦 已实现的MCP服务器

### 1️⃣ bio-sequence-mcp - 序列分析服务器

**功能**:
- DNA/RNA序列分析（GC含量、翻译、互补）
- 蛋白质序列分析（分子量、等电点、疏水性）
- ORF查找
- 序列翻译
- 反向互补序列生成

**使用示例**:
```
用户: 分析这个DNA序列: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

AI: [调用 bio-sequence.analyze_dna]
    结果:
    - 长度: 42 bp
    - GC含量: 52.38%
    - 翻译: MAIVMGR*KGAR*
    - 反向互补: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
```

**位置**: `mcp-servers/bio-sequence-mcp/`

---

### 2️⃣ bio-structure-mcp - 蛋白质结构服务器

**功能**:
- PDB文件解析
- 二级结构分析（α螺旋、β折叠）
- 结构几何指标（质心、尺寸）
- 序列提取
- ESM-Fold结构预测

**使用示例**:
```
用户: 解析这个PDB文件
[上传1AKI.pdb]

AI: [调用 bio-structure.parse_pdb]
    结构信息:
    - 标题: HEN EGG-WHITE LYSOZYME
    - 链数: 1
    - 残基数: 129
    - 原子数: 1001
    - 分辨率: 1.50 Å
```

**位置**: `mcp-servers/bio-structure-mcp/`

---

### 3️⃣ bio-database-mcp - 数据库查询服务器

**功能**:
- PubMed文献搜索
- NCBI核酸/蛋白质数据库
- 基因信息查询
- BLAST序列比对
- UniProt详细信息

**使用示例**:
```
用户: 用BLAST搜索这个序列: MVHLTPEEKSAVTALWGKVN

AI: [调用 bio-database.run_blast]
    BLAST结果:
    - 程序: blastp
    - Top hit: Hemoglobin subunit alpha
    - E-value: 2e-50
    - Identity: 98.6%
```

**位置**: `mcp-servers/bio-database-mcp/`

---

### 4️⃣ bio-design-mcp - 生物设计服务器

**功能**:
- 蛋白质序列反推 DNA 设计
- mRNA 序列设计与表达级别优化
- 靶点可成药性/安全性/竞争格局基础评估
- 设计模块可用性检查

**使用示例**:
```
用户: 给这个蛋白设计一个偏 E.coli 表达的DNA序列: MKT

AI: [调用 bio-design.design_dna]
    结果:
    - DNA序列: ATGAAGACC
    - 长度: 9 bp
    - GC含量: 33.33%
```

**位置**: `mcp-servers/bio-design-mcp/`

---

### 5️⃣ bio-lab-mcp - 工作区实验编排服务器

**功能**:
- 列出工作区项目与统一入口状态
- 读取项目配置、自检摘要和所需工具
- 调用项目级 `validate_project.py`
- 调用工作区级 `workspace-validate` smoke test
- 读取项目 pipeline 步骤说明

**使用示例**:
```
用户: 看一下这个工作区有哪些可运行项目

AI: [调用 bio-lab.list_workspace_projects]
    结果:
    - ai_design_playground
    - test_env_validation
    - yeast_rnaseq_demo
```

**位置**: `mcp-servers/bio-lab-mcp/`

---

## 🧱 预留槽位

以下目录当前是规划位，不是已实现服务：

- 当前没有新的 MCP 预留目录；后续如果接入真实 lab automation 或 protocol 平台，再单独加位。

---

## 🚀 快速开始

### 安装

```bash
# 方法1: 自动安装（推荐）
cd /home/vimalinx/Projects/bio_studio/mcp-servers
bash install-all.sh

# 方法2: 手动安装
pip install -r mcp-requirements.txt
```

### 配置Claude Code

1. 找到Claude Code配置文件:
   - **Linux**: `~/.config/claude-code/config.json`
   - **macOS**: `~/Library/Application Support/Claude Code/config.json`

2. 在当前仓库里生成最新配置:

```bash
python render_claude_config.py --write claude-config.json
cat claude-config.json
```

3. 添加MCP服务器配置（见 `claude-config.json`）:

```json
{
  "mcpServers": {
    "bio-design": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-design-mcp/design_server.py"]
    },
    "bio-lab": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-lab-mcp/lab_server.py"]
    },
    "bio-sequence": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-sequence-mcp/sequence_server.py"]
    },
    "bio-structure": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-structure-mcp/structure_server.py"]
    },
    "bio-database": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-database-mcp/database_server.py"]
    }
  }
}
```

4. 若使用 `bio-database-mcp`，建议先设置:

```bash
export BIO_STUDIO_ENTREZ_EMAIL="you@example.com"
export BIO_STUDIO_NCBI_API_KEY="optional_ncbi_api_key"
```

5. 重启Claude Code

### 验证安装

在Claude Code中测试：

```
用户: 检查MCP服务器是否正常工作

AI: 我可以访问以下生物信息学工具:
    ✅ bio-design-mcp (设计与靶点评估)
    ✅ bio-lab-mcp (项目发现与验证)
    ✅ bio-sequence-mcp (序列分析)
    ✅ bio-structure-mcp (结构分析)
    ✅ bio-database-mcp (数据库查询)
```

---

## 💡 实际应用场景

### 场景1: 快速序列分析

```
用户: 我有一个DNA序列，帮我分析一下
      ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

AI: [自动调用 bio-sequence.analyze_dna]
    分析完成！

    基本信息:
    - 长度: 42 bp
    - GC含量: 52.38%

    翻译结果:
    - 蛋白质序列: MAIVMGR*KGAR*
    - 长度: 14 aa

    其他信息:
    - 互补序列: TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
    - 反向互补: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT

    建议: 这个序列包含终止密码子，可能是基因片段
```

---

### 场景2: 文献调研

```
用户: 帮我找CRISPR-Cas9的最新文献，只看2024年的

AI: [调用 bio-database.search_pubmed]
    找到 1,234 篇2024年发表的CRISPR-Cas9相关文献

    Top 10:
    1. "CRISPR-Cas9 genome editing advances in 2024"
       期刊: Nature Biotechnology
       作者: Smith J, et al.
       年份: 2024

    2. "Clinical applications of CRISPR-Cas9"
       期刊: Science
       ...

    需要我详细分析某篇文献吗？或者下载摘要？
```

---

### 场景3: 蛋白质结构预测

```
用户: 预测这个蛋白质的3D结构:
      MVHLTPEEKSAVTALWGKVN

AI: [调用 bio-structure.predict_structure_esm]
    正在使用ESM-Fold v1预测结构...

    ✅ 预测完成！

    结果摘要:
    - 序列长度: 18 aa
    - 模型: ESM-Fold v1
    - 置信度: High (pLDDT > 90)

    [返回PDB格式坐标]

    后续分析建议:
    1. 查看二级结构
    2. 计算结构指标
    3. 可视化3D结构

    需要我继续分析吗？
```

---

### 场景4: 完整研究流程

```
用户: 帮我研究TP53基因，包括序列、结构、功能

AI: [多步骤工作流]

    步骤1: 查询基因信息
    [调用 bio-database.get_gene_info]
    → TP53: 肿瘤抑制基因，位于17号染色体

    步骤2: 获取蛋白质序列
    [调用 bio-database.search_protein]
    → 找到 UniProt P04637 (393 aa)

    步骤3: 分析序列特征
    [调用 bio-sequence.analyze_protein]
    → 分子量: 43.7 kDa, 等电点: 6.5

    步骤4: 预测3D结构
    [调用 bio-structure.predict_structure_esm]
    → ESM-Fold生成结构模型

    步骤5: 搜索相关文献
    [调用 bio-database.search_pubmed]
    → 找到 15,000+ 相关文献

    ✅ 完整报告已生成！

    需要详细查看哪个部分？
```

---

## 🔧 开发新MCP服务器

### 模板

```python
#!/usr/bin/env python3
"""My Bio MCP Server"""

import asyncio
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

server = Server("my-bio-mcp")

@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="my_tool",
            description="工具描述",
            inputSchema={
                "type": "object",
                "properties": {
                    "input": {"type": "string", "description": "输入参数"}
                },
                "required": ["input"]
            }
        )
    ]

@server.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    # 实现工具逻辑
    result = process(arguments)
    return [TextContent(type="text", result)]

async def main():
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream)

if __name__ == "__main__":
    asyncio.run(main())
```

### 最佳实践

1. **清晰的工具描述** - 帮助AI理解何时使用
2. **结构化输出** - 使用JSON格式
3. **错误处理** - 捕获并返回错误信息
4. **文档完整** - README包含使用示例

---

## 📚 相关文档

- [bio-sequence-mcp README](./bio-sequence-mcp/README.md)
- [bio-structure-mcp README](./bio-structure-mcp/README.md)
- [bio-database-mcp README](./bio-database-mcp/README.md)

---

## 🎓 学习资源

### MCP协议
- [MCP规范](https://modelcontextprotocol.io/)
- [Claude Code文档](https://docs.anthropic.com/claude-code)

### 生物信息学
- [BioPython文档](https://biopython.org/)
- [NCBI资源](https://www.ncbi.nlm.nih.gov/)

---

## 🗺️ 路线图

### v0.2 (计划中)
- [ ] bio-design-mcp - DNA/RNA设计工具
- [ ] bio-lab-mcp - 实验流程自动化
- [ ] bio-variant-mcp - 变异检测

### v0.3 (未来)
- [ ] 本地BLAST支持
- [ ] 结构可视化集成
- [ ] 批量处理模式

---

## ❓ FAQ

**Q: MCP服务器能离线使用吗？**
A: bio-sequence-mcp可以，bio-database-mcp需要网络

**Q: 如何添加新的数据库？**
A: 在bio-database-mcp中添加新工具函数

**Q: 可以同时运行多个MCP服务器吗？**
A: 可以，Claude Code自动管理

**Q: 性能如何？**
A: 本地计算<1秒，数据库查询取决于网络

---

## 🤝 贡献

欢迎贡献新的MCP服务器！

1. 在 `mcp-servers/` 创建新目录
2. 实现MCP服务器
3. 添加文档和示例
4. 更新此README

---

## 📄 许可证

MIT

---

**开始你的AI驱动的生物信息学研究之旅！** 🧬🚀
