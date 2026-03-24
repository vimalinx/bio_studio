# Bio Design MCP Server

生物设计类 MCP server，封装当前工作区里的 DNA 设计、mRNA 优化和靶点评估脚本。

## 当前能力

1. `design_dna`
   - 输入蛋白质序列
   - 生成面向目标物种的 DNA 序列
   - 返回 GC、重复序列、限制性位点等基础分析
   - 可选附带 PCR 引物建议

2. `design_mrna`
   - 输入蛋白质序列
   - 按 `low` / `medium` / `high` 表达目标生成 mRNA
   - 返回稳定性、GC 和折叠能等分析结果

3. `analyze_target`
   - 对候选靶点做基础可成药性分析
   - 返回安全性、疾病相关性和竞争格局摘要

4. `get_server_capabilities`
   - 查看当前绑定到哪些本地设计脚本
   - 检查模块是否可导入

## 脚本来源

当前 server 直接复用这些本地脚本：

- `tools/scripts/dna_design.py`
- `tools/scripts/mrna_optimize.py`
- `tools/scripts/target_analysis.py`

## 配置

建议在仓库根目录重新渲染 Claude 配置：

```bash
python mcp-servers/render_claude_config.py --write mcp-servers/claude-config.json
```

最小配置项如下：

```json
{
  "mcpServers": {
    "bio-design": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-design-mcp/design_server.py"]
    }
  }
}
```

## 示例

### DNA 设计

```
用户: 给这个蛋白设计一个 E.coli 表达用的 DNA 序列: MKT

AI: [调用 design_dna]
    - DNA序列: ATGAAGACC
    - 长度: 9 bp
    - GC含量: 33.33%
```

### mRNA 设计

```
用户: 按中等表达强度给这个蛋白生成 mRNA: MKT

AI: [调用 design_mrna]
    - mRNA序列: GCCACCAUGAAGACCAATAAA
    - 稳定性评分: 0.7
```

### 靶点评估

```
用户: 帮我看一下 TP53 作为肿瘤靶点值不值得继续分析

AI: [调用 analyze_target]
    - 输出综合评分、风险等级、疾病相关性和竞争格局摘要
```

## 注意事项

- `analyze_target` 会访问外部 API，结果依赖网络和第三方服务状态。
- 当前实现偏向工作区内脚本封装，不是严格标准化的设计平台 API。
