# Bio Sequence MCP Server

DNA/RNA/蛋白质序列分析MCP服务器，为Bio Studio提供核心序列操作功能。

## 功能特性

### 🔬 序列分析工具

1. **analyze_dna** - DNA序列分析
   - GC含量计算
   - 翻译为蛋白质
   - 互补序列生成
   - 反向互补序列

2. **analyze_rna** - RNA序列分析
   - GC含量计算
   - 翻译为蛋白质
   - 反向转录为DNA

3. **analyze_protein** - 蛋白质序列分析
   - 分子量计算
   - 等电点预测
   - 不稳定性指数
   - 氨基酸组成
   - 疏水性(GRAVY)

4. **find_orfs** - 查找开放阅读框
   - 三框翻译
   - 起始/终止密码子识别
   - 可配置最小长度

5. **translate_sequence** - 序列翻译
   - 支持多种遗传密码表
   - 在终止密码子处停止

6. **reverse_complement** - 反向互补序列生成

## 安装

```bash
# 安装依赖
pip install -e .

# 或使用
pip install mcp biopython
```

## 使用方法

### 配置Claude Code

在Claude Code配置文件中添加：

```json
{
  "mcpServers": {
    "bio-sequence": {
      "command": "python",
      "args": ["/path/to/bio_studio/mcp-servers/bio-sequence-mcp/sequence_server.py"]
    }
  }
}
```

### 示例对话

```
用户: 分析这个DNA序列: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

AI: [调用 analyze_dna 工具]
    分析结果:
    - 长度: 42 bp
    - GC含量: 52.38%
    - 翻译: MAIVMGR*KGAR*
    - 反向互补: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
```

```
用户: 查找这个序列中的ORF: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

AI: [调用 find_orfs 工具]
    发现 2 个ORF:
    - 框架+1: 起始0, 终止39, 长度39
      蛋白质: MAIVMGR*KGAR*
```

```
用户: 这个蛋白质的分子量和等电点? MVHLTPEEKSAVTALWGKVN

AI: [调用 analyze_protein 工具]
    - 分子量: 1785.03 Da
    - 等电点: 8.52
    - 疏水性: -0.65 (亲水)
```

## 技术细节

### 依赖
- **mcp**: Model Context Protocol SDK
- **biopython**: 生物信息学Python库

### 实现要点
- 使用asyncio异步处理
- 标准MCP协议实现
- 完整的错误处理
- JSON格式输出

## 扩展开发

可以添加的功能：
- FASTA文件解析
- 序列比对 (BLAST, pairwise)
- 序列motif搜索
- 限制性酶切位点分析
- PCR引物设计

## 许可证

MIT
