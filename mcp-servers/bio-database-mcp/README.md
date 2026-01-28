# Bio Database MCP Server

生物数据库查询MCP服务器，连接NCBI、UniProt等生物数据库。

## 功能特性

### 🔗 数据库查询工具

1. **search_pubmed** - PubMed文献搜索
   - 科学论文检索
   - 支持布尔查询
   - 返回标题、作者、期刊信息

2. **search_nucleotide** - 核酸序列数据库
   - GenBank查询
   - RefSeq搜索
   - 序列元数据提取

3. **search_protein** - 蛋白质数据库
   - NCBI蛋白质
   - UniProt集成
   - 功能注释

4. **get_gene_info** - 基因信息查询
   - 基因符号查询
   - 染色体定位
   - 功能描述

5. **run_blast** - BLAST序列比对
   - NCBI在线BLAST
   - 支持blastp/blastn
   - 同源性搜索

6. **get_uniprot_info** - UniProt详细信息
   - 蛋白质功能
   - 亚细胞定位
   - 分子量计算

## 安装

```bash
pip install -e .
```

**重要**: 首次使用前需要配置NCBI邮箱（在代码中设置）：
```python
from Bio import Entrez
Entrez.email = "your-email@example.com"
```

## 使用方法

### 配置Claude Code

```json
{
  "mcpServers": {
    "bio-database": {
      "command": "python",
      "args": ["/path/to/bio_studio/mcp-servers/bio-database-mcp/database_server.py"]
    }
  }
}
```

### 示例对话

```
用户: 搜索CRISPR相关的最新文献

AI: [调用 search_pubmed]
    找到 15,234 篇相关文献:
    1. "CRISPR-Cas9 genome editing" - Nature, 2024
    2. "Base editing with CRISPR" - Science, 2024
    ...
```

```
用户: 查询TP53基因的信息

AI: [调用 get_gene_info]
    基因信息:
    - 名称: TP53
    - 描述: Tumor protein p53
    - 染色体: 17
    - 位置: 17p13.1
    - 生物: 人类
    - 别名: P53
```

```
用户: 搜索血红蛋白的蛋白质信息

AI: [调用 search_protein]
    找到 5 个匹配:
    1. Hemoglobin subunit alpha (P69905)
       - 生物: Homo sapiens
       - 长度: 142 aa
    ...
```

```
用户: 用BLAST搜索这个序列: MVHLTPEEKSAVTALWGKVN

AI: [调用 run_blast]
    BLAST结果:
    - 程序: blastp
    - 数据库: nr
    - 发现 23 个匹配
    Top hit: Hemoglobin subunit alpha
      E-value: 2e-50
      Identity: 98.6%
```

```
用户: 查看UniProt P04637的详细信息

AI: [调用 get_uniprot_info]
    蛋白质信息:
    - 登录号: P04637
    - 名称: Tumor suppressor p53
    - 基因: TP53
    - 长度: 393 aa
    - 分子量: 43.7 kDa
    - 定位: Nucleus, Cytoplasm
```

## 数据库API限制

### NCBI Entrez
- 无API密钥: 3请求/秒
- 有API密钥: 10请求/秒
- 需要提供邮箱地址

### UniProt REST
- 无严格限制
- 建议合理使用

### BLAST
- 单次查询最长1000aa
- 可能需要等待30秒-2分钟

## 技术细节

### 依赖
- **biopython**: Bio.Entrez, Bio.Blast
- **requests**: UniProt API调用
- **mcp**: MCP协议

### 网络要求
需要稳定的互联网连接访问：
- https://eutils.ncbi.nlm.nih.gov
- https://blast.ncbi.nlm.nih.gov
- https://rest.uniprot.org

## 最佳实践

1. **使用具体的关键词**
   ```
   好: "CRISPR AND 2024[PDAT]"
   差: "CRISPR"
   ```

2. **设置合理的max_results**
   - 文献搜索: 10-20
   - 序列搜索: 5-10

3. **利用BLAST限制**
   - 使用序列片段而非全长
   - 选择合适的数据库

## 常见问题

**Q: BLAST超时怎么办？**
A: 序列太长，使用片段或本地BLAST

**Q: 如何获取NCBI API密钥？**
A: 在 https://www.ncbi.nlm.nih.gov/account/ 注册

**Q: 搜索不到结果？**
A: 检查拼写，尝试同义词或更宽泛的关键词

## 未来扩展

- 本地BLAST数据库支持
- Ensembl数据库集成
- PDB结构数据库
- 药物数据库查询

## 许可证

MIT
