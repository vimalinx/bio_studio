# 酵母菌基因组学习 - 快速开始

## 🚀 3分钟开始

```bash
# 1. 进入项目目录
cd ~/bio_studio/projects/yeast_genome_learning

# 2. 安装数据库（需要5-10分钟）
bash scripts/01_setup_database.sh

# 3. 验证安装
bash scripts/02_verify_install.sh

# 4. 提取你的第一个基因
bash scripts/03_extract_gene.sh ACT1
```

## 📚 你会学到什么

| 技能 | 描述 | 重要性 |
|------|------|--------|
| **基因组数据库** | 理解参考基因组的结构和组织 | ⭐⭐⭐⭐⭐ |
| **序列格式** | FASTA、GFF、GTF 文件格式 | ⭐⭐⭐⭐⭐ |
| **序列工具** | samtools、BLAST 等 | ⭐⭐⭐⭐ |
| **基因注释** | 理解基因功能注释 | ⭐⭐⭐⭐ |
| **数据分析** | 序列统计、比对 | ⭐⭐⭐ |

## 🎯 为什么从酵母开始？

| 特点 | 酵母菌 | 人类 |
|------|--------|------|
| 基因组大小 | 12.1 Mb | 3,200 Mb |
| 下载时间 | ~5分钟 | ~2小时 |
| 分析速度 | 秒级 | 分钟级 |
| 研究深度 | 最深入 | 很深入 |
| 验证难度 | 容易 | 困难 |

**结论**: 用酵母掌握方法，再迁移到人类基因组！

## 📊 验证你的学习成果

### ✅ 第1个验证：ACT1 基因
```bash
bash scripts/03_extract_gene.sh ACT1
```

**预期结果**:
- 系统名称: YFL039C
- 染色体: VI
- 基因长度: 1,128 bp
- 蛋白质长度: 375 aa

### ✅ 第2个验证：基因组统计
```bash
# 在验证脚本中会自动显示：
# - 染色体数量: 17
# - 基因组大小: 12.1 Mb
# - 基因数量: ~6,000
```

### ✅ 第3个验证：数据完整性
```bash
bash scripts/02_verify_install.sh
```

所有检查都应该显示 ✅

## 💡 推荐学习顺序

### Day 1: 安装和验证
- [ ] 运行数据库安装
- [ ] 验证安装成功
- [ ] 提取 ACT1 基因
- [ ] 查看结果文件

### Day 2: 探索数据库
- [ ] 了解目录结构
- [ ] 查看不同文件格式
- [ ] 提取3-5个不同基因
- [ ] 比较它们的大小

### Day 3: 序列分析
- [ ] 计算GC含量
- [ ] 查找开放阅读框（ORF）
- [ ] 序列统计分析
- [ ] 可视化序列分布

### Day 4: 基因注释
- [ ] 理解GFF格式
- [ ] 阅读SGD功能注释
- [ ] 学习基因命名规则
- [ ] GO功能分类

### Day 5: BLAST搜索
- [ ] 使用BLAST搜索
- [ ] 理解E-value和Score
- [ ] 多序列比对
- [ ] 解读比对结果

## 🛠️ 可用工具

### 已安装（Bio Studio环境）
```bash
# 序列操作
samtools     # SAM/BAM操作
seqkit       # FASTA/Q处理

# 序列比对
blastn       # 核苷酸BLAST
blastp       # 蛋白质BLAST

# 其他
makeblastdb  # 创建BLAST数据库
```

### 使用示例
```bash
# 查看序列
samtools faidx data/sequence/genomic.fna NC_001133.9

# 统计序列
samtools faidx data/sequence/genomic.fna | awk '{sum+=$2} END {print sum}'

# 搜索基因
grep 'ACT1' data/annotation/SGD_features.tab
```

## 📖 更多资源

### 完整文档
- 项目README: `~/bio_studio/projects/yeast_genome_learning/README.md`
- 技能文档: `~/bio_studio/.claude/skills/yeast_database/SKILL.md`

### 外部资源
- SGD官网: https://www.yeastgenome.org
- YeastBook: https://www.yeastgenome.org/help/yeast-book
- NCBI教程: https://www.ncbi.nlm.nih.gov/education/

## ❓ 常见问题

**Q: 数据库需要多久下载？**
A: 5-10分钟，取决于网速（~200MB）

**Q: 需要多少空间？**
A: 完整数据库约200MB

**Q: 可以删除数据重新下载吗？**
A: 可以，直接删除 `data/` 目录，重新运行安装脚本

**Q: 酵母菌和人类基因组分析方法一样吗？**
A: 是的！方法完全相同，只是数据规模不同

**Q: 学完酵母后可以做什么？**
A:
- 分析人类基因组
- 比较不同物种
- 做实际研究项目
- 学习宏基因组

## 🎓 完成后的下一步

1. ✅ **巩固基础**: 用酵母练习所有基础操作
2. ✅ **迁移学习**: 下载人类基因组，用相同方法分析
3. ✅ **进阶分析**: 尝试RNA-seq、变异检测等
4. ✅ **实际项目**: 分析你感兴趣的物种或基因

---

**开始你的基因组学之旅吧！** 🍺✨
