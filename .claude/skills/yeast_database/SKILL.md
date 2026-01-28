# 酵母菌基因组数据库管理技能

## 技能描述
管理酵母菌（*Saccharomyces cerevisiae*）基因组数据库的下载、验证、查询和分析。

## 文件位置
`~/bio_studio/projects/yeast_genome_learning/`

## 可用脚本

### 1. 数据库安装
**脚本**: `scripts/01_setup_database.sh`

**功能**:
- 下载酵母菌参考基因组（R64-1-1）
- 下载基因注释（GFF格式）
- 下载蛋白质序列
- 下载SGD功能注释
- 创建BLAST数据库

**使用**:
```bash
cd ~/bio_studio/projects/yeast_genome_learning
bash scripts/01_setup_database.sh
```

**预计时间**: 5-10分钟
**占用空间**: ~200MB

**输出**:
```
data/
├── sequence/genomic.fna      # 基因组序列
├── annotation/genomic.gff    # 基因注释
├── proteins/protein.faa      # 蛋白质序列
├── annotation/SGD_features.tab  # 功能注释
└── blastdb/                  # BLAST数据库
```

---

### 2. 验证安装
**脚本**: `scripts/02_verify_install.sh`

**功能**:
- 检查所有必需文件
- 统计基因组信息
- 验证ACT1基因（标准测试）
- 检查BLAST数据库

**使用**:
```bash
bash scripts/02_verify_install.sh
```

**验证点**:
- ✅ 基因组大小: ~12.1 Mb
- ✅ 染色体数量: 17 (16 + mt)
- ✅ 基因数量: ~6,000
- ✅ 蛋白数量: ~6,000
- ✅ ACT1基因存在

---

### 3. 提取基因序列
**脚本**: `scripts/03_extract_gene.sh <GENE_NAME>`

**功能**:
- 查找基因信息
- 提取DNA序列
- 提取蛋白质序列
- 显示序列统计

**使用**:
```bash
# 提取ACT1基因
bash scripts/03_extract_gene.sh ACT1

# 列出可用基因
bash scripts/03_extract_gene.sh
```

**输出**:
```
results/<GENE_NAME>.fa           # DNA序列
results/<GENE_NAME>_protein.faa  # 蛋白质序列
```

**常用基因**:
- `ACT1` - 肌动蛋白（经典测试基因）
- `ADH1` - 酒精脱氢酶
- `HIS3` - 组氨酸合成
- `URA3` - 尿嘧啶合成

---

## 数据库内容

### 基因组信息
- **物种**: *Saccharomyces cerevisiae* (酿酒酵母)
- **菌株**: S288C
- **版本**: R64-1-1 (Release 64)
- **基因组大小**: 12,147,813 bp
- **染色体**: 16条核染色体 + 1个线粒体
- **基因**: ~6,607个蛋白编码基因
- **蛋白质**: ~6,031个蛋白质

### 文件格式说明

#### FASTA (`.fna`, `.faa`)
```
>序列标识符
ATCGATCGATCG...
```

#### GFF/GTF (`.gff`, `.gtf`)
```
染色体  来源  类型  起始  结束  分数  链  相位  属性
```

#### SGD Features Tab
```
染色体  起始  结束  基因类型  系统名  标准名  ...
```

---

## 常用分析命令

### 使用 samtools
```bash
# 查看所有染色体
samtools faidx data/sequence/genomic.fna

# 提取特定区域
samtools faidx data/sequence/genomic.fna NC_001133.9:1000-2000

# 统计基因组大小
samtools faidx data/sequence/genomic.fna | awk '{sum+=$2} END {print sum}'
```

### 查询基因信息
```bash
# 从GFF查找
grep 'ACT1' data/annotation/genomic.gff

# 从SGD表查找
grep 'ACT1' data/annotation/SGD_features.tab

# 查找特定染色体基因
awk '$1=="NC_001133.9"' data/annotation/SGD_features.tab
```

### 统计分析
```bash
# 统计基因数量
grep -w 'gene' data/annotation/genomic.gff | wc -l

# 按染色体统计
awk '$3=="gene" {print $1}' data/annotation/genomic.gff | sort | uniq -c

# 统计不同类型
awk '{print $3}' data/annotation/genomic.gff | sort | uniq -c | sort -rn
```

### BLAST搜索
```bash
# 核苷酸搜索
blastn -db data/blastdb/yeast_genome -query query.fa -out results.txt

# 蛋白质搜索
blastp -db data/blastdb/yeast_protein -query query.faa -out results.txt

# 详细输出格式
blastn -db data/blastdb/yeast_genome -query query.fa -outfmt 7
```

---

## 学习路径

### 第1周: 数据库基础
- [x] 安装数据库
- [x] 验证安装
- [x] 理解文件格式
- [ ] 提取多个基因序列

### 第2周: 序列分析
- [ ] 计算GC含量
- [ ] 查找ORF
- [ ] 序列统计
- [ ] 绘制序列分布图

### 第3周: 基因注释
- [ ] 理解GFF格式
- [ ] 解读SGD注释
- [ ] GO功能分类
- [ ] 基因命名规则

### 第4周: BLAST分析
- [ ] 序列比对
- [ ] 理解E-value
- [ ] 多序列比对
- [ ] 进化分析

---

## 参考资源

### 官方数据库
- **SGD**: https://www.yeastgenome.org
- **NCBI**: https://www.ncbi.nlm.nih.gov/genome/?term=Saccharomyces+cerevisiae
- **Ensembl Fungi**: https://fungi.ensembl.org/Saccharomyces_cerevisiae/

### 学习文档
- **项目README**: `~/bio_studio/projects/yeast_genome_learning/README.md`
- **Bio Studio规范**: `~/bio_studio/README.md`
- **避坑指南**: `~/bio_studio/docs/BEST_PRACTICES.md`

---

## 故障排除

### 下载失败
**问题**: wget 下载中断
**解决**:
```bash
# 脚本已使用 -c 参数，支持断点续传
# 直接重新运行即可
bash scripts/01_setup_database.sh
```

### 索引创建失败
**问题**: samtools 未安装
**解决**:
```bash
conda activate bio
conda install -c bioconda samtools
```

### BLAST未安装
**问题**: makeblastdb 命令不存在
**解决**:
```bash
conda activate bio
conda install -c bioconda blast
```

---

## 更新日志

**v1.0 (2026-01-25)**:
- ✅ 初始版本
- ✅ 数据库下载脚本
- ✅ 验证脚本
- ✅ 基因提取脚本
