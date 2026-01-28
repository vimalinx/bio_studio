# 共享数据说明

这个目录包含在多个项目中共享的数据资源。

## 目录结构

```
shared_data/
├── references/      # 标准参考序列与索引
│   ├── ecoli_K12.fa
│   ├── ecoli_K12.fna
│   ├── ecoli_K12_index.*.bt2
│   ├── ebov_zaire.fa
│   ├── ebov_zaire.fa.fai
│   └── ebov_zaire_index.*.bt2
│
├── databases/      # 共享数据库（按需添加）
│   ├── nr/
│   ├── nt/
│   └── uniprot/
│
└── tools/          # 外部工具索引（按需添加）
    └── index.json
```

## 添加共享数据

### 添加参考基因组

```bash
# 下载并解压参考基因组
cd shared_data/references
wget http://example.com/hg38.fa.gz
gunzip hg38.fa.gz

# 为BWA创建索引
bwa index hg38.fa
```

### 添加BLAST数据库

```bash
cd shared_data/databases/nr
makeblastdb -in nr.fa -dbtype nucl -parse_seqids -title "NR Database"
```

## 注意事项

- 参考基因组应该定期更新
- BLAST数据库占用空间大，注意存储空间
- 不同的项目可能需要不同的参考序列版本

## 文档

- 每个参考基因组应该有对应的 `.fai` (FAI索引)
- BLAST数据库需要 `.nin`, `.nhr`, `.nsq` 等索引文件
- BWA索引会自动生成 `.bwt`, `.pac`, `.ann`, `.amb` 等文件
