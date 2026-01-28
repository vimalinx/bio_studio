# 酵母菌基因组数据库 - 展示网页

## 🌐 在线查看

### 自动打开
网页应该已经在你的浏览器中自动打开！

### 手动打开（如果未自动打开）

```bash
# 方法1: 使用默认浏览器
xdg-open ~/bio_studio/projects/yeast_genome_learning/www/index.html

# 方法2: 在文件管理器中打开
# 导航到: ~/bio_studio/projects/yeast_genome_learning/www/
# 双击: index.html

# 方法3: 在浏览器地址栏输入
file:///home/vimalinx/bio_studio/projects/yeast_genome_learning/www/index.html
```

---

## 📊 网页内容

### 📈 可视化展示
- **基因组概览**: 6个核心统计数据卡片
- **基因功能分类**: 7大类功能的进度条可视化
- **染色体分布**: 16个染色体的3D样式卡片
- **核心生命活动**: 6大系统的功能模块
- **重要基因详解**: 4个经典基因详细介绍
- **酵母 vs 人类**: 对比表格
- **学习建议**: 4周学习路径

### 🎨 设计特点
- ✨ 渐变色彩设计
- 📱 响应式布局（手机/平板/电脑）
- 🎭 悬停动画效果
- 🎯 清晰的数据可视化
- 🌈 渐变进度条

---

## 📂 项目文件结构

```
~/bio_studio/projects/yeast_genome_learning/
├── www/
│   └── index.html          # 🌐 展示网页（本文件）
├── data/                    # 原始数据库
│   ├── sequence/          # 基因组序列
│   ├── annotation/        # 基因注释
│   ├── proteins/          # 蛋白质序列
│   └── blastdb/           # BLAST数据库
├── scripts/                # 分析脚本
│   ├── 01_setup_database.sh
│   ├── 02_verify_install.sh
│   ├── 03_extract_gene.sh
│   └── analyze_genes.sh
├── results/                # 分析结果
│   ├── all_genes.txt      # 6600个基因列表
│   ├── gene_overview.txt  # 可视化概览
│   └── gene_functions_explained.txt
└── logs/                   # 运行日志
```

---

## 🎯 网页功能清单

### 1️⃣ 数据展示
- ✅ 基因组统计（6个数字卡片）
- ✅ 基因功能分类（7个进度条）
- ✅ 染色体分布（17个染色体卡片）
- ✅ 重要基因详解（4个基因）

### 2️⃣ 教育内容
- ✅ 核心生命活动（6大系统）
- ✅ 与人类对比（表格）
- ✅ 基因命名规则（图解）
- ✅ 实用命令（代码块）
- ✅ 学习路径（4周计划）

### 3️⃣ 资源链接
- ✅ SGD官网
- ✅ 基因查询页面
- ✅ 基因本体数据库
- ✅ 本地数据路径

---

## 💡 使用建议

### 作为学习工具
1. **浏览可视化** - 理解基因分布和功能
2. **查询基因** - 使用 all_genes.txt 搜索特定基因
3. **学习命名** - 理解酵母基因命名规则
4. **对比人类** - 了解基因组大小的差异

### 作为分享展示
1. **课堂展示** - 投影讲解酵母基因组
2. **报告生成** - 截图用于文档或PPT
3. **学习记录** - 记录你的学习进度

### 作为开发模板
1. **修改数据** - 更新自己的分析结果
2. **自定义样式** - 修改CSS颜色和布局
3. **添加功能** - 集成更多交互功能

---

## 🛠️ 自定义和扩展

### 修改数据
直接编辑 index.html 中的数字和文本：
```html
<div class="stat-item">
    <div class="number">6,600</div>  <!-- 改成你的数字 -->
    <div class="label">总基因数</div>
</div>
```

### 添加新图表
使用现有的样式类：
```html
<div class="card">
    <h2>📊 你的新图表</h2>
    <!-- 你的内容 -->
</div>
```

### 自定义颜色
修改 CSS 渐变：
```css
.fill-transcription {
    background: linear-gradient(90deg, #你的颜色1, #你的颜色2);
}
```

---

## 📈 数据来源

### 主要数据
- **基因组**: Ensembl Release 112
- **注释**: Ensembl GFF3
- **蛋白质**: Ensembl FASTA
- **SGD**: SGD 功能注释表

### 分析工具
- **samtools**: 序列索引
- **makeblastdb**: BLAST数据库
- **awk/bash**: 数据处理

---

## 🔄 更新数据

### 重新生成网页
如果更新了分析结果，重新运行：

```bash
cd ~/bio_studio/projects/yeast_genome_learning
bash scripts/analyze_genes.sh

# 网页会自动引用最新的结果文件
```

### 手动更新
1. 更新 results/ 中的数据文件
2. 刷新浏览器页面即可看到新数据

---

## 🐛 故障排除

### 网页无法打开
```bash
# 检查文件是否存在
ls -l ~/bio_studio/projects/yeast_genome_learning/www/index.html

# 检查权限
chmod +r ~/bio_studio/projects/yeast_genome_learning/www/index.html
```

### 样式显示异常
- 使用现代浏览器（Chrome、Firefox、Edge）
- 清除浏览器缓存
- 确保JavaScript已启用

---

## 📚 相关文档

- **项目README**: ~/bio_studio/projects/yeast_genome_learning/README.md
- **快速开始**: ~/bio_studio/projects/yeast_genome_learning/QUICKSTART.md
- **技能文档**: ~/bio_studio/.claude/skills/yeast_database/SKILL.md

---

## 🎓 下一步

### 学习方向
1. **深入特定基因** - 研究ACT1、ADH1等
2. **代谢途径** - 分析糖酵解、TCA循环
3. **进化比较** - 与其他酵母菌株对比
4. **迁移学习** - 应用到人类基因组

### 实践项目
1. 提取特定功能的所有基因
2. 分析染色体基因密度差异
3. 构建基因互作网络图
4. 进行序列进化分析

---

**🍺 享受你的基因组学习之旅！**

生成时间: $(date)
版本: v1.0
