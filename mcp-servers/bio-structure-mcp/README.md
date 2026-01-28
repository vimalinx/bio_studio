# Bio Structure MCP Server

蛋白质结构分析MCP服务器，为Bio Studio提供结构生物学分析功能。

## 功能特性

### 🧬 结构分析工具

1. **parse_pdb** - 解析PDB文件
   - 提取链信息
   - 统计残基和原子数量
   - 显示分辨率和实验方法

2. **analyze_secondary_structure** - 二级结构分析
   - α螺旋识别
   - β折叠识别
   - DSSP算法分析（需安装DSSP）

3. **calculate_structure_metrics** - 结构几何指标
   - 质心计算
   - 结构尺寸
   - 原子坐标统计

4. **extract_sequence** - 序列提取
   - 从PDB提取蛋白质序列
   - 多链支持

5. **predict_structure_esm** - ESM-Fold结构预测
   - 使用Meta的ESM-Fold模型
   - 从序列直接预测3D结构

## 安装

```bash
# 基础安装
pip install -e .

# 包含ESM结构预测
pip install -e ".[esm]"
```

## 使用方法

### 配置Claude Code

```json
{
  "mcpServers": {
    "bio-structure": {
      "command": "python",
      "args": ["/path/to/bio_studio/mcp-servers/bio-structure-mcp/structure_server.py"]
    }
  }
}
```

### 示例对话

```
用户: 解析这个PDB文件:
[上传1AKI.pdb文件]

AI: [调用 parse_pdb 工具]
    结构信息:
    - 标题: HEN EGG-WHITE LYSOZYME
    - 链数: 1
    - 残基数: 129
    - 原子数: 1001
    - 分辨率: 1.50 Å
    - 实验方法: X-RAY DIFFRACTION
```

```
用户: 分析这个蛋白质的二级结构
[PDB内容]

AI: [调用 analyze_secondary_structure 工具]
    二级结构分析:
    - 序列长度: 129
    - 螺旋: 45.2%
    - 折叠: 18.3%
    - 卷曲: 36.5%
```

```
用户: 预测这个蛋白质的结构: MVHLTPEEKSAVTALWGKVN

AI: [调用 predict_structure_esm 工具]
    使用ESM-Fold v1预测结构...
    返回完整的PDB格式坐标
```

## 技术细节

### 依赖
- **mcp**: MCP协议实现
- **biopython**: PDB解析
- **numpy**: 数值计算
- **fair-esm** (可选): 蛋白质结构预测

### ESM-Fold使用
```python
# 需要GPU（推荐）
pip install fair-esm torch

# CPU模式（较慢）
export CUDA_VISIBLE_DEVICES=""
```

### DSSP安装（可选）
用于高级二级结构分析：
```bash
# Ubuntu/Debian
sudo apt-get install dssp

# macOS
brew install dssp

# 或使用Python版
pip install dssp
```

## 扩展功能

可以添加：
- 配体结合位点分析
- 分子对接结果解析
- 结构比对
- RMSD计算
- B因子分析
- 氢键检测

## 注意事项

- 大型PDB文件可能需要较长处理时间
- ESM-Fold预测需要大量内存（建议8GB+）
- DSSP需要单独安装

## 许可证

MIT
