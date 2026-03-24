# evo2-7b-hf

这个目录为 Evo2 权重转换后的 HuggingFace 兼容模型目录预留。

## 背景

仓库中包含 `tools/scripts/convert_evo2_weights.py` 和 `tools/scripts/evo2_quantized.py`，说明这里被设计成 Evo2 本地 HF 目录的目标位置。

## 当前状态

- 当前为空，表示该工作面尚未启用或尚未填充权重文件

## 约定

- 如果未来启用，应在这里放置转换后的 `config.json`、权重分片和 tokenizer 相关文件
- 不要把无关模型文件混放到这里
