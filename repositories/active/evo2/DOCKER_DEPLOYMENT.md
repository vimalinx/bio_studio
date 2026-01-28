# Evo2 Docker 部署指南

## 硬件信息

- GPU: NVIDIA RTX 5070 (12GB VRAM)
- CUDA: 13.1
- 驱动: 590.48.01

## 推荐模型

| 模型 | 参数量 | 上下文长度 | 显存需求 | 推荐状态 |
|------|--------|-----------|---------|----------|
| evo2_1b_base | 1B | 8K | ~8GB | ✅ 推荐 |
| evo2_7b_base | 7B | 8K | ~16GB | ⚠️ 可能 OOM |
| evo2_7b | 7B | 1M | ~40GB | ❌ 需要多 GPU |
| evo2_40b | 40B | 1M | ~200GB | ❌ 需要多 GPU |

## 快速开始

### 方法 1: 使用部署脚本

```bash
cd /media/vimalinx/Data/bio_studio/evo2

# 使用默认 1B 模型部署
./deploy.sh

# 或指定模型
./deploy.sh evo2_1b_base
```

### 方法 2: 使用 Docker Compose

```bash
cd /media/vimalinx/Data/bio_studio/evo2

# 构建镜像
docker compose build

# 运行容器
docker compose run --rm evo2 bash

# 在容器内运行测试
python -m evo2.test.test_evo2_generation --model_name evo2_1b_base
```

### 方法 3: 使用 Docker 原生命令

```bash
cd /media/vimalinx/Data/bio_studio/evo2

# 构建镜像
docker build -t evo2 .

# 运行容器
docker run -it --rm \
    --gpus all \
    --shm-size=16g \
    -v ./huggingface_cache:/root/.cache/huggingface \
    -v $(pwd):/workdir \
    evo2 bash
```

## 运行示例代码

进入容器后：

```bash
# 运行测试
python -m evo2.test.test_evo2_generation --model_name evo2_1b_base

# 运行示例脚本
python examples/run_evo2.py
```

## Python API 使用

```python
from evo2 import Evo2

# 加载模型
model = Evo2('evo2_1b_base')

# 生成 DNA 序列
output = model.generate(
    prompt_seqs=["ACGT"],
    n_tokens=100,
    temperature=1.0,
    top_k=4
)
print(output.sequences[0])
```

## 目录结构

```
evo2/
├── deploy.sh              # 部署脚本
├── docker-compose.yml     # Docker Compose 配置
├── Dockerfile             # Docker 镜像定义
├── examples/              # 示例代码
│   └── run_evo2.py       # 使用示例
└── huggingface_cache/    # 模型缓存目录 (自动创建)
```

## 注意事项

1. **首次运行**: 模型会从 HuggingFace 下载，需要一些时间
2. **显存不足**: 如果遇到 OOM，尝试使用更小的模型或减少批次大小
3. **FP8 支持**: 40B 和 1B 模型需要 FP8 支持以获得最佳精度
