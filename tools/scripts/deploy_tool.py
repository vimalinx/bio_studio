#!/usr/bin/env python3
"""
Bio Studio - 工具部署助手
自动检测已部署的工具并生成对应的Skill文档
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime


class ToolDeployer:
    """工具部署和Skill生成助手"""

    def __init__(self):
        self.bio_studio = Path.home() / "bio_studio"
        self.repos_dir = self.bio_studio / "repositories" / "active"
        self.skills_dir = self.bio_studio / ".claude" / "skills"
        self.tools_index = self.bio_studio / "tools" / "external" / "index.json"

    def scan_repositories(self):
        """扫描所有已部署的仓库"""
        repos = []
        for repo_path in self.repos_dir.iterdir():
            if repo_path.is_dir() and not repo_path.name.startswith('.'):
                repo_info = self.analyze_repository(repo_path)
                if repo_info:
                    repos.append(repo_info)
        return repos

    def analyze_repository(self, repo_path):
        """分析仓库信息"""
        repo_name = repo_path.name
        info = {
            "name": repo_name,
            "path": str(repo_path),
            "description": "",
            "commands": [],
            "dependencies": [],
            "docs": [],
            "has_skill": False,
            "skill_name": ""
        }

        # 检查是否已有Skill
        skill_path = self.skills_dir / repo_name.replace("-", "_")
        info["has_skill"] = skill_path.exists()
        info["skill_name"] = skill_path.name if info["has_skill"] else ""

        # 查找README
        for readme in ["README.md", "readme.md", "README.rst"]:
            readme_file = repo_path / readme
            if readme_file.exists():
                info["docs"].append(str(readme_file))
                # 提取描述
                try:
                    with open(readme_file) as f:
                        lines = f.readlines()[:10]  # 读前10行
                        for line in lines:
                            if line.strip() and not line.startswith('#'):
                                info["description"] = line.strip()
                                break
                except:
                    pass
                break

        # 查找可执行文件和脚本
        for item in repo_path.rglob("*"):
            if item.is_file():
                # Python脚本
                if item.suffix == '.py' and item.name not in ['__init__.py', 'setup.py']:
                    info["commands"].append({
                        "type": "python",
                        "name": item.stem,
                        "path": str(item),
                        "run": f"python {item.relative_to(repo_path)}"
                    })
                # Shell脚本
                elif item.suffix in ['.sh', '.bash']:
                    info["commands"].append({
                        "type": "shell",
                        "name": item.stem,
                        "path": str(item),
                        "run": f"bash {item.relative_to(repo_path)}"
                    })
                # 可执行文件
                elif os.access(item, os.X_OK):
                    info["commands"].append({
                        "type": "executable",
                        "name": item.name,
                        "path": str(item),
                        "run": f"./{item.relative_to(repo_path)}"
                    })

        # 查找requirements.txt或environment.yml
        if (repo_path / "requirements.txt").exists():
            info["dependencies"].append("requirements.txt")
        if (repo_path / "environment.yml").exists():
            info["dependencies"].append("environment.yml")
        if (repo_path / "setup.py").exists():
            info["dependencies"].append("setup.py")

        return info if info["commands"] or info["docs"] else None

    def generate_skill_template(self, repo_info):
        """生成Skill模板"""
        skill_name = repo_info["name"].replace("-", "_")
        skill_dir = self.skills_dir / skill_name

        # 创建Skill目录
        skill_dir.mkdir(exist_ok=True)

        # 生成SKILL.md
        skill_content = f"""---
name: {skill_name}
description: {repo_info.get('description', repo_info['name'])}
tools:
  - Bash
  - Read
  - Write
---

# {skill_name} Skill

## 工具描述

{repo_info.get('description', '自动化工具: ' + repo_info['name'])}

## 安装位置

```
{repo_info['path']}
```

## 可用命令

"""

        # 添加命令
        if repo_info["commands"]:
            for cmd in repo_info["commands"][:10]:  # 最多显示10个命令
                skill_content += f"""
### `{cmd['run']}`

**类型**: {cmd['type']}
**路径**: `{cmd['path']}`
"""
                if cmd['type'] == 'python':
                    skill_content += f"""
**使用示例**:
```bash
cd {repo_info['path']}
{cmd['run']} --help
```
"""

        # 添加依赖信息
        if repo_info["dependencies"]:
            skill_content += f"""
## 依赖安装

"""
            for dep in repo_info["dependencies"]:
                dep_file = Path(repo_info['path']) / dep
                skill_content += f"""
### {dep}

```bash
cd {repo_info['path']}
"""
                if dep == "requirements.txt":
                    skill_content += "pip install -r requirements.txt\n"
                elif dep == "environment.yml":
                    skill_content += "conda env create -f environment.yml\n"
                elif dep == "setup.py":
                    skill_content += "pip install -e .\n"
                skill_content += "```\n"

        # 添加文档链接
        if repo_info["docs"]:
            skill_content += f"""
## 文档

"""
            for doc in repo_info["docs"]:
                skill_content += f"- `{doc}`\n"

        # 添加使用说明
        skill_content += f"""
## 典型使用场景

### 场景1: 基本使用

**用户**: "使用 {skill_name} 分析这个文件"

**AI应该**:
1. 检查输入文件格式
2. 选择合适的命令: `在此处填写推荐命令`
3. 运行分析
4. 解读结果

### 场景2: 高级选项

**用户**: "使用 {skill_name} 的特定参数"

**AI应该**:
1. 查看 `--help` 了解所有选项
2. 根据需求选择参数
3. 执行并验证结果

## 注意事项

- ⚠️ 首次使用需要安装依赖
- ⚠️ 检查输入文件格式要求
- ⚠️ 注意输出文件位置

## 最后更新

**时间**: {datetime.now().strftime('%Y-%m-%d %H:%M')}
**仓库**: {repo_info['name']}
**状态**: 已部署
"""

        # 写入SKILL.md
        skill_file = skill_dir / "SKILL.md"
        with open(skill_file, 'w') as f:
            f.write(skill_content)

        return skill_file

    def update_index(self, repos):
        """更新工具索引"""
        index = {
            "last_updated": datetime.now().isoformat(),
            "total_tools": len(repos),
            "tools": []
        }

        for repo in repos:
            index["tools"].append({
                "name": repo["name"],
                "path": repo["path"],
                "has_skill": repo["has_skill"],
                "skill_name": repo["skill_name"],
                "description": repo.get("description", ""),
                "commands_count": len(repo["commands"])
            })

        # 保存索引
        self.tools_index.parent.mkdir(exist_ok=True)
        with open(self.tools_index, 'w') as f:
            json.dump(index, f, indent=2)

        return index

    def deploy_interactive(self, repo_path=None):
        """交互式部署"""
        if repo_path:
            # 部署单个仓库
            repo_info = self.analyze_repository(Path(repo_path))
            if not repo_info:
                print("❌ 无法分析仓库")
                return

            print(f"\n🔍 检测到工具: {repo_info['name']}")
            print(f"📍 位置: {repo_info['path']}")
            print(f"📝 描述: {repo_info.get('description', '无')}")
            print(f"🔧 命令数: {len(repo_info['commands'])}")

            if repo_info["has_skill"]:
                print(f"✅ 已有Skill: {repo_info['skill_name']}")
                choice = input("\n是否更新Skill? (y/n): ").lower()
                if choice != 'y':
                    return
            else:
                print(f"⚠️  尚未创建Skill")
                choice = input("\n是否创建Skill? (y/n): ").lower()
                if choice != 'y':
                    return

            # 生成Skill
            skill_file = self.generate_skill_template(repo_info)
            print(f"✅ Skill已生成: {skill_file}")
            print(f"📝 请编辑 {skill_file} 添加详细使用说明")

        else:
            # 扫描所有仓库
            print("\n🔍 扫描已部署的仓库...")
            repos = self.scan_repositories()

            if not repos:
                print("❌ 未找到任何仓库")
                return

            print(f"\n✅ 找到 {len(repos)} 个仓库:\n")

            for i, repo in enumerate(repos, 1):
                status = "✅" if repo["has_skill"] else "⚠️ "
                print(f"{i}. {status} {repo['name']}")
                print(f"   路径: {repo['path']}")
                print(f"   命令: {len(repo['commands'])} 个")
                print(f"   描述: {repo.get('description', '无')[:50]}")
                print()

            # 更新索引
            self.update_index(repos)
            print(f"✅ 索引已更新: {self.tools_index}")

            # 询问是否为没有Skill的工具创建
            no_skill = [r for r in repos if not r["has_skill"]]
            if no_skill:
                choice = input(f"\n为 {len(no_skill)} 个无Skill的工具创建模板? (y/n): ").lower()
                if choice == 'y':
                    for repo in no_skill:
                        skill_file = self.generate_skill_template(repo)
                        print(f"✅ {repo['name']}: {skill_file}")


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="Bio Studio 工具部署助手")
    parser.add_argument("--repo", "-r", help="指定仓库路径")
    parser.add_argument("--scan", "-s", action="store_true", help="扫描所有仓库")

    args = parser.parse_args()

    deployer = ToolDeployer()

    if args.repo:
        deployer.deploy_interactive(args.repo)
    elif args.scan:
        deployer.deploy_interactive()
    else:
        # 默认交互模式
        print("🧬 Bio Studio - 工具部署助手\n")
        print("1. 扫描所有仓库并生成Skill")
        print("2. 为特定仓库生成Skill")
        print("3. 仅更新索引")

        choice = input("\n选择操作 (1/2/3): ").strip()

        if choice == "1":
            deployer.deploy_interactive()
        elif choice == "2":
            repo = input("输入仓库路径: ").strip()
            deployer.deploy_interactive(repo)
        elif choice == "3":
            repos = deployer.scan_repositories()
            deployer.update_index(repos)
            print("✅ 索引已更新")
        else:
            print("无效选择")


if __name__ == "__main__":
    main()
