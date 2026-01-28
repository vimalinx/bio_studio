import subprocess
import sys
import datetime
import os
import socket
from pathlib import Path

# Add conda bin to PATH manually for this script
CONDA_BIN = "/home/vimalinx/miniforge3/envs/bio/bin"
ENV_PYTHON = Path(CONDA_BIN) / "python"
os.environ["PATH"] = f"{CONDA_BIN}:{os.environ.get('PATH', '')}"

def get_cmd_version(cmd, args="--version"):
    try:
        result = subprocess.run([cmd, args], capture_output=True, text=True, timeout=5)
        # 尝试从stdout或stderr获取第一行
        output = result.stdout.strip() if result.stdout.strip() else result.stderr.strip()
        return output.split('\n')[0] if output else "Installed (Version unknown)"
    except FileNotFoundError:
        return "Not Found"
    except Exception as e:
        return f"Error: {e}"

def run_env_python(code):
    python_exec = str(ENV_PYTHON) if ENV_PYTHON.exists() else sys.executable
    return subprocess.run([python_exec, "-c", code], capture_output=True, text=True, timeout=5)

def get_python_version():
    try:
        res = run_env_python("import sys; print(sys.version.split()[0])")
        if res.returncode == 0 and res.stdout.strip():
            return res.stdout.strip().splitlines()[0]
    except Exception:
        pass
    return sys.version.split()[0]

def get_python_pkg_version(package):
    code = (
        "import importlib.util, importlib\n"
        f"pkg = {package!r}\n"
        "spec = importlib.util.find_spec(pkg)\n"
        "if spec is None:\n"
        "    print('Not Installed')\n"
        "else:\n"
        "    module = importlib.import_module(pkg)\n"
        "    ver = getattr(module, '__version__', None)\n"
        "    ver_str = str(ver).strip().lower() if ver is not None else ''\n"
        "    print(ver if ver_str not in ('', 'none') else 'Installed (No __version__)')\n"
    )
    res = run_env_python(code)
    if res.returncode == 0 and res.stdout.strip():
        return res.stdout.strip().splitlines()[0]

    # 回退到pip
    try:
        python_exec = str(ENV_PYTHON) if ENV_PYTHON.exists() else sys.executable
        res = subprocess.run([python_exec, "-m", "pip", "show", package],
                             capture_output=True, text=True, timeout=5)
        for line in res.stdout.split('\n'):
            if line.startswith("Version:"):
                return line.split(":")[1].strip()
        return "Installed (Pip check)"
    except Exception:
        return "Error checking"

def get_hostname():
    output = subprocess.getoutput("hostname").strip()
    if not output or "not found" in output.lower() or "command not found" in output.lower():
        return socket.gethostname()
    return output

def generate_report():
    report = []
    report.append("# Bio Studio Environment Report")
    report.append(f"**Date**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"**Host**: {get_hostname()}")
    report.append("")
    
    # 1. Python Environment
    report.append("## 🐍 Python Environment")
    report.append(f"- **Python**: {get_python_version()}")
    report.append("| Package | Version | Purpose |")
    report.append("|---|---|---|")
    
    py_pkgs = [
        ("Bio", "Bioinformatics (Biopython)"),
        ("torch", "Deep Learning"),
        ("pandas", "Data Analysis"),
        ("numpy", "Numerics"),
        ("scipy", "Scientific Computing"),
        ("matplotlib", "Visualization"),
        ("seaborn", "Visualization"),
        ("sklearn", "Machine Learning"),
        ("pysam", "NGS Data Handling"),
        ("scanpy", "Single Cell"),
        ("esm", "Protein Language Models"),
        ("jupyter", "Interactive Notebooks")
    ]
    
    for pkg, purpose in py_pkgs:
        ver = get_python_pkg_version(pkg)
        report.append(f"| `{pkg}` | {ver} | {purpose} |")
    
    report.append("")
    
    # 2. Bioinformatics Tools
    report.append("## 🧬 Bioinformatics Tools (CLI)")
    report.append("| Tool | Version Detected | Category |")
    report.append("|---|---|---|")
    
    bio_tools = [
        ("fastqc", "--version", "QC"),
        ("multiqc", "--version", "QC Report"),
        ("fastp", "--version", "QC/Trimming"),
        ("seqkit", "version", "Sequence Utils"),
        ("bwa", "", "Alignment (DNA)"), # bwa outputs version to stderr without args usually
        ("bowtie2", "--version", "Alignment (DNA)"),
        ("STAR", "--version", "Alignment (RNA)"),
        ("hisat2", "--version", "Alignment (RNA)"),
        ("samtools", "--version", "BAM Handling"),
        ("bcftools", "--version", "Variant Calling"),
        ("bedtools", "--version", "Genomic Arithmetic"),
        ("featureCounts", "-v", "Quantification"),
        ("muscle", "-version", "MSA"),
        ("mafft", "--version", "MSA"),
        ("iqtree", "--version", "Phylogenetics"),
        ("hmmsearch", "-h", "HMMER (HMM Search)"), 
        ("RNAfold", "--version", "RNA Structure (ViennaRNA)"),
        ("prodigal", "-v", "Gene Prediction")
    ]
    
    for tool, arg, cat in bio_tools:
        # BWA special case
        if tool == "bwa":
            # bwa writes to stderr
            try:
                res = subprocess.run(["bwa"], capture_output=True, text=True)
                # Look for Version: in stderr
                ver = "Installed"
                for line in res.stderr.split('\n'):
                    if "Version:" in line:
                        ver = line.strip()
                        break
            except:
                ver = "Not Found"
        else:
            ver = get_cmd_version(tool, arg)
            
        # Clean up version string (keep it short)
        if len(ver) > 50:
            ver = ver[:47] + "..."
            
        report.append(f"| `{tool}` | {ver} | {cat} |")

    report.append("")
    report.append("## 📂 System Info")
    active_env = subprocess.getoutput("echo $CONDA_DEFAULT_ENV").strip() or "unknown"
    report.append(f"- **Conda Env (active)**: {active_env}")
    if ENV_PYTHON.exists():
        report.append("- **Conda Env (target)**: bio")
    
    import shutil
    total, used, free = shutil.disk_usage(".")
    percent = (used / total) * 100
    report.append(f"- **Disk Usage**: {percent:.1f}% used ({free // (2**30)} GB free)")
    
    return "\n".join(report)

if __name__ == "__main__":
    content = generate_report()
    
    output_path = Path("docs/ENVIRONMENT.md")
    output_path.parent.mkdir(exist_ok=True)
    
    with open(output_path, "w") as f:
        f.write(content)
        
    print(f"✅ Report generated at: {output_path.absolute()}")
    print(content)
