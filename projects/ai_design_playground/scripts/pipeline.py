#!/usr/bin/env python3
"""ai_design_playground: in-silico toy sequence pipeline."""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TypedDict


def load_config_module():
    cfg_path = Path(__file__).with_name("config.py")
    spec = importlib.util.spec_from_file_location(
        "ai_design_playground_config", cfg_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Failed to load config.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cfg = load_config_module()
VALIDATOR_PATH = Path(__file__).with_name("validate_project.py")


def load_workspace_env():
    script_dir = Path(__file__).resolve().parent
    for candidate_root in script_dir.parents:
        env_path = candidate_root / "lib" / "workspace_env.py"
        if not env_path.exists():
            continue
        spec = importlib.util.spec_from_file_location("bio_workspace_env", env_path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return None


def load_template_runtime():
    script_dir = Path(__file__).resolve().parent
    for candidate_root in script_dir.parents:
        runtime_path = candidate_root / "lib" / "template_runtime.py"
        if not runtime_path.exists():
            continue
        spec = importlib.util.spec_from_file_location("bio_template_runtime", runtime_path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return None


WORKSPACE_ENV = load_workspace_env()
TEMPLATE_RUNTIME = load_template_runtime()


@dataclasses.dataclass(frozen=True)
class FastaRecord:
    id: str
    seq: str
    description: str = ""


class ProteinSummary(TypedDict):
    id: str
    length_aa: int
    invalid_aa: int
    invalid_frac: float


class RNAFoldResult(TypedDict):
    id: str
    length_nt: int
    structure: str | None
    mfe: float | None


def require_tools(tool_names: list[str]) -> None:
    if WORKSPACE_ENV is not None:
        WORKSPACE_ENV.ensure_workspace_path()
    missing = [t for t in tool_names if shutil.which(t) is None]
    if missing:
        raise RuntimeError(f"Missing required tools on PATH: {', '.join(missing)}")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_fasta(records: list[FastaRecord], out_fa: Path, wrap: int = 80) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as f:
        for r in records:
            header = r.id if not r.description else f"{r.id} {r.description}"
            f.write(f">{header}\n")
            s = r.seq.strip().replace(" ", "").replace("\n", "").upper()
            for i in range(0, len(s), wrap):
                f.write(s[i : i + wrap] + "\n")


def parse_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_id: str | None = None
    current_desc = ""
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal current_id, current_desc, seq_parts
        if current_id is None:
            return
        seq = "".join(seq_parts).replace(" ", "").replace("\n", "").upper()
        records.append(FastaRecord(id=current_id, seq=seq, description=current_desc))
        current_id = None
        current_desc = ""
        seq_parts = []

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Invalid FASTA header in {path}")
            parts = header.split(None, 1)
            current_id = parts[0]
            current_desc = parts[1] if len(parts) > 1 else ""
        else:
            if current_id is None:
                raise ValueError(f"FASTA sequence appears before header in {path}")
            seq_parts.append(line)
    flush()
    return records


def gc_frac(dna: str) -> float:
    s = dna.upper()
    if not s:
        return 0.0
    gc = sum(1 for ch in s if ch in ("G", "C"))
    return gc / len(s)


def generate_toy_dna(length_bp: int, gc_target: float, rng: random.Random) -> str:
    p_gc = min(0.95, max(0.05, gc_target))
    bases_gc = ("G", "C")
    bases_at = ("A", "T")
    seq: list[str] = []
    for _ in range(length_bp):
        if rng.random() < p_gc:
            idx = rng.randrange(2)
            seq.append(bases_gc[idx])
        else:
            idx = rng.randrange(2)
            seq.append(bases_at[idx])
    return "".join(seq)


def contains_any_site(dna: str, sites: dict[str, str]) -> list[str]:
    upper = dna.upper()
    hits: list[str] = []
    for name, site in sites.items():
        if site.upper() in upper:
            hits.append(name)
    return hits


def generate_toy_records(
    n: int, length_bp: int, gc_target: float, seed: int
) -> list[FastaRecord]:
    rng = random.Random(seed)
    records: list[FastaRecord] = []
    attempts = 0
    while len(records) < n:
        attempts += 1
        if attempts > n * 50:
            raise RuntimeError("Too many failed attempts generating toy DNA")
        dna = generate_toy_dna(length_bp, gc_target, rng)
        if contains_any_site(dna, cfg.AVOID_RESTRICTION_SITES):
            continue
        record_id = f"toy_{len(records) + 1:03d}"
        records.append(
            FastaRecord(
                id=record_id, seq=dna, description=f"len={length_bp} gc~{gc_target:.2f}"
            )
        )
    return records


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True, check=True
    )


def run_seqkit_stats(fasta: Path, out_txt: Path) -> None:
    cp = run(["seqkit", "stats", str(fasta)])
    write_text(out_txt, cp.stdout)


def run_prodigal(dna_fa: Path, out_dir: Path, mode: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    proteins = out_dir / "prodigal_proteins.faa"
    genes = out_dir / "prodigal_genes.fna"
    out = out_dir / "prodigal.out"
    cmd = [
        "prodigal",
        "-i",
        str(dna_fa),
        "-a",
        str(proteins),
        "-d",
        str(genes),
        "-o",
        str(out),
        "-f",
        "gff",
    ]
    if mode:
        cmd.extend(["-p", mode])
    run(cmd)
    return {
        "proteins": str(proteins),
        "genes": str(genes),
        "out": str(out),
    }


def pick_top_proteins(proteins_faa: Path, top_k: int) -> list[FastaRecord]:
    proteins = parse_fasta(proteins_faa)
    proteins.sort(key=lambda r: len(r.seq), reverse=True)
    return proteins[:top_k]


def protein_basic_metrics(seq: str) -> dict[str, int | float]:
    aa = seq.replace("*", "").upper()
    length = len(aa)
    if length == 0:
        return {"length_aa": 0, "invalid_aa": 0, "invalid_frac": 0.0}
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = sum(1 for ch in aa if ch not in valid)
    return {
        "length_aa": length,
        "invalid_aa": invalid,
        "invalid_frac": round(invalid / length, 6),
    }


def esm_contacts_score(
    seqs: list[FastaRecord],
    model_name: str,
) -> dict[str, dict[str, float | None]]:
    import importlib

    torch = importlib.import_module("torch")
    esm = importlib.import_module("esm")

    load_fn = getattr(esm.pretrained, "load_model_and_alphabet_hub", None)
    if load_fn is None:
        load_fn = getattr(esm.pretrained, "load_model_and_alphabet", None)
    if load_fn is None:
        raise RuntimeError("esm.pretrained loader function not found")
    model, alphabet = load_fn(model_name)
    model = model.eval()
    if torch.cuda.is_available():
        model = model.to("cuda")
    batch_converter = alphabet.get_batch_converter()

    results: dict[str, dict[str, float | None]] = {}
    for rec in seqs:
        raw = rec.seq.upper().replace(" ", "").replace("\n", "")
        raw = raw.replace("*", "")
        valid = set("ACDEFGHIKLMNPQRSTVWY")
        clean = "".join(ch if ch in valid else "X" for ch in raw)
        if not clean:
            results[rec.id] = {"esm_contacts_mean": None}
            continue

        data = [(rec.id, clean)]
        _, _, tokens = batch_converter(data)
        if torch.cuda.is_available():
            tokens = tokens.to("cuda")
        with torch.no_grad():
            out = model(tokens, return_contacts=True)
        contacts = out.get("contacts")
        if contacts is None:
            results[rec.id] = {"esm_contacts_mean": None}
            continue
        c = contacts[0]
        results[rec.id] = {"esm_contacts_mean": round(float(c.mean().item()), 6)}
    return results


def rnafold_mfe(rna_seq: str) -> dict[str, float | str | None]:
    cp = subprocess.run(
        ["RNAfold", "--noPS"],
        input=f">seq\n{rna_seq}\n",
        text=True,
        capture_output=True,
        check=True,
    )
    lines = [ln.strip() for ln in cp.stdout.splitlines() if ln.strip()]
    if not lines:
        return {"structure": None, "mfe": None}

    best_parts: list[str] | None = None
    for ln in lines:
        parts = ln.split()
        if not parts:
            continue
        last = parts[-1]
        if last.startswith("(") and last.endswith(")"):
            best_parts = parts

    if best_parts is None:
        return {"structure": None, "mfe": None}

    structure = best_parts[0]
    mfe = None
    try:
        mfe = float(best_parts[-1].strip("()"))
    except ValueError:
        mfe = None
    return {"structure": structure, "mfe": mfe}


def run_shared_ai_design_analysis(
    dna_fa: Path,
    dna_records: list[FastaRecord],
) -> dict[str, object] | None:
    if TEMPLATE_RUNTIME is None:
        return None
    rna_sequences = [
        (rec.id, rec.seq.replace("T", "U"))
        for rec in dna_records[: min(3, len(dna_records))]
    ]
    return TEMPLATE_RUNTIME.run_ai_design_playground_analysis(
        cfg,
        dna_fa,
        rna_sequences,
    )


def run_validation() -> int:
    spec = importlib.util.spec_from_file_location("project_validation", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载验证脚本: {VALIDATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_validation()


def print_steps() -> None:
    print("可用步骤:")
    print("  generate_toy_sequences: 生成 toy DNA 并写入 data/raw/toy_dna.fa")
    print("  run_local_analysis: 运行 seqkit / prodigal / RNAfold 等本地分析")
    print("  write_reports: 输出 report.json 和 report.md")


def main() -> None:
    parser = argparse.ArgumentParser(description="In-silico toy sequence pipeline")
    parser.add_argument(
        "--steps", action="store_true", help="只列出可用步骤，不执行"
    )
    parser.add_argument(
        "--validate", action="store_true", help="运行项目级自检并退出"
    )
    _ = parser.add_argument(
        "--n", type=int, default=cfg.N_SEQUENCES, help="Number of toy DNA sequences"
    )
    _ = parser.add_argument(
        "--length", type=int, default=cfg.DNA_LENGTH_BP, help="Toy DNA length (bp)"
    )
    _ = parser.add_argument(
        "--gc", type=float, default=cfg.GC_TARGET, help="Target GC fraction (0-1)"
    )
    _ = parser.add_argument(
        "--seed", type=int, default=cfg.RANDOM_SEED, help="RNG seed"
    )
    parser.add_argument(
        "--top-proteins",
        type=int,
        default=5,
        help="How many predicted proteins to summarize",
    )
    _ = parser.add_argument(
        "--use-esm-contacts",
        action="store_true",
        default=cfg.USE_ESM_CONTACTS_DEFAULT,
        help="Run ESM contact prediction for top proteins (downloads model if missing)",
    )
    _ = parser.add_argument(
        "--esm-model", default=cfg.ESM_MODEL_NAME_DEFAULT, help="ESM model name"
    )
    args = parser.parse_args()

    if args.steps:
        print_steps()
        return

    if args.validate:
        exit_code = run_validation()
        report_path = Path(cfg.LOGS_DIR) / "validation_report.json"
        print(f"项目验证完成，请查看: {report_path}")
        sys.exit(exit_code)

    require_tools(["seqkit", "prodigal", "RNAfold"])

    raw_dir = Path(cfg.RAW_DIR)
    processed_dir = Path(cfg.PROCESSED_DIR)
    results_dir = Path(cfg.RESULTS_DIR)
    for d in (raw_dir, processed_dir, results_dir):
        d.mkdir(parents=True, exist_ok=True)

    dna_fa = raw_dir / "toy_dna.fa"
    seqkit_txt = results_dir / "seqkit_stats.txt"
    report_json = results_dir / "report.json"
    report_md = results_dir / "report.md"

    dna_records = generate_toy_records(args.n, args.length, args.gc, args.seed)
    write_fasta(dna_records, dna_fa)
    shared_analysis = run_shared_ai_design_analysis(dna_fa, dna_records)
    if shared_analysis is None:
        run_seqkit_stats(dna_fa, seqkit_txt)
        prodigal_out = run_prodigal(dna_fa, processed_dir, cfg.PRODIGAL_MODE)
        rnafold_results: list[RNAFoldResult] = []
        for rec in dna_records[: min(3, len(dna_records))]:
            rna = rec.seq.replace("T", "U")
            fold = rnafold_mfe(rna)
            structure = fold.get("structure")
            if not isinstance(structure, str):
                structure = None
            mfe = fold.get("mfe")
            if not isinstance(mfe, float):
                mfe = None
            rnafold_results.append(
                {
                    "id": rec.id,
                    "length_nt": len(rna),
                    "structure": structure,
                    "mfe": mfe,
                }
            )
    else:
        seqkit_txt = Path(shared_analysis["seqkit_stats"])
        prodigal_out = shared_analysis["prodigal"]
        rnafold_results = shared_analysis["rnafold"]

    proteins_faa = Path(prodigal_out["proteins"])
    top_proteins = pick_top_proteins(proteins_faa, args.top_proteins)

    protein_summaries: list[ProteinSummary] = []
    for rec in top_proteins:
        seq = rec.seq.replace("*", "")
        metrics = protein_basic_metrics(seq)
        protein_summaries.append(
            {
                "id": rec.id,
                "length_aa": int(metrics["length_aa"]),
                "invalid_aa": int(metrics["invalid_aa"]),
                "invalid_frac": float(metrics["invalid_frac"]),
            }
        )

    esm_scores: dict[str, dict[str, float | None]] = {}
    if args.use_esm_contacts and top_proteins:
        esm_scores = esm_contacts_score(top_proteins, args.esm_model)

    gc_values = [gc_frac(r.seq) for r in dna_records]
    mean_gc = sum(gc_values) / len(gc_values)
    stdev_gc = math.sqrt(sum((x - mean_gc) ** 2 for x in gc_values) / len(gc_values))

    summary = {
        "safety_notice": cfg.SAFETY_NOTICE,
        "generated": {
            "n": args.n,
            "length_bp": args.length,
            "gc_target": args.gc,
            "seed": args.seed,
            "output_fasta": str(dna_fa),
            "observed_gc_mean": round(mean_gc, 6),
            "observed_gc_stdev": round(stdev_gc, 6),
        },
        "seqkit_stats": str(seqkit_txt),
        "prodigal": prodigal_out,
        "top_proteins": protein_summaries,
        "esm_contacts": esm_scores,
        "rnafold": rnafold_results,
    }

    write_text(report_json, json.dumps(summary, indent=2, ensure_ascii=False))

    md_lines: list[str] = []
    md_lines.append("# AI Design Playground (in-silico)\n")
    md_lines.append("This generates toy DNA and runs local analyses.\n")
    md_lines.append(f"\n**Safety**: {cfg.SAFETY_NOTICE}\n")
    md_lines.append("\n## Generated Toy DNA\n")
    md_lines.append(f"- n: {args.n}\n")
    md_lines.append(f"- length: {args.length} bp\n")
    md_lines.append(f"- target GC: {args.gc:.2f}\n")
    md_lines.append(f"- observed GC mean: {mean_gc:.3f}\n")
    md_lines.append(f"- observed GC stdev: {stdev_gc:.3f}\n")
    md_lines.append(f"- fasta: `{dna_fa}`\n")
    md_lines.append("\n## ORF Calling (Prodigal)\n")
    md_lines.append(f"- proteins: `{prodigal_out['proteins']}`\n")
    md_lines.append(f"- genes: `{prodigal_out['genes']}`\n")
    md_lines.append("\n## Top Predicted Proteins\n")
    for p in protein_summaries:
        extra = ""
        if (
            p["id"] in esm_scores
            and esm_scores[p["id"]].get("esm_contacts_mean") is not None
        ):
            extra = f" | esm_contacts_mean={esm_scores[p['id']]['esm_contacts_mean']}"
        md_lines.append(
            f"- {p['id']}: len={p.get('length_aa')} invalid_frac={p.get('invalid_frac')}{extra}\n"
        )
    md_lines.append("\n## RNAfold (first 3 toy sequences as RNA)\n")
    for r in rnafold_results:
        md_lines.append(f"- {r['id']}: mfe={r.get('mfe')}\n")
    md_lines.append("\n## Outputs\n")
    md_lines.append(f"- report: `{report_json}`\n")
    md_lines.append(f"- report: `{report_md}`\n")

    write_text(report_md, "".join(md_lines))

    print(f"OK: wrote {report_md}")
    print(f"OK: wrote {report_json}")


if __name__ == "__main__":
    main()
