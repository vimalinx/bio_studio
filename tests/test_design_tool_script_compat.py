from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DNA_SCRIPT = ROOT / "tools" / "scripts" / "dna_design.py"
MRNA_SCRIPT = ROOT / "tools" / "scripts" / "mrna_optimize.py"


def _install_fake_import_dependencies() -> None:
    bio_module = types.ModuleType("Bio")
    bio_seq_module = types.ModuleType("Bio.Seq")
    bio_sequtils_module = types.ModuleType("Bio.SeqUtils")
    bio_melting_module = types.ModuleType("Bio.SeqUtils.MeltingTemp")
    bio_seqrecord_module = types.ModuleType("Bio.SeqRecord")
    pandas_module = types.ModuleType("pandas")
    numpy_module = types.ModuleType("numpy")
    click_module = types.ModuleType("click")

    class Seq(str):
        pass

    class SeqRecord:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    def gc_fraction(sequence) -> float:
        return 0.5

    def tm_wallace(sequence) -> float:
        return 60.0

    bio_seq_module.Seq = Seq
    bio_sequtils_module.gc_fraction = gc_fraction
    bio_melting_module.Tm_Wallace = tm_wallace
    bio_seqrecord_module.SeqRecord = SeqRecord
    bio_module.SeqIO = types.SimpleNamespace(read=lambda *args, **kwargs: None, write=lambda *args, **kwargs: None)

    sys.modules["Bio"] = bio_module
    sys.modules["Bio.Seq"] = bio_seq_module
    sys.modules["Bio.SeqUtils"] = bio_sequtils_module
    sys.modules["Bio.SeqUtils.MeltingTemp"] = bio_melting_module
    sys.modules["Bio.SeqRecord"] = bio_seqrecord_module
    sys.modules["pandas"] = pandas_module
    sys.modules["numpy"] = numpy_module
    sys.modules["click"] = click_module


def _load_module(path: Path, name: str):
    _install_fake_import_dependencies()
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dna_design_script_imports_without_legacy_gc_symbol() -> None:
    module = _load_module(DNA_SCRIPT, "bio_studio_dna_design_import_test")

    assert hasattr(module, "DNADesigner")


def test_mrna_optimize_script_imports_without_legacy_gc_symbol() -> None:
    module = _load_module(MRNA_SCRIPT, "bio_studio_mrna_optimize_import_test")

    assert hasattr(module, "MRNAOptimizer")
