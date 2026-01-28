#!/usr/bin/env python3
"""
Demo: Using Biomni tools directly as a library
"""

import sys
import subprocess
from pathlib import Path

# 1. 动态定位并添加 Biomni
try:
    root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
    BIOMNI_PATH = Path(root) / "repositories/active/Biomni"
    sys.path.append(str(BIOMNI_PATH))
    print(f"✅ Loaded Biomni from: {BIOMNI_PATH}")
except Exception as e:
    print(f"❌ Failed to locate Biomni: {e}")
    sys.exit(1)

# 2. 导入工具
try:
    from biomni.tool.molecular_biology import (
        annotate_open_reading_frames,
        find_restriction_sites
    )
except ImportError as e:
    print(f"❌ ImportError: {e}")
    print("Ensure you have 'biopython', 'pandas', 'beautifulsoup4' installed.")
    sys.exit(1)

# Test Sequence
# Contains:
# - Start Codon (ATG)
# - Some coding sequence
# - EcoRI site (GAATTC) at the end
dna_seq = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGAATTC"

print(f"\n🧬 Analyzing Sequence (Length: {len(dna_seq)} bp):")
print(f"  {dna_seq}")

# 3. Find ORFs
print("\n[1] Finding ORFs (min_length=30)...")
try:
    orfs_result = annotate_open_reading_frames(dna_seq, min_length=30)
    orfs = orfs_result.get('orfs', [])
    print(f"  Found {len(orfs)} ORFs.")
    for i, orf in enumerate(orfs):
        print(f"  - ORF {i+1} ({orf.strand}): {orf.start}-{orf.end} | AA: {orf.aa_sequence}")
except Exception as e:
    print(f"  ❌ Error finding ORFs: {e}")

# 4. Find Restriction Sites
print("\n[2] Finding Restriction Sites (EcoRI)...")
try:
    # EcoRI recognition: G^AATTC
    sites_result = find_restriction_sites(dna_seq, ["EcoRI"], is_circular=False)
    ecori_info = sites_result["restriction_sites"].get("EcoRI")
    
    if ecori_info and ecori_info['sites']:
        print(f"  ✅ EcoRI found at positions: {ecori_info['sites']}")
        print(f"  Recognition Sequence: {ecori_info['recognition_sequence']}")
    else:
        print("  No EcoRI sites found.")
except Exception as e:
    print(f"  ❌ Error finding restriction sites: {e}")
