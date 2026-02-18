from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

# ----------------------------------------------------------------------
# FASTA I/O
# ----------------------------------------------------------------------


def read_fasta(path: str) -> Tuple[str, str]:
    
    header = None
    seq_lines: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                
                continue
            seq_lines.append(line.upper())
    if header is None:
        header = "sequence"
    return header, "".join(seq_lines)


def write_fasta(path: str, header: str, seq: str, width: int = 60) -> None:
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


# ----------------------------------------------------------------------
# Design parser
# ----------------------------------------------------------------------


@dataclass
class PlasmidDesign:
    enzymes: List[str]
    antibiotics: List[str]


def parse_design(path: str) -> PlasmidDesign:
   
    enzymes: List[str] = []
    antibiotics: List[str] = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [x.strip() for x in line.split(",", 1)]
            if len(parts) != 2:
                print(f"[WARN] Could not parse line in design file: {line}")
                continue
            left, right = parts
            if left.startswith("Multiple_Cloning_Site"):
                enzymes.append(right)
            elif left.startswith("Antibiotic_marker"):
                antibiotics.append(right)
            else:
                print(f"[WARN] Unknown design directive '{left}', skipping.")
    return PlasmidDesign(enzymes=enzymes, antibiotics=antibiotics)


# ----------------------------------------------------------------------
# Markers loader
# ----------------------------------------------------------------------


def load_markers(path: str) -> Dict[str, str]:
    
    markers: Dict[str, str] = {}
    with open(path) as fh:
        header = fh.readline().strip()
        if not header:
            print("[WARN] markers.tab is empty.")
            return markers
        cols = header.split("\t")
        col_index = {name: i for i, name in enumerate(cols)}
        if "antibiotic_name" not in col_index or "sequence" not in col_index:
            print("[WARN] markers.tab missing required columns; "
                  "expected 'antibiotic_name' and 'sequence'.")
            return markers

        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(cols):
                print(f"[WARN] Skipping malformed markers line: {line}")
                continue
            ab_name = parts[col_index["antibiotic_name"]].strip()
            seq = parts[col_index["sequence"]].strip().upper()
            if not ab_name or not seq:
                continue
            markers[ab_name] = seq
    return markers


# ----------------------------------------------------------------------
# ORI finder
# ----------------------------------------------------------------------


def find_ori(seq: str, window: int = 800, step: int = 100) -> Tuple[int, int]:
    
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0, 0
    if window > n:
        window = n

    best_start = 0
    best_score = -1.0

    for start in range(0, max(1, n - window + 1), step):
        win = seq[start:start + window]
        if not win:
            continue
        at_count = sum(1 for b in win if b in ("A", "T"))
        at_fraction = at_count / len(win)
        if at_fraction > best_score:
            best_score = at_fraction
            best_start = start

    return best_start, min(best_start + window, n)


# ----------------------------------------------------------------------
# Restriction sites and builders
# ----------------------------------------------------------------------

RE_SITE: Dict[str, str] = {
    # Common pUC19 MCS sites
    "EcoRI": "GAATTC",
    "SacI": "GAGCTC",
    "KpnI": "GGTACC",
    "SmaI": "CCCGGG",
    "BamHI": "GGATCC",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "HindIII": "AAGCTT",
    # extend if needed
}


def build_mcs(enzymes: List[str]) -> str:
    
    pieces: List[str] = []
    for name in enzymes:
        if name not in RE_SITE:
            print(f"[WARN] Unknown restriction enzyme '{name}', skipping.")
            continue
        pieces.append(RE_SITE[name])
    return "".join(pieces)


def mutate_motif(seq: str, motif: str) -> str:
    
    s = list(seq)
    m = motif
    mlen = len(motif)
    i = 0
    while i <= len(s) - mlen:
        if "".join(s[i:i + mlen]) == m:
            j = i + 1  # mutate second base
            s[j] = "C" if s[j] == "A" else "A"
            i += mlen
        else:
            i += 1
    return "".join(s)


def build_replication_core() -> str:
   
    # Start codon, 100 As, stop codon
    return "ATG" + "A" * 100 + "TAA"


def build_plasmid(
    input_seq: str,
    design: PlasmidDesign,
    markers: Dict[str, str],
    remove_sites: Optional[List[str]] = None,
) -> str:
    
    ori_start, ori_end = find_ori(input_seq)
    ori_block = input_seq[ori_start:ori_end]

    replication_core = build_replication_core()

    plasmid_parts: List[str] = [ori_block, replication_core]

    # Markers
    for ab_name in design.antibiotics:
        if ab_name not in markers:
            print(f"[WARN] Antibiotic '{ab_name}' not in markers.tab, skipping.")
            continue
        # Add small spacer before each marker
        plasmid_parts.append("AAAA")
        plasmid_parts.append(markers[ab_name])

    # MCS
    mcs_seq = build_mcs(design.enzymes)
    if mcs_seq:
        plasmid_parts.append("AAAA")
        plasmid_parts.append(mcs_seq)

    plasmid_seq = "".join(plasmid_parts)

    # Remove specified restriction sites if requested
    if remove_sites:
        for name in remove_sites:
            if name not in RE_SITE:
                print(f"[WARN] Cannot remove site for unknown enzyme '{name}'.")
                continue
            motif = RE_SITE[name]
            if motif in plasmid_seq:
                plasmid_seq = mutate_motif(plasmid_seq, motif)

    return plasmid_seq


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def infer_sites_to_remove(input_path: str, design: PlasmidDesign) -> List[str]:
    
    lower = input_path.lower()
    sites: List[str] = []
    if "puc19" in lower and "EcoRI" not in design.enzymes:
        sites.append("EcoRI")
    return sites


def main(argv: Optional[List[str]] = None) -> None:
    argv = argv or sys.argv[1:]
    if len(argv) != 4:
        print("Usage: python plasmid_maker.py Input.Fa Design.txt markers.tab Output.Fa")
        raise SystemExit(1)

    input_fa, design_txt, markers_tab, output_fa = argv

    header, in_seq = read_fasta(input_fa)
    design = parse_design(design_txt)
    markers = load_markers(markers_tab)

    remove_sites = infer_sites_to_remove(input_fa, design)

    out_seq = build_plasmid(
        in_seq,
        design=design,
        markers=markers,
        remove_sites=remove_sites,
    )

    write_fasta(output_fa, "designed_plasmid", out_seq)
    print(f"[INFO] Wrote designed plasmid to {output_fa}")


if __name__ == "__main__":
    main()
