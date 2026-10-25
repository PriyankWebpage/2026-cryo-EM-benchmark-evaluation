#!/usr/bin/env python3
"""Grade every predicted model against its answer key.

One row per model file and reference (the deposited model; also the symmetry assembly where
there is one).

US-align (after superposition):
    TM_ref        TM-score normalised by the reference, all chains taken together and no
                  assumption about chain or fragment order (-mm 5 -ter 1 -split 0), so a model
                  built in fragments is not penalised for that
    TM_ref_mm1    TM-score with whole chains paired one to one (-mm 1 -ter 0), for comparison
    aligned_pct   aligned residues in percent of the reference residues
    seq_identity  sequence identity over the aligned residues
Residue scores (in place, no superposition: model and answer key share the map frame). Built and
reference C-alpha atoms are paired one to one by increasing distance, up to 3 A:
    recall            paired reference residues / reference residues
    precision_all     paired model residues / all built residues
    precision_in_ref  the same over "visible" built residues only, i.e. those whose sequence
                      position is modelled in the answer key
    completeness      paired residues with the correct amino acid / reference residues
    aa_accuracy       paired residues with the correct amino acid / paired residues
    ca_rmsd, bb_rmsd  RMSD over the paired residues (C-alpha; backbone N, CA, C, O)
    F1_1A/2A/3A       F1 of recall and precision with the pairing cutoff at 1, 2 and 3 A
    breaks_per_100    consecutive C-alpha atoms more than 4 A apart, per 100 built residues
A missing or empty model is a failure: TM-score, recall, precision, completeness and F1 are 0,
RMSD is left blank.
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np

from common import (MODEL_FILES, PAIR_CUTOFF, TsvWriter, add_target_arguments, backbone_rmsd,
                    chain_breaks, chain_indices, f1, greedy_pairs, load_targets, model_path,
                    positions_in_reference, reference_paths, residues)

COLUMNS = ["pdb", "tool", "model", "reference", "final", "failure", "ref_res", "built_res",
           "built_in_ref", "fragments", "TM_ref", "TM_ref_mm1", "aligned_pct", "seq_identity",
           "recall", "precision_all", "precision_in_ref", "completeness", "aa_accuracy",
           "ca_rmsd", "bb_rmsd", "F1_1A", "F1_2A", "F1_3A", "breaks_per_100", "file"]
ZERO_SCORES = dict(built_in_ref=0, TM_ref=0.0, TM_ref_mm1=0.0, aligned_pct=0.0, recall=0.0,
                   precision_all=0.0, precision_in_ref=0.0, completeness=0.0,
                   F1_1A=0.0, F1_2A=0.0, F1_3A=0.0)


def usalign(exe, model, ref, options):
    """(TM-score normalised by the reference, identity over the alignment, aligned length)."""
    cmd = [exe, str(model), str(ref)] + options + ["-mol", "prot", "-outfmt", "2"]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=3600).stdout
        fields = [line for line in out.splitlines() if line and not line.startswith("#")][0]
        fields = fields.split("\t")
        return float(fields[3]), float(fields[7]), int(fields[10])
    except (subprocess.TimeoutExpired, IndexError, ValueError):
        return float("nan"), float("nan"), 0


def grade(model_file, ref_file, usalign_exe):
    """Scores of one model against one reference, as a dict of formatted values."""
    ref = residues(ref_file)
    row = {"ref_res": len(ref)}
    if not model_file.exists() or model_file.stat().st_size == 0:
        row.update(failure="yes", built_res=0, fragments=0, **ZERO_SCORES)
        return row
    mod = residues(model_file)
    row.update(failure="no" if mod else "yes", built_res=len(mod),
               fragments=len(chain_indices(mod)))
    if not mod:
        row.update(ZERO_SCORES)
        return row

    tm5, identity, aligned = usalign(usalign_exe, model_file, ref_file,
                                     ["-mm", "5", "-ter", "1", "-split", "0"])
    tm1, _, _ = usalign(usalign_exe, model_file, ref_file, ["-mm", "1", "-ter", "0"])
    row.update(TM_ref=f"{tm5:.4f}", TM_ref_mm1=f"{tm1:.4f}",
               aligned_pct=f"{100 * aligned / len(ref):.1f}", seq_identity=f"{identity:.3f}")

    ref_xyz = np.array([r["ca"] for r in ref])
    mod_xyz = np.array([r["ca"] for r in mod])
    visible = positions_in_reference(mod, ref)
    row["built_in_ref"] = int(visible.sum())
    pairs = greedy_pairs(ref_xyz, mod_xyz, PAIR_CUTOFF)
    n = len(pairs)
    same_aa = sum(1 for i, j, _ in pairs if ref[i]["name1"] == mod[j]["name1"])
    n_visible = sum(1 for _, j, _ in pairs if visible[j])
    row["recall"] = f"{100 * n / len(ref_xyz):.1f}"
    row["precision_all"] = f"{100 * n / len(mod_xyz):.1f}"
    row["precision_in_ref"] = f"{100 * n_visible / visible.sum():.1f}" if visible.any() else ""
    row["completeness"] = f"{100 * same_aa / len(ref_xyz):.1f}"
    row["aa_accuracy"] = f"{100 * same_aa / n:.1f}" if n else ""
    if n:
        row["ca_rmsd"] = f"{np.sqrt(np.mean([d * d for _, _, d in pairs])):.2f}"
        bb = backbone_rmsd(pairs, ref, mod)
        row["bb_rmsd"] = f"{bb:.2f}" if bb is not None else ""
    for cutoff in (1.0, 2.0, 3.0):
        row[f"F1_{int(cutoff)}A"] = f"{f1(ref_xyz, mod_xyz, cutoff):.1f}"
    row["breaks_per_100"] = f"{100 * chain_breaks(mod) / len(mod):.2f}"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="data directory (see README)")
    parser.add_argument("--usalign", default="USalign", help="US-align executable")
    parser.add_argument("--out", type=Path, default=Path("results/grades.tsv"))
    add_target_arguments(parser)
    args = parser.parse_args()

    out = TsvWriter(args.out, COLUMNS)
    n_rows = 0
    for target in load_targets(args.targets_file, args.targets):
        pdb = target["pdb"]
        for tool, label, file_name, final in MODEL_FILES:
            model_file = model_path(args.data, pdb, tool, file_name)
            for ref_name, ref_file in reference_paths(args.data, pdb):
                row = grade(model_file, ref_file, args.usalign)
                row.update(pdb=pdb, tool=tool, model=label, reference=ref_name,
                           final="yes" if final else "no",
                           file=model_file.relative_to(args.data).as_posix())
                out.write(row)
                n_rows += 1
                print(f"{pdb} {tool:12} {label:9} vs {ref_name:9} "
                      f"TM {row.get('TM_ref', ''):>7} recall {row.get('recall', ''):>5} "
                      f"precision {row.get('precision_all', ''):>5}", flush=True)
    out.close()
    print(f"wrote {args.out}: {n_rows} rows")


if __name__ == "__main__":
    main()
