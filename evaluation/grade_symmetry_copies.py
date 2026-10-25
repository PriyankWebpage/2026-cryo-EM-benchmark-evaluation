#!/usr/bin/env python3
"""Per-copy grading for targets whose map holds more copies of the protein than the answer key.

For 26DK (about 47 copies in the map), 9UBC (2) and 9I8C (12) the deposited answer key is a
single copy. Graded in place against that one copy, a tool that built every copy gets a very low
precision, and a tool that built one copy at another position gets completeness 0. Here each
model is graded against the copy it built best, without moving the answer key:

  * Best copy: each model fragment of at least 15 residues is aligned by sequence to the answer
    key and superposed on it (Kabsch fit of the aligned C-alpha atoms). The superposition that
    brings the most model residues within 3 A of the answer key is applied to the whole model.
    Recall, completeness, amino-acid accuracy, RMSD and F1 are then computed as in
    grade_models.py, against that copy.
  * Precision cannot be judged against a one-copy answer key, so the map is used instead: the
    percentage of built C-alpha atoms at or above the recommended EMDB contour level
    (in_density_pct; in_density_visible_pct over the "visible" residues only, i.e. those whose
    sequence position is modelled in the answer key). A model outside the density scores near 0.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from common import (MIN_IDENTITY, MODEL_FILES, PAIR_CUTOFF, PER_COPY_TARGETS, TsvWriter,
                    add_target_arguments, align, backbone_rmsd, chain_indices, f1, greedy_pairs,
                    in_density, load_targets, map_path, model_path, positions_in_reference,
                    read_map, reference_paths, residues, sequence)

COLUMNS = ["pdb", "tool", "ref_res", "built_res", "fragments", "in_density_pct", "built_visible",
           "in_density_visible_pct", "copies_est", "best_fragment", "recall", "completeness",
           "aa_accuracy", "ca_rmsd", "bb_rmsd", "F1_1A", "F1_2A", "F1_3A", "note"]
MIN_FRAGMENT = 15      # residues; shorter fragments are not used to place the copy
MIN_ALIGNED = 10       # aligned residues needed for a superposition
COPY_RADIUS = 6.0      # A; model residues this close to the answer key belong to the graded copy


def kabsch(p, q):
    """Rotation r and translation t with r @ p + t ~ q (least squares)."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, qc - r @ pc


def best_copy(mod, ref, ref_xyz, mod_xyz):
    """(fragment name, rotation, translation) of the superposition matching most residues."""
    ref_seq = sequence(ref)
    best = None
    for name, idx in chain_indices(mod).items():
        if len(idx) < MIN_FRAGMENT:
            continue
        fragment = [mod[i] for i in idx]
        pairs, identity = align(sequence(fragment), ref_seq)
        if len(pairs) < MIN_ALIGNED or identity < MIN_IDENTITY:
            continue
        rot, shift = kabsch(np.array([fragment[i]["ca"] for i, _ in pairs]),
                            np.array([ref[j]["ca"] for _, j in pairs]))
        n = len(greedy_pairs(ref_xyz, mod_xyz @ rot.T + shift, PAIR_CUTOFF))
        if best is None or n > best[0]:
            best = (n, name, rot, shift)
    return None if best is None else best[1:]


def grade(mod, ref, grid, level):
    ref_xyz = np.array([r["ca"] for r in ref])
    mod_xyz = np.array([r["ca"] for r in mod])
    visible = positions_in_reference(mod, ref)
    row = {"in_density_pct": f"{in_density(grid, level, mod_xyz):.1f}",
           "built_visible": int(visible.sum()),
           "in_density_visible_pct": (f"{in_density(grid, level, mod_xyz[visible]):.1f}"
                                      if visible.any() else "0.0"),
           "copies_est": f"{len(mod) / len(ref):.1f}"}
    placed = best_copy(mod, ref, ref_xyz, mod_xyz)
    if placed is None:
        row.update(recall=0.0, completeness=0.0, F1_1A=0.0, F1_2A=0.0, F1_3A=0.0,
                   note="no fragment aligns to the reference")
        return row
    name, rot, shift = placed
    moved = mod_xyz @ rot.T + shift
    pairs = greedy_pairs(ref_xyz, moved, PAIR_CUTOFF)
    same_aa = sum(1 for i, j, _ in pairs if ref[i]["name1"] == mod[j]["name1"])
    row.update(best_fragment=name, recall=f"{100 * len(pairs) / len(ref_xyz):.1f}",
               completeness=f"{100 * same_aa / len(ref_xyz):.1f}")
    if pairs:
        row["aa_accuracy"] = f"{100 * same_aa / len(pairs):.1f}"
        row["ca_rmsd"] = f"{np.sqrt(np.mean([d * d for _, _, d in pairs])):.2f}"
        bb = backbone_rmsd(pairs, ref, mod, transform=(rot, shift))
        row["bb_rmsd"] = f"{bb:.2f}" if bb is not None else ""
    # F1 per copy: precision counted over the model residues on the graded copy only
    on_copy = moved[cKDTree(ref_xyz).query(moved)[0] < COPY_RADIUS]
    for cutoff in (1.0, 2.0, 3.0):
        row[f"F1_{int(cutoff)}A"] = f"{f1(ref_xyz, on_copy, cutoff):.1f}"
    row["note"] = f"{len(on_copy)} model residues on the best copy"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="data directory (see README)")
    parser.add_argument("--out", type=Path, default=Path("results/symmetry_copy_grades.tsv"))
    add_target_arguments(parser)
    args = parser.parse_args()

    out = TsvWriter(args.out, COLUMNS)
    for target in load_targets(args.targets_file, args.targets or PER_COPY_TARGETS):
        pdb = target["pdb"]
        grid = read_map(map_path(args.data, target["emd"]))
        level = float(target["contour_level"])
        ref = residues(reference_paths(args.data, pdb)[0][1])
        for tool, _, file_name, final in MODEL_FILES:
            if not final:
                continue
            path = model_path(args.data, pdb, tool, file_name)
            mod = residues(path) if path.exists() else []
            row = {"pdb": pdb, "tool": tool, "ref_res": len(ref), "built_res": len(mod),
                   "fragments": len(chain_indices(mod))}
            if mod:
                row.update(grade(mod, ref, grid, level))
            else:
                row.update(in_density_pct=0.0, in_density_visible_pct=0.0, recall=0.0,
                           completeness=0.0, F1_1A=0.0, F1_2A=0.0, F1_3A=0.0, note="no model")
            out.write(row)
            print(f"{pdb} {tool:12} in density {row['in_density_pct']:>5} % "
                  f"completeness {row['completeness']:>5} bb RMSD {row.get('bb_rmsd', ''):>5}",
                  flush=True)
    out.close()
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
