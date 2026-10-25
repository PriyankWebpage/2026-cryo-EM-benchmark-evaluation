#!/usr/bin/env python3
"""Every table of the report, computed from the grade tables.

Rules:
  * Each tool's answer is its final model; raw / unrefined outputs are reported separately.
  * Headline averages are over the 14 maps common to all tools: 26DK is left out for every tool,
    because MICA's 26DK model lies outside the map. The 15-map averages are reported as well.
  * Precision counts "visible" residues only: built residues whose sequence position is
    modelled in the deposited answer key.
      - normal targets: share of those residues paired within 3 A with an answer-key residue;
      - per-copy targets (26DK, 9UBC, 9I8C: more copies in the map than in the answer key):
        share of those residues whose C-alpha lies in the map density at the EMDB contour.
  * Per-copy targets: completeness, recall, RMSD, F1 and amino-acid accuracy come from the copy
    the tool built best: the in-place values if they are at least as complete, otherwise the
    best-copy values from grade_symmetry_copies.py.
  * A model with fewer than half of its C-alpha atoms in the density is misplaced: completeness,
    recall, precision and F1 are 0; RMSD and amino-acid accuracy are undefined (MICA on 26DK).
"""
import argparse
import json
from pathlib import Path

from common import (MODEL_FILES, PER_COPY_TARGETS, TOOL_NAMES, TOOLS, add_target_arguments,
                    load_targets, read_tsv)

HEADLINE_EXCLUDED = "26DK"
GRADE_KEYS = ("TM_ref", "completeness", "recall", "aa_accuracy", "ca_rmsd", "bb_rmsd", "F1_1A",
              "F1_2A", "F1_3A", "precision_all", "precision_in_ref", "aligned_pct",
              "breaks_per_100", "TM_ref_mm1")
COPY_KEYS = ("completeness", "recall", "aa_accuracy", "ca_rmsd", "bb_rmsd", "F1_1A", "F1_2A",
             "F1_3A")


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def mean(values):
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else None


def final_scores(grade, copy):
    """Scores used in the report for one target and tool; `copy` is the per-copy row or None."""
    row = {k: num(grade[k]) for k in GRADE_KEYS}
    row.update(precision_visible=row["precision_in_ref"], copy=False, misplaced=False)
    if copy is None:
        return row
    row["copy"] = True
    if num(copy["in_density_pct"]) < 50:
        row.update(completeness=0.0, recall=0.0, F1_1A=0.0, F1_2A=0.0, F1_3A=0.0,
                   precision_all=0.0, precision_visible=0.0,
                   aa_accuracy=None, ca_rmsd=None, bb_rmsd=None, misplaced=True)
        return row
    row["precision_visible"] = num(copy["in_density_visible_pct"])
    row["precision_all"] = num(copy["in_density_pct"])
    if num(copy["completeness"]) > (row["completeness"] or 0):
        for k in COPY_KEYS:
            row[k] = num(copy[k])
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--grades", type=Path, default=Path("results/grades.tsv"))
    parser.add_argument("--copies", type=Path, default=Path("results/symmetry_copy_grades.tsv"))
    parser.add_argument("--openstructure", type=Path,
                        default=Path("results/openstructure_scores.tsv"))
    parser.add_argument("--phenix", type=Path, default=Path("results/phenix_scores.tsv"))
    parser.add_argument("--out-json", type=Path, default=Path("results/report_tables.json"))
    parser.add_argument("--out-md", type=Path, default=Path("results/report_tables.md"))
    add_target_arguments(parser)
    args = parser.parse_args()

    order = [t["pdb"] for t in load_targets(args.targets_file, args.targets)]
    grades = read_tsv(args.grades)
    # lookup tables keyed by (target, tool), final models against the deposited answer key
    final_grades = {(r["pdb"], r["tool"]): r for r in grades
                    if r["final"] == "yes" and r["reference"] == "deposited"}
    copy_grades = {(r["pdb"], r["tool"]): r for r in read_tsv(args.copies)}
    ost = {(r["pdb"], r["tool"]): r for r in read_tsv(args.openstructure)
           if r["model"] == "final" and r["reference"] == "deposited"}
    phenix = {(r["pdb"], r["tool"]): r for r in read_tsv(args.phenix)}

    scores = {(p, t): final_scores(final_grades[(p, t)],
                                   copy_grades.get((p, t)) if p in PER_COPY_TARGETS else None)
              for p in order for t in TOOLS}
    target_sets = {"14": [p for p in order if p != HEADLINE_EXCLUDED], "15": order}
    out = {"rules": __doc__, "targets14": target_sets["14"], "tables": {}}
    tables = out["tables"]

    for n, targets in target_sets.items():
        tables[f"headline{n}"] = [
            [TOOL_NAMES[t]] + [mean([scores[(p, t)][k] for p in targets])
                               for k in ("TM_ref", "completeness", "precision_visible", "bb_rmsd")]
            for t in TOOLS]
        tables[f"secondary{n}"] = [
            [TOOL_NAMES[t]] + [mean([scores[(p, t)][k] for p in targets])
                               for k in ("recall", "aa_accuracy", "ca_rmsd", "F1_1A", "F1_3A",
                                         "precision_all", "aligned_pct", "breaks_per_100")]
            for t in TOOLS]
        # QS-score and DockQ judge chain-to-chain interfaces: multi-chain answer keys only
        multi = [p for p in targets if (num(ost[(p, "modelangelo")]["ref_chains"]) or 0) > 1]
        tables[f"ost{n}"] = [
            [TOOL_NAMES[t], mean([num(ost[(p, t)]["lddt"]) for p in targets]),
             mean([num(ost[(p, t)]["bb_lddt"]) for p in targets]),
             mean([num(ost[(p, t)]["qs_global"]) for p in multi]),
             mean([num(ost[(p, t)]["dockq_ave"]) for p in multi])]
            for t in TOOLS]
        out[f"ost_multichain{n}"] = len(multi)
        tables[f"phenix{n}"] = [
            [TOOL_NAMES[t]] + [mean([num(phenix[(p, t)][k]) for p in targets])
                               for k in ("cc_found_pct", "CC_mask", "clashscore",
                                         "molprobity_score", "rama_outliers_pct",
                                         "rama_favored_pct", "rota_outliers_pct")]
            for t in TOOLS]

    # literal grading: in place against the one-copy answer key, precision over all built residues
    tables["literal15"] = [
        [TOOL_NAMES[t]] + [mean([num(final_grades[(p, t)][k]) for p in order])
                           for k in ("TM_ref", "completeness", "precision_all", "bb_rmsd")]
        for t in TOOLS]

    # per-target tables; * marks per-copy values
    formats = {"TM_ref": "{:.3f}", "completeness": "{:.1f}", "precision_visible": "{:.1f}",
               "bb_rmsd": "{:.2f}"}
    for k, fmt in formats.items():
        rows = []
        for p in order:
            cells = [p]
            for t in TOOLS:
                r = scores[(p, t)]
                cell = "n/a" if r[k] is None else fmt.format(r[k])
                if r["misplaced"] and k == "completeness":   # also show the superposed value
                    cell += f" ({float(copy_grades[(p, t)]['completeness']):.1f})"
                if r["copy"] and k != "TM_ref":
                    cell += "*"
                cells.append(cell)
            rows.append(cells)
        tables[f"per_target_{k}"] = rows

    tables["percopy"] = []
    for p in (p for p in PER_COPY_TARGETS if p in order):
        for t in TOOLS:
            c, r = copy_grades[(p, t)], scores[(p, t)]
            if r["misplaced"]:
                completeness = f"0 in place ({float(c['completeness']):.1f} superposed)"
                rmsd = f"n/a ({float(c['bb_rmsd']):.2f} superposed)"
                used = "misplaced"
            else:
                completeness = f"{r['completeness']:.1f}"
                rmsd = "n/a" if r["bb_rmsd"] is None else f"{r['bb_rmsd']:.2f}"
                best = num(c["completeness"]) > num(final_grades[(p, t)]["completeness"])
                used = "best copy" if best else "in place"
            tables["percopy"].append([p, TOOL_NAMES[t], f"{int(c['built_res']):,}",
                                      str(round(float(c["copies_est"]))),
                                      f"{float(c['in_density_visible_pct']):.1f}",
                                      completeness, rmsd, used])

    # alternative outputs (raw / unrefined), literal grading, mean over all targets
    tables["alternative15"] = []
    for tool, label, _, final in MODEL_FILES:
        if final:
            continue
        rows = [r for r in grades if r["tool"] == tool and r["model"] == label
                and r["reference"] == "deposited" and r["pdb"] in order]
        tables["alternative15"].append(
            [f"{TOOL_NAMES[tool]}, {label}"]
            + [mean([num(r[k]) for r in rows])
               for k in ("TM_ref", "completeness", "precision_all", "recall", "bb_rmsd")])

    tables["per_target_lddt"] = [
        [p] + [("" if num(ost[(p, t)]["lddt"]) is None else f"{num(ost[(p, t)]['lddt']):.3f}")
               for t in TOOLS]
        for p in order]

    others = target_sets["14"]
    out["facts"] = {
        "precision_all_vs_visible_cryozeta14": [
            mean([scores[(p, "cryozeta")]["precision_all"] for p in others]),
            mean([scores[(p, "cryozeta")]["precision_visible"] for p in others])],
        "mica_other14_cc_mask_min": min(num(phenix[(p, "mica")]["CC_mask"]) for p in others),
        "mica_other14_precision_visible_range": [
            min(scores[(p, "mica")]["precision_visible"] for p in others),
            max(scores[(p, "mica")]["precision_visible"] for p in others)],
        "mica_26DK_CC_mask": num(phenix[("26DK", "mica")]["CC_mask"]) if "26DK" in order else None,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(out, fh, indent=1)

    lines = ["# Report tables (generated by evaluation/make_report_tables.py)\n"]

    def add_table(title, header, rows, decimals):
        lines.append(f"\n## {title}\n")
        lines.append("| " + " | ".join(header) + " |")
        lines.append("|" + "---|" * len(header))
        for row in rows:
            lines.append("| " + " | ".join(
                x if isinstance(x, str) else ("" if x is None else f"{x:.{decimals}f}")
                for x in row) + " |")

    primary = ["Tool", "TM-score", "Completeness %", "Precision %", "Backbone RMSD (A)"]
    add_table("Headline, 14 maps", primary, tables["headline14"], 3)
    add_table("Headline, 15 maps", primary, tables["headline15"], 3)
    add_table("Secondary scores, 14 maps",
              ["Tool", "Recall %", "Amino-acid accuracy %", "C-alpha RMSD (A)", "F1 at 1 A %",
               "F1 at 3 A %", "Precision, all built %", "Aligned length %", "Breaks / 100 res"],
              tables["secondary14"], 2)
    add_table(f"OpenStructure, 14 maps (QS-score and DockQ over the {out['ost_multichain14']} "
              "multi-chain targets)", ["Tool", "lDDT", "Backbone lDDT", "QS-score", "DockQ"],
              tables["ost14"], 3)
    add_table("Phenix, 14 maps",
              ["Tool", "chain_comparison found %", "CC_mask", "Clashscore", "MolProbity score",
               "Rama outliers %", "Rama favored %", "Rotamer outliers %"],
              tables["phenix14"], 2)
    add_table("Literal grading, 15 maps (one-copy answer key, precision over all built residues)",
              ["Tool", "TM-score", "Completeness %", "Precision, all built %",
               "Backbone RMSD (A)"], tables["literal15"], 3)
    add_table("Alternative outputs, 15 maps, literal grading",
              ["Output", "TM-score", "Completeness %", "Precision, all built %", "Recall %",
               "Backbone RMSD (A)"], tables["alternative15"], 3)
    add_table("Per-copy targets",
              ["Target", "Tool", "Built residues", "Copies", "Visible residues in density %",
               "Completeness %", "Backbone RMSD (A)", "Copy used"], tables["percopy"], 1)
    names = [TOOL_NAMES[t] for t in TOOLS]
    for k, title in (("TM_ref", "TM-score"), ("completeness", "completeness %"),
                     ("precision_visible", "precision %"), ("bb_rmsd", "backbone RMSD (A)")):
        add_table(f"Per target: {title} (* per-copy value)", ["Target"] + names,
                  tables[f"per_target_{k}"], 1)
    add_table("Per target: lDDT", ["Target"] + names, tables["per_target_lddt"], 3)
    lines.append("\n## Facts\n```\n" + json.dumps(out["facts"], indent=1) + "\n```")
    args.out_md.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out_json} and {args.out_md}")


if __name__ == "__main__":
    main()
