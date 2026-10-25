#!/usr/bin/env python3
"""OpenStructure scores for every predicted model: lDDT, backbone lDDT, QS-score, DockQ, ICS, IPS.

Each model is compared with the protein chains of its answer key (deposited model; also the
symmetry assembly where there is one) by `ost compare-structures`. Models are first passed
through merge_fragments.py, because OpenStructure maps chains one to one and would otherwise
score only one fragment of a fragmented model. The answer key is reduced to its protein chains:
36PP also contains three DNA chains that no tool was asked to build.

The run can be interrupted: rows already in the output file are skipped on the next start. The
JSON written by OpenStructure for every pair is kept in the work directory.
"""
import argparse
import json
import shlex
import subprocess
import time
from pathlib import Path

import gemmi

from common import (MODEL_FILES, TsvWriter, add_target_arguments, is_amino_acid, load_targets,
                    model_path, read_tsv, reference_paths)
from merge_fragments import merge_fragments

SCORES = ["lddt", "bb_lddt", "qs_global", "qs_best", "dockq_ave", "dockq_wave", "ics", "ips"]
COLUMNS = (["pdb", "tool", "model", "reference", "status"] + SCORES
           + ["model_chains", "ref_chains", "mapped_chains", "seconds", "note"])


def protein_only(src, dst):
    """Copy of the answer key with protein chains only (ligands, waters, nucleic acids removed)."""
    if not dst.exists():
        st = gemmi.read_structure(str(src))
        st.setup_entities()
        st.remove_ligands_and_waters()
        st.remove_hydrogens()
        for chain in list(st[0]):
            if not any(is_amino_acid(r.name) for r in chain):
                st[0].remove_chain(chain.name)
        st.setup_entities()
        st.make_mmcif_document().write_file(str(dst))
    return dst


def compare_structures(ost, model, ref, out_json, timeout=1800):
    cmd = shlex.split(ost) + ["compare-structures", "-m", str(model), "-r", str(ref),
                              "-o", str(out_json), "--lddt", "--bb-lddt", "--qs-score",
                              "--dockq", "--ics", "--ips", "-ft"]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT"}
    try:
        with open(out_json) as fh:
            return json.load(fh)
    except (OSError, ValueError) as err:
        return {"status": "NO_OUTPUT", "exception": str(err)[:100]}


def score(model_file, ref_file, stem, work, ost):
    """One output row (without the identifying columns) for one model and reference."""
    row = {}
    try:
        merged = work / "merged" / f"{stem}.cif"
        merge_fragments(model_file, ref_file, merged)
        st = gemmi.read_structure(str(merged))
        n_res, n_chains = sum(len(chain) for chain in st[0]), len(st[0])
    except Exception as err:   # report the failure in the table and carry on with the next model
        n_res, n_chains = 0, 0
        row["note"] = f"merge failed: {str(err)[:60]}"
    ref = protein_only(ref_file, work / "references" / ref_file.name)
    result = (compare_structures(ost, merged, ref, work / "json" / f"{stem}.json")
              if n_res else {"status": "EMPTY"})
    row["status"] = result.get("status", "?")
    for key in SCORES:
        value = result.get(key)
        row[key] = f"{value:.3f}" if isinstance(value, (int, float)) else ""
    row["model_chains"] = n_chains
    row["ref_chains"] = len(result.get("reference_chains") or [])
    row["mapped_chains"] = len(result.get("chain_mapping") or {})
    if result.get("status") != "SUCCESS":
        row["note"] = (row.get("note", "") + " " + str(result.get("exception", ""))[:80]).strip()
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="data directory (see README)")
    parser.add_argument("--ost", default="conda run -n ost ost",
                        help="command that starts OpenStructure (default: %(default)s)")
    parser.add_argument("--work", type=Path, default=Path("results/work/openstructure"),
                        help="directory for merged models, references and JSON output")
    parser.add_argument("--out", type=Path, default=Path("results/openstructure_scores.tsv"))
    add_target_arguments(parser)
    args = parser.parse_args()

    for sub in ("merged", "references", "json"):
        (args.work / sub).mkdir(parents=True, exist_ok=True)
    done = set()
    if args.out.exists():
        done = {(r["pdb"], r["tool"], r["model"], r["reference"]) for r in read_tsv(args.out)}
    out = TsvWriter(args.out, COLUMNS, append=True)
    targets = load_targets(args.targets_file, args.targets)
    # final models first, alternative outputs afterwards
    ordered = sorted(MODEL_FILES, key=lambda m: not m[3])
    for tool, label, file_name, _ in ordered:
        for target in targets:
            pdb = target["pdb"]
            model_file = model_path(args.data, pdb, tool, file_name)
            for ref_name, ref_file in reference_paths(args.data, pdb):
                if (pdb, tool, label, ref_name) in done:
                    continue
                start = time.time()
                if not model_file.exists() or model_file.stat().st_size == 0:
                    row = dict(status="MISSING_MODEL", lddt=0, bb_lddt=0, qs_global=0,
                               dockq_ave=0)
                else:
                    row = score(model_file, ref_file, f"{pdb}_{tool}_{label}_{ref_name}",
                                args.work, args.ost)
                    row["seconds"] = int(time.time() - start)
                row.update(pdb=pdb, tool=tool, model=label, reference=ref_name)
                out.write(row)
                print(f"{pdb} {tool:12} {label:9} vs {ref_name:9} {row['status']:8} "
                      f"lDDT {row['lddt']:>6} QS {row['qs_global']:>6} "
                      f"DockQ {row['dockq_ave']:>6} {row.get('seconds', '')}s", flush=True)
    out.close()
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
