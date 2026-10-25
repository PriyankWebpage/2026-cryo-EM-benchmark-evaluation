#!/usr/bin/env python3
"""Phenix scores for the final model of every tool.

  phenix.chain_comparison   C-alpha match in place against the deposited answer key
                            (% of answer-key residues found, sequence match, fragments, ...)
  phenix.map_model_cc       model-to-map correlation (CC_mask, CC_box, CC_peaks, CC_volume) at
                            the resolution of the deposited map
  phenix.molprobity         geometry (clashscore, Ramachandran, rotamers, MolProbity score)

Phenix 1.20 does not read the tools' minimal mmCIF files, so models and answer keys are first
converted to PDB format with gemmi (atoms unchanged). The run can be interrupted: rows already in
the output file are skipped on the next start. The Phenix logs are kept in the work directory.
"""
import argparse
import gzip
import itertools
import re
import shlex
import shutil
import string
import subprocess
import time
from pathlib import Path

import gemmi

from common import (FINAL_MODEL_FILES, TsvWriter, add_target_arguments, load_targets, map_path,
                    model_path, read_tsv, reference_paths)

CHAIN_COMPARISON = ["cc_close_rmsd", "cc_close_n", "cc_far_n", "cc_forward", "cc_reverse",
                    "cc_mixed", "cc_found_pct", "cc_ca_score", "cc_seq_match_pct",
                    "cc_seq_score", "cc_mean_length", "cc_fragments", "cc_bad_connections"]
COLUMNS = (["pdb", "tool", "status", "ref_res", "model_res", "cc_found_pct", "cc_close_rmsd",
            "cc_close_n", "cc_far_n", "cc_forward", "cc_reverse", "cc_mixed", "cc_ca_score",
            "cc_seq_match_pct", "cc_seq_score", "cc_mean_length", "cc_fragments",
            "cc_bad_connections", "CC_mask", "CC_box", "CC_peaks", "CC_volume",
            "clashscore", "rama_outliers_pct", "rama_favored_pct", "rota_outliers_pct",
            "cbeta_dev", "rms_bonds", "rms_angles", "molprobity_score", "seconds", "note"])


def run_phenix(env, command, log, cwd, timeout=3600):
    """Run a Phenix command inside the Phenix environment; return its output (also saved to log)."""
    script = f"source {shlex.quote(str(env))} >/dev/null 2>&1; {command}"
    try:
        result = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                                timeout=timeout, cwd=cwd)
        output = result.stdout + "\n" + result.stderr
    except subprocess.TimeoutExpired:
        output = "TIMEOUT"
    Path(log).write_text(output)
    return output


def first_number(pattern, text):
    match = re.search(pattern, text)
    return match.group(1) if match else ""


def parse_chain_comparison(text):
    """Summary row of phenix.chain_comparison (target = answer key, query = model): label followed
    by RMSD, N close, N far, forward, reverse, mixed, % found, CA score, sequence match %,
    sequence score, mean length, fragments, bad connections."""
    scores = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) >= 14 and all(re.fullmatch(r"-?\d+(\.\d+)?", x) for x in fields[1:14]):
            scores = dict(zip(CHAIN_COMPARISON, fields[1:14]))
    scores["ref_res"] = first_number(r"Original residues in target:\s*(\d+)", text)
    scores["model_res"] = first_number(r"Original residues in target:\s*\d+\s+In query:\s*(\d+)",
                                       text)
    return scores


def parse_map_model_cc(text):
    return {key: first_number(rf"{key}\s*[:=]\s*([-\d.]+)", text)
            for key in ("CC_mask", "CC_box", "CC_peaks", "CC_volume")}


def parse_molprobity(text):
    return {"clashscore": first_number(r"Clashscore\s*=\s*([\d.]+)", text),
            "rama_outliers_pct": first_number(r"Ramachandran outliers\s*=\s*([\d.]+)\s*%", text),
            "rama_favored_pct": first_number(r"(?m)^\s*favored\s*=\s*([\d.]+)\s*%", text),
            "rota_outliers_pct": first_number(r"Rotamer outliers\s*=\s*([\d.]+)\s*%", text),
            "cbeta_dev": first_number(r"C-beta deviations\s*=\s*(\d+)", text),
            "rms_bonds": first_number(r"RMS\(bonds\)\s*=\s*([\d.]+)", text),
            "rms_angles": first_number(r"RMS\(angles\)\s*=\s*([\d.]+)", text),
            "molprobity_score": first_number(r"MolProbity score\s*=\s*([\d.]+)", text)}


def to_pdb(src, dst):
    """PDB-format copy of a structure. PDB chain ids hold at most two characters, so longer names
    (ModelAngelo fragments such as "Aaa") are replaced by unused one- or two-character ids; chain
    names play no role in the three Phenix scores."""
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    names = [chain.name for chain in st[0]]
    if any(len(n) > 2 for n in names):
        chars = string.ascii_uppercase + string.ascii_lowercase + string.digits
        pool = list(chars) + ["".join(p) for p in itertools.product(chars, repeat=2)]
        used = {n for n in names if len(n) <= 2}
        free = (c for c in pool if c not in used)
        for chain in st[0]:
            if len(chain.name) > 2:
                chain.name = next(free)
    st.write_pdb(str(dst))
    return dst


def plain_map(path, work):
    """Path of an uncompressed copy of the map (Phenix does not read gzip files)."""
    if path.suffix != ".gz":
        return path
    out = work / path.stem
    if not out.exists():
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return out


def score(model_file, ref_pdb, map_file, resolution, folder, env):
    """One output row (without pdb, tool and seconds) for one final model."""
    folder.mkdir(parents=True, exist_ok=True)
    model = to_pdb(model_file, folder / "model.pdb")
    q = shlex.quote
    row, notes = {}, []
    text = run_phenix(env, f"phenix.chain_comparison {q(str(ref_pdb))} {q(str(model))} "
                           "unique_part_of_target_only=False",
                      folder / "chain_comparison.log", folder)
    row.update(parse_chain_comparison(text))
    if not row.get("cc_found_pct"):
        notes.append("chain_comparison unparsed")
    text = run_phenix(env, f"phenix.map_model_cc {q(str(model))} {q(str(map_file))} "
                           f"resolution={resolution}",
                      folder / "map_model_cc.log", folder)
    row.update(parse_map_model_cc(text))
    if not row["CC_mask"]:
        notes.append("map_model_cc unparsed")
    text = run_phenix(env, f"phenix.molprobity {q(str(model))} output.prefix=mp",
                      folder / "molprobity.log", folder)
    row.update(parse_molprobity(text))
    if not row["molprobity_score"]:
        notes.append("molprobity unparsed")
    row["status"] = "PARTIAL" if notes else "OK"
    row["note"] = "; ".join(notes)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="data directory (see README)")
    parser.add_argument("--phenix-env", type=Path, required=True,
                        help="phenix_env.sh of the Phenix installation")
    parser.add_argument("--work", type=Path, default=Path("results/work/phenix"),
                        help="directory for converted files, maps and Phenix logs")
    parser.add_argument("--out", type=Path, default=Path("results/phenix_scores.tsv"))
    add_target_arguments(parser)
    args = parser.parse_args()

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    done = {(r["pdb"], r["tool"]) for r in read_tsv(args.out)} if args.out.exists() else set()
    out = TsvWriter(args.out, COLUMNS, append=True)
    for target in load_targets(args.targets_file, args.targets):
        pdb = target["pdb"]
        if all((pdb, tool) in done for tool, _ in FINAL_MODEL_FILES):
            continue
        map_file = plain_map(map_path(args.data, target["emd"]).resolve(), work)
        ref_pdb = to_pdb(reference_paths(args.data, pdb)[0][1], work / f"{pdb}_reference.pdb")
        for tool, file_name in FINAL_MODEL_FILES:
            if (pdb, tool) in done:
                continue
            start = time.time()
            model_file = model_path(args.data, pdb, tool, file_name)
            if not model_file.exists() or model_file.stat().st_size == 0:
                row = {"status": "MISSING_MODEL"}
            else:
                row = score(model_file, ref_pdb, map_file, target["resolution"],
                            work / f"{pdb}_{tool}", args.phenix_env)
            row.update(pdb=pdb, tool=tool, seconds=int(time.time() - start))
            out.write(row)
            print(f"{pdb} {tool:12} {row['status']:8} CC_mask {row.get('CC_mask', ''):>6} "
                  f"clashscore {row.get('clashscore', ''):>6} "
                  f"MolProbity {row.get('molprobity_score', ''):>5} "
                  f"found {row.get('cc_found_pct', ''):>6} {row['seconds']}s", flush=True)
    out.close()
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
