#!/usr/bin/env python3
"""Build the answer keys from the PDB.

For each target the entry is downloaded from RCSB and written as
    <PDB>_deposited_protein.cif   the deposited model with ligands and waters removed
    <PDB>_assembly_protein.cif    biological assembly 1 expanded by its symmetry operators,
                                  only when it has more than one copy (9I8C)
The Zenodo record already contains the answer keys used for the published grades. Entries can be
revised after release, so files rebuilt later may differ from those.
"""
import argparse
import urllib.request
from pathlib import Path

import gemmi

from common import add_target_arguments, is_amino_acid, load_targets

URL = "https://files.rcsb.org/download/{}.cif"


def protein_residues(st):
    return sum(1 for chain in st[0] for res in chain if is_amino_acid(res.name))


def make_references(pdb, out_dir, cache_dir):
    entry = cache_dir / f"{pdb}.cif"
    if not entry.exists() or entry.stat().st_size == 0:
        urllib.request.urlretrieve(URL.format(pdb), entry)
    st = gemmi.read_structure(str(entry))
    st.setup_entities()
    st.remove_ligands_and_waters()
    n_deposited = protein_residues(st)
    st.make_mmcif_document().write_file(str(out_dir / f"{pdb}_deposited_protein.cif"))

    assembly = next((a for a in st.assemblies if a.name == "1"), None)
    copies = max(1, sum(len(g.operators) for g in assembly.generators)) if assembly else 1
    if copies == 1:
        return f"{pdb}\tdeposited {n_deposited} protein residues\tassembly 1 = deposited"
    st.transform_to_assembly("1", gemmi.HowToNameCopiedChain.AddNumber)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(out_dir / f"{pdb}_assembly_protein.cif"))
    return (f"{pdb}\tdeposited {n_deposited} protein residues\t"
            f"assembly 1: {protein_residues(st)} ({copies} copies)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True,
                        help="data directory; answer keys go to DATA/reference_models")
    parser.add_argument("--cache", type=Path, default=None,
                        help="where downloaded entries are kept (default: DATA/pdb_entries)")
    add_target_arguments(parser)
    args = parser.parse_args()

    out_dir = args.data / "reference_models"
    cache_dir = args.cache or args.data / "pdb_entries"
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    for target in load_targets(args.targets_file, args.targets):
        print(make_references(target["pdb"], out_dir, cache_dir), flush=True)


if __name__ == "__main__":
    main()
