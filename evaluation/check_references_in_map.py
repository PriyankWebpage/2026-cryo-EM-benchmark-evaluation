#!/usr/bin/env python3
"""Check that every answer key sits in its map's density.

For each C-alpha atom of the answer key, the map value at its position (trilinear interpolation)
is compared with the recommended EMDB contour level. Reported per target: the percentage of
C-alpha atoms at or above the contour, the median map value at the C-alpha atoms relative to the
contour, and the chance level (percentage of all grid points above the contour). A correct answer
key lies far above the chance level.
"""
import argparse
from pathlib import Path

import gemmi
import numpy as np

from common import (TsvWriter, add_target_arguments, is_amino_acid, load_targets, map_path,
                    read_map, reference_paths)

COLUMNS = ["pdb", "reference", "contour_level", "chance_pct", "ca_atoms", "in_density_pct",
           "median_over_contour"]


def ca_positions(path):
    st = gemmi.read_structure(str(path))
    positions = []
    for chain in st[0]:
        for res in chain:
            atom = res.find_atom("CA", "*")
            if atom is not None and is_amino_acid(res.name):
                positions.append(atom.pos)
    return positions


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="data directory (see README)")
    parser.add_argument("--out", type=Path, default=Path("results/reference_check.tsv"))
    add_target_arguments(parser)
    args = parser.parse_args()

    out = TsvWriter(args.out, COLUMNS)
    for target in load_targets(args.targets_file, args.targets):
        grid = read_map(map_path(args.data, target["emd"]))
        level = float(target["contour_level"])
        chance = 100.0 * np.mean(np.asarray(grid) >= level)
        for ref_name, ref_file in reference_paths(args.data, target["pdb"]):
            positions = ca_positions(ref_file)
            values = np.array([grid.interpolate_value(p) for p in positions])
            row = dict(pdb=target["pdb"], reference=ref_name, contour_level=level,
                       chance_pct=f"{chance:.1f}", ca_atoms=len(positions),
                       in_density_pct=f"{100 * np.mean(values >= level):.1f}",
                       median_over_contour=f"{np.median(values) / level:.2f}")
            out.write(row)
            print(f"{row['pdb']:5} {ref_name:9} C-alpha in density {row['in_density_pct']:>5} % "
                  f"(chance level {row['chance_pct']:>4} %)",
                  flush=True)
    out.close()
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
