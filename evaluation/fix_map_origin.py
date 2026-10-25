#!/usr/bin/env python3
"""Move models into the map frame where the tool ignored the map's start offset.

A map grid does not always start at index (0, 0, 0): the header field NSTART gives the index of
the first grid point on each axis. A program that ignores it writes coordinates shifted by
NSTART x voxel size. Of the 15 maps only EMD-80558 (26DK) has a non-zero start,
(-80, -80, -80) voxels of 1.66 A, and ModelAngelo ignores it, so both of its 26DK models are
off by 132.8 A on each axis. This script translates every model listed in common.FRAME_FIXED by
the offset read from the map header and writes it next to the original with the suffix
"_frame_fixed". Atoms are only translated; the original files are left untouched.
"""
import argparse
from pathlib import Path

import gemmi

from common import (DEFAULT_TARGETS, FRAME_FIXED, MODEL_FILES, load_targets, map_path,
                    model_path)


def start_offset(map_file):
    """Offset (x, y, z) in A of the first grid point, from the header fields NSTART and MAPC/R/S."""
    ccp4 = gemmi.read_ccp4_map(str(map_file))
    nstart = [ccp4.header_i32(i) for i in (5, 6, 7)]              # column, row, section
    axis = [ccp4.header_i32(i) - 1 for i in (17, 18, 19)]         # axis of column, row, section
    sampling = [ccp4.header_i32(i) for i in (8, 9, 10)]           # grid intervals along x, y, z
    cell = [ccp4.header_float(i) for i in (11, 12, 13)]           # cell lengths along x, y, z
    offset = [0.0, 0.0, 0.0]
    for start, ax in zip(nstart, axis):
        offset[ax] = start * cell[ax] / sampling[ax]
    # the header stores the cell in single precision; round to the 0.001 A of coordinate files
    return [round(x, 3) for x in offset]


def translate(src, dst, offset):
    st = gemmi.read_structure(str(src))
    shift = gemmi.Transform(gemmi.Mat33(), gemmi.Vec3(*offset))
    for model in st:
        model.transform_pos_and_adp(shift)
    if dst.suffix == ".pdb":
        st.write_pdb(str(dst))
    else:
        st.make_mmcif_document().write_file(str(dst))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="data directory (see README)")
    parser.add_argument("--targets-file", type=Path, default=DEFAULT_TARGETS,
                        help="target table (default: config/targets.tsv)")
    args = parser.parse_args()

    targets = {t["pdb"]: t for t in load_targets(args.targets_file)}
    for pdb, tool in sorted(FRAME_FIXED):
        offset = start_offset(map_path(args.data, targets[pdb]["emd"]))
        for model_tool, _, file_name, _ in MODEL_FILES:
            if model_tool != tool:
                continue
            src = args.data / "predicted_models" / tool / pdb / file_name
            dst = model_path(args.data, pdb, tool, file_name)
            if not src.exists():
                print(f"{pdb} {tool}: {src.name} not found, skipped")
                continue
            translate(src, dst, offset)
            print(f"{pdb} {tool}: {src.name} -> {dst.name}, shifted by "
                  f"({offset[0]:.2f}, {offset[1]:.2f}, {offset[2]:.2f}) A")


if __name__ == "__main__":
    main()
