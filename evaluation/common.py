"""Definitions and structure helpers shared by the evaluation scripts.

All scripts read one data directory (DATA) laid out like the Zenodo record:

    DATA/maps/emd_<EMDB>.map[.gz]                      cryo-EM maps given to the tools
    DATA/reference_models/<PDB>_deposited_protein.cif  answer key (deposited model, no ligands
                                                       or waters)
    DATA/reference_models/<PDB>_assembly_protein.cif   biological assembly 1, where it has
                                                       several copies
    DATA/predicted_models/<tool>/<PDB>/<file>          models built by the tools (MODEL_FILES)
"""
import csv
import re
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

REPO_DIR = Path(__file__).resolve().parent.parent
DEFAULT_TARGETS = REPO_DIR / "config" / "targets.tsv"

TOOLS = ["modelangelo", "cryozeta", "cryoatom", "emodelx", "mica"]
TOOL_NAMES = {"modelangelo": "ModelAngelo", "cryozeta": "CryoZeta", "cryoatom": "CryoAtom v1",
              "emodelx": "EModelX(+AF)", "mica": "MICA"}

# Output files of each tool: (tool, label, file name, is the tool's final answer).
# "raw" (ModelAngelo, CryoAtom: before pruning) and "unrefined" (EModelX, MICA: before Phenix
# refinement) are alternative outputs, graded for information only.
MODEL_FILES = [
    ("modelangelo", "final", "final_pruned.cif", True),
    ("modelangelo", "raw", "raw_unpruned.cif", False),
    ("cryozeta", "final", "final.cif", True),
    ("cryoatom", "final", "final_pruned.cif", True),
    ("cryoatom", "raw", "raw_unpruned.cif", False),
    ("emodelx", "final", "final_refined.pdb", True),
    ("emodelx", "unrefined", "unrefined.pdb", False),
    ("mica", "final", "final_refined.pdb", True),
    ("mica", "unrefined", "unrefined.pdb", False),
]
FINAL_MODEL_FILES = [(tool, name) for tool, _, name, final in MODEL_FILES if final]

# Models written without the map's start offset (header field NSTART). They are graded after
# fix_map_origin.py has moved them into the map frame (file name suffix "_frame_fixed").
FRAME_FIXED = {("26DK", "modelangelo")}

# Maps that hold more copies of the protein than the deposited answer key (26DK: about 47,
# 9UBC: 2, 9I8C: 12). These targets are also graded per copy (grade_symmetry_copies.py).
PER_COPY_TARGETS = ["26DK", "9UBC", "9I8C"]

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
PAIR_CUTOFF = 3.0   # A; a built C-alpha within this distance of an answer-key C-alpha is a match
MIN_IDENTITY = 30.0  # %; below this a model chain is not considered the same protein


def load_targets(path=DEFAULT_TARGETS, only=None):
    """Targets from the TSV file, in file order; `only` restricts them to the given PDB ids."""
    with open(path, newline="") as fh:
        targets = list(csv.DictReader(fh, delimiter="\t"))
    if only:
        unknown = set(only) - {t["pdb"] for t in targets}
        if unknown:
            raise SystemExit(f"unknown target(s): {', '.join(sorted(unknown))}")
        targets = [t for t in targets if t["pdb"] in only]
    return targets


def add_target_arguments(parser):
    """Command-line options shared by all scripts that loop over targets."""
    parser.add_argument("--targets-file", type=Path, default=DEFAULT_TARGETS,
                        help="target table (default: config/targets.tsv)")
    parser.add_argument("--targets", nargs="+", metavar="PDB",
                        help="grade only these targets (default: all)")


def model_path(data_dir, pdb, tool, file_name):
    """Path of a predicted model; the frame-corrected file where the tool ignored the map offset."""
    if (pdb, tool) in FRAME_FIXED:
        name = Path(file_name)
        file_name = f"{name.stem}_frame_fixed{name.suffix}"
    return Path(data_dir) / "predicted_models" / tool / pdb / file_name


def reference_paths(data_dir, pdb):
    """[(name, path)]: the deposited answer key, plus the symmetry assembly where there is one."""
    folder = Path(data_dir) / "reference_models"
    refs = [("deposited", folder / f"{pdb}_deposited_protein.cif")]
    assembly = folder / f"{pdb}_assembly_protein.cif"
    if assembly.exists():
        refs.append(("assembly", assembly))
    return refs


def map_path(data_dir, emd):
    """The target's map, plain or gzip-compressed."""
    for name in (f"emd_{emd}.map", f"emd_{emd}.map.gz"):
        path = Path(data_dir) / "maps" / name
        if path.exists():
            return path
    raise FileNotFoundError(f"no map for EMD-{emd} in {Path(data_dir) / 'maps'}")


def read_map(path):
    """Map as a gemmi grid in the model frame (axis order, origin and start offset applied)."""
    ccp4 = gemmi.read_ccp4_map(str(path))
    ccp4.setup(float("nan"))
    return ccp4.grid


def in_density(grid, level, points):
    """Percentage of points whose interpolated map value is at or above the contour level."""
    values = np.array([grid.interpolate_value(gemmi.Position(*p)) for p in points])
    return 100.0 * np.mean(values >= level)


def is_amino_acid(name):
    info = gemmi.find_tabulated_residue(name)
    return info is not None and info.is_amino_acid()


def residues(path):
    """Amino-acid residues that have a C-alpha, in file order.

    Each residue is a dict: chain, seqnum, name1 (one-letter code), ca (coordinates) and
    bb (coordinates of the backbone atoms N, CA, C, O that are present).
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    out = []
    for chain in st[0]:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            ca = res.find_atom("CA", "*")
            if ca is None:
                continue
            bb = {}
            for name in ("N", "CA", "C", "O"):
                atom = res.find_atom(name, "*")
                if atom is not None:
                    bb[name] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
            out.append(dict(chain=chain.name, seqnum=res.seqid.num,
                            name1=info.one_letter_code.upper(),
                            ca=np.array([ca.pos.x, ca.pos.y, ca.pos.z]), bb=bb))
    return out


def chain_indices(res):
    """Chain name -> indices of its residues in `res`."""
    chains = {}
    for i, r in enumerate(res):
        chains.setdefault(r["chain"], []).append(i)
    return chains


def sequence(res):
    return "".join(r["name1"] if r["name1"] in AMINO_ACIDS else "X" for r in res)


def align(query, target):
    """Global sequence alignment (gemmi, BLOSUM62).

    Returns the aligned position pairs (i in query, j in target) and the identity in percent
    of the query length.
    """
    al = gemmi.align_string_sequences(list(query), list(target), [], gemmi.AlignmentScoring("b"))
    pairs, i, j = [], 0, 0
    for count, op in re.findall(r"(\d+)([MID])", al.cigar_str()):
        n = int(count)
        if op == "M":
            pairs += [(i + k, j + k) for k in range(n)]
            i += n
            j += n
        elif op == "I":   # query only
            i += n
        else:             # target only
            j += n
    return pairs, al.calculate_identity(1)


def positions_in_reference(model_res, ref_res):
    """Flags over the model residues: True where the residue's sequence position is modelled in
    the reference ("visible" residues). Each model chain is aligned to the reference chain with
    the highest sequence identity; chains below MIN_IDENTITY match nothing.
    """
    ref_seqs = [sequence([ref_res[i] for i in idx]) for idx in chain_indices(ref_res).values()]
    flags = np.zeros(len(model_res), dtype=bool)
    for idx in chain_indices(model_res).values():
        seq = sequence([model_res[i] for i in idx])
        best = None
        for ref_seq in ref_seqs:
            pairs, identity = align(seq, ref_seq)
            if best is None or identity > best[1]:
                best = (pairs, identity)
        if best is not None and best[1] >= MIN_IDENTITY:
            for i, _ in best[0]:
                flags[idx[i]] = True
    return flags


def greedy_pairs(ref_xyz, model_xyz, cutoff):
    """One-to-one pairing of reference and model points by increasing distance, up to `cutoff`.

    Returns a list of (reference index, model index, distance).
    """
    if len(ref_xyz) == 0 or len(model_xyz) == 0:
        return []
    k = min(len(model_xyz), 8)
    dist, idx = cKDTree(model_xyz).query(ref_xyz, k=k, distance_upper_bound=cutoff)
    dist = dist.reshape(len(ref_xyz), -1)
    idx = idx.reshape(len(ref_xyz), -1)
    candidates = sorted((dist[i, c], i, idx[i, c])
                        for i in range(len(ref_xyz)) for c in range(k) if np.isfinite(dist[i, c]))
    used_ref, used_model, pairs = set(), set(), []
    for d, i, j in candidates:
        if i in used_ref or j in used_model:
            continue
        used_ref.add(i)
        used_model.add(j)
        pairs.append((i, j, d))
    return pairs


def chain_breaks(res):
    """Consecutive C-alpha atoms more than 4 A apart within a chain, in residue-number order."""
    breaks = 0
    for idx in chain_indices(res).values():
        chain = sorted((res[i] for i in idx), key=lambda r: r["seqnum"])
        breaks += sum(np.linalg.norm(a["ca"] - b["ca"]) > 4.0 for a, b in zip(chain, chain[1:]))
    return int(breaks)


def backbone_rmsd(pairs, ref_res, model_res, transform=None):
    """RMSD over the backbone atoms (N, CA, C, O) present in both residues of each pair.

    `transform` = (rotation, translation) is applied to the model atoms first.
    """
    sq = []
    for i, j, _ in pairs:
        for name in ("N", "CA", "C", "O"):
            if name in ref_res[i]["bb"] and name in model_res[j]["bb"]:
                x = model_res[j]["bb"][name]
                if transform is not None:
                    x = transform[0] @ x + transform[1]
                sq.append(np.sum((ref_res[i]["bb"][name] - x) ** 2))
    return float(np.sqrt(np.mean(sq))) if sq else None


def f1(ref_xyz, model_xyz, cutoff):
    """F1 (percent) of recall and precision at one distance cutoff."""
    n = len(greedy_pairs(ref_xyz, model_xyz, cutoff))
    recall, precision = n / len(ref_xyz), n / max(len(model_xyz), 1)
    return 100 * 2 * precision * recall / (precision + recall) if precision + recall else 0.0


class TsvWriter:
    """Writes rows (dicts) to a TSV file as they are produced, so partial results survive."""

    def __init__(self, path, columns, append=False):
        self.columns = columns
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        new = not (append and path.exists())
        self.fh = open(path, "w" if new else "a")
        if new:
            self.write(dict(zip(columns, columns)))

    def write(self, row):
        self.fh.write("\t".join(str(row.get(c, "")) for c in self.columns) + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
