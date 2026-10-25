#!/usr/bin/env python3
"""Prepare a predicted model for OpenStructure by joining its fragments, without moving any atom.

OpenStructure maps model chains to reference chains one to one, so a correct model built in
several fragments (ModelAngelo, CryoAtom) would have only one fragment scored. Here every
fragment is assigned to the answer-key chain it sits on, and fragments of the same answer-key
chain are joined into one chain with that chain's name.

How a fragment is assigned (in the map frame, nothing is moved):
  1. In place: for each residue, the nearest answer-key C-alpha within 3 A gives the chain and
     residue number. The fragment goes to the chain that most of its residues sit on and is
     renumbered with one constant offset so that its numbering agrees with the answer key.
  2. Otherwise (an extra symmetry copy, or a misplaced fragment): the answer-key chain with the
     most similar sequence, keeping the tool's residue numbers.
Only in-place fragments are joined, and only with each other. If an in-place fragment repeats
residue numbers already present in its chain for fewer than half of its residues (a few
residues built twice), those residues are dropped from the later, smaller fragment; a larger
overlap means a second copy in place, which becomes a separate chain <name><n>.
Rule-2 fragments of an answer-key chain that already has in-place fragments are dropped: they are
not part of the answer key, and keeping them would let OpenStructure's chain mapping pick a stray
copy instead of the in-place chain. Rule-2 fragments are kept, one chain each, only when nothing
was built in place for that answer-key chain (a model built on another symmetry copy). Residues
that match no answer-key chain, ligands, waters and hydrogens are dropped.
"""
import argparse
from collections import Counter

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from common import AMINO_ACIDS, MIN_IDENTITY, PAIR_CUTOFF, is_amino_acid


def protein_chains(st):
    """[(chain name, amino-acid residues)] of the first model, skipping chains without any."""
    chains = []
    for chain in st[0]:
        res = [r for r in chain if is_amino_acid(r.name)]
        if res:
            chains.append((chain.name, res))
    return chains


def one_letter_sequence(res):
    codes = (gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in res)
    return "".join(c if c in AMINO_ACIDS else "X" for c in codes)


def ca_position(res):
    atom = res.find_atom("CA", "*")
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z]) if atom else None


def identity(query, target):
    al = gemmi.align_string_sequences(list(query), list(target), [], gemmi.AlignmentScoring("b"))
    return al.calculate_identity(1)


def renumbered_copy(res, number):
    new = gemmi.Residue()
    new.name = res.name
    new.seqid = gemmi.SeqId(number, " ")
    for atom in res:
        new.add_atom(atom)
    return new


def assign_fragments(model_chains, ref_chains):
    """Fragments as (in_place, answer-key chain, residue numbers, residues); unmatched residues."""
    points, labels, ref_seqs = [], [], {}
    for name, res in ref_chains:
        ref_seqs[name] = one_letter_sequence(res)
        for r in res:
            xyz = ca_position(r)
            if xyz is not None:
                points.append(xyz)
                labels.append((name, r.seqid.num))
    tree = cKDTree(np.array(points))

    fragments, unmatched = [], []
    for _, res in model_chains:
        nearest = []   # per residue: (chain, number) of the answer-key C-alpha within 3 A, or None
        for r in res:
            xyz = ca_position(r)
            if xyz is None:
                nearest.append(None)
                continue
            dist, i = tree.query(xyz)
            nearest.append(labels[i] if dist <= PAIR_CUTOFF else None)
        hits = Counter(n[0] for n in nearest if n is not None)
        if hits and hits.most_common(1)[0][1] >= max(3, 0.5 * len(res)):
            chain = hits.most_common(1)[0][0]
            # the tools number residues by sequence position: keep that numbering, shifted by the
            # median offset to the answer-key numbering over the residues that sit on the chain
            offsets = [n[1] - r.seqid.num for n, r in zip(nearest, res)
                       if n is not None and n[0] == chain]
            offset = int(round(float(np.median(offsets)))) if offsets else 0
            fragments.append((True, chain, [r.seqid.num + offset for r in res], res))
        else:
            seq = one_letter_sequence(res)
            scores = {name: identity(seq, ref_seq) for name, ref_seq in ref_seqs.items()}
            chain = max(scores, key=scores.get)
            if scores[chain] < MIN_IDENTITY:
                unmatched += res
                continue
            fragments.append((False, chain, [r.seqid.num for r in res], res))
    return fragments, unmatched


def join_fragments(fragments):
    """Chain name -> {residue number: residue}; also returns the number of trimmed residues."""
    in_place_chains = {chain for in_place, chain, _, _ in fragments if in_place}
    chains, trimmed = {}, 0
    # in-place fragments first, largest first, so the largest keeps residues built twice
    for in_place, target, numbers, res in sorted(fragments, key=lambda f: (not f[0], -len(f[3]))):
        if not in_place and target in in_place_chains:
            continue
        name, k = target, 2
        if in_place:
            repeated = [n in chains.get(name, {}) for n in numbers]
            if 0 < sum(repeated) < 0.5 * len(numbers):
                res = [r for r, rep in zip(res, repeated) if not rep]
                numbers = [n for n, rep in zip(numbers, repeated) if not rep]
                trimmed += sum(repeated)
        while (not in_place and name in chains) or any(n in chains.get(name, {}) for n in numbers):
            name = f"{target}{k}"
            k += 1
        chains.setdefault(name, {})
        for r, number in zip(res, numbers):
            chains[name][number] = renumbered_copy(r, number)
    return chains, trimmed


def merge_fragments(model_file, ref_file, out_file):
    """Write the merged model to `out_file` (mmCIF) and return a one-line summary."""
    ref = gemmi.read_structure(str(ref_file))
    ref.setup_entities()
    mod = gemmi.read_structure(str(model_file))
    mod.setup_entities()
    mod.remove_ligands_and_waters()
    mod.remove_hydrogens()

    ref_chains = protein_chains(ref)
    fragments, unmatched = assign_fragments(protein_chains(mod), ref_chains)
    in_place_chains = {chain for in_place, chain, _, _ in fragments if in_place}
    dropped = len(unmatched) + sum(len(res) for in_place, chain, _, res in fragments
                                   if not in_place and chain in in_place_chains)
    chains, trimmed = join_fragments(fragments)

    new = gemmi.Structure()
    new.cell = mod.cell
    new.spacegroup_hm = mod.spacegroup_hm
    model = gemmi.Model("1")
    ref_names = [name for name, _ in ref_chains]
    for name in ref_names + sorted(n for n in chains if n not in ref_names):
        if name in chains:
            chain = gemmi.Chain(name)
            for number in sorted(chains[name]):
                chain.add_residue(chains[name][number])
            model.add_chain(chain)
    new.add_model(model)
    new.setup_entities()
    new.assign_label_seq_id()
    new.make_mmcif_document().write_file(str(out_file))

    n_fragments = len(fragments) + (1 if unmatched else 0)
    n_in_place = sum(1 for in_place, *_ in fragments if in_place)
    return (f"{n_fragments} fragments ({n_in_place} placed in the map frame) -> "
            f"{len(model)} chains ({', '.join(ch.name for ch in model)}); "
            f"{sum(len(ch) for ch in model)} residues kept, "
            f"{dropped} residues dropped (extra copies / unmatched), "
            f"{trimmed} overlapping residues trimmed")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("model", help="predicted model (mmCIF or PDB)")
    parser.add_argument("reference", help="answer key (mmCIF)")
    parser.add_argument("out", help="merged model (mmCIF)")
    args = parser.parse_args()
    print(merge_fragments(args.model, args.reference, args.out))


if __name__ == "__main__":
    main()
