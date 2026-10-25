# Evaluation of five cryo-EM model-building tools on 15 structures

This repository contains the code used to grade the atomic models built by five de novo
cryo-EM model-building tools for 15 cryo-EM maps released in 2026. Each model is compared with
the structure deposited in the PDB for the same map (the answer key).

The input data, all predicted models and the published grade tables are on Zenodo:
[10.5281/zenodo.23199153](https://doi.org/10.5281/zenodo.23199153).

| Tool | Version evaluated |
|---|---|
| ModelAngelo | 1.0.18 |
| CryoZeta | commit `c544c6bb` |
| CryoAtom | v1.0.0 |
| EModelX(+AF) | commit `53713cf8` |
| MICA | commit `0bddecdb` |

Grading follows the *2026 Cryo-EM Model-Building Benchmark: Evaluation Plan* (University of
Missouri, version of 27 September 2026), with the additions described under
[Grading rules](#grading-rules).

## Scores

| Score | Program | What it measures |
|---|---|---|
| **TM-score** | US-align | Overall fold agreement after superposition (0 to 1), normalised by the answer key; all chains taken together, fragment order ignored (`-mm 5 -ter 1 -split 0`) |
| **Completeness** | residue script | Answer-key residues with a built C-alpha within 3 Å **and** the correct amino acid, % |
| **Backbone precision** | residue script | Built residues with a C-alpha within 3 Å of an answer-key C-alpha, % |
| **Backbone RMSD** | residue script | RMSD of N, CA, C, O over the matched residues, Å |
| Recall, amino-acid accuracy, C-alpha RMSD, F1 at 1/2/3 Å, chain breaks | residue script | See the docstring of `evaluation/grade_models.py` |
| Aligned length, sequence identity | US-align | Over the TM-score alignment |
| lDDT, backbone lDDT, QS-score, DockQ, ICS, IPS | OpenStructure | Local accuracy and chain interfaces |
| chain_comparison, CC_mask, MolProbity | Phenix | C-alpha match, model-to-map fit, geometry |

`grades.tsv` also reports `TM_ref_mm1`, the TM-score with whole chains paired one to one
(`-mm 1 -ter 0`), for comparison only. For heavily fragmented models, US-align's chain pairing in
that mode can give different values on repeated runs of the same files; the headline TM-score
(`-mm 5`) is stable.

The four scores in bold are the primary scores. All residue scores are computed **in place**:
the model and the answer key share the map's coordinate frame, so nothing is superposed. Built and
answer-key C-alpha atoms are paired one to one by increasing distance, up to 3 Å.

## Grading rules

- **Final model only.** Each tool is graded on the single final model it is documented to deliver.
  The raw (ModelAngelo, CryoAtom) and unrefined (EModelX, MICA) outputs are graded as well, but
  reported separately.
- **Map start offset (26DK).** The 26DK map grid starts at (-80, -80, -80) voxels. ModelAngelo
  ignores this header field, so its models are 132.8 Å off on each axis.
  `fix_map_origin.py` translates them back into the map frame before grading.
- **Maps with more copies than the answer key (26DK, 9UBC, 9I8C).** Here the deposited answer key
  holds one copy of the protein, but the map holds about 47, 2 and 12 copies. These targets are
  also graded per copy (`grade_symmetry_copies.py`): completeness and RMSD against the copy the
  tool built best, and precision as the share of built C-alpha atoms inside the map density at the
  recommended EMDB contour level.
- **Misplaced models.** A model with fewer than half of its C-alpha atoms in the density scores 0
  for completeness, recall, precision and F1 (this applies to MICA on 26DK).
- **Headline averages.** These cover the 14 maps common to all tools: 26DK is left out for every
  tool. Averages over all 15 maps are reported as well.
- **Precision over visible residues.** In the headline, precision counts only built residues
  whose sequence position is modelled in the answer key. CryoZeta, for example, builds the full
  sequence, including disordered parts that the deposited model leaves out.
- **Protein only.** No tool was asked to build the three DNA chains of 36PP. The residue
  scores use amino acids only, and OpenStructure is given protein-only answer keys.
- **Fragmented models and OpenStructure.** OpenStructure pairs chains one to one, so
  `merge_fragments.py` first joins the fragments that sit on the same answer-key chain. This
  step renumbers residues but moves no atoms.

## Requirements

- Python 3 (tested with 3.13) with the packages in `requirements.txt` (gemmi, NumPy, SciPy):
  `pip install -r requirements.txt`
- [US-align](https://github.com/pylelab/USalign) (version 20260920 was used), compiled with
  `g++ -O3 -ffast-math -std=c++17 -o USalign USalign.cpp`
- [OpenStructure](https://openstructure.org) 2.12.0, for example in a conda environment:
  `conda create -n ost -c conda-forge -c bioconda openstructure=2.12.0`
- [Phenix](https://phenix-online.org) 1.20.1-4487

## Data

Download the Zenodo record and unpack it into one directory (`data/` below). The scripts expect
this layout:

```
data/
  maps/emd_<EMDB>.map.gz                      EMDB maps given to the tools
  reference_models/<PDB>_deposited_protein.cif   answer keys (deposited models, no ligands or waters)
  reference_models/9I8C_assembly_protein.cif     12-copy biological assembly of 9I8C
  predicted_models/<tool>/<PDB>/...              models built by the five tools
```

The 15 targets (PDB and EMDB identifiers, resolution, recommended contour level) are listed in
`config/targets.tsv`. The answer keys can also be rebuilt from the PDB with
`evaluation/make_references.py`. Entries can be revised after release, so the Zenodo copies are
the exact files used for the published grades.

## Running the evaluation

Run the whole pipeline with:

```bash
DATA=data USALIGN=/path/to/USalign PHENIX_ENV=/path/to/phenix/phenix_env.sh ./run_all.sh
```

Or run the steps one by one. The scripts are in `evaluation/`; each has `--help`, and `--targets`
limits a run to some targets.

| Step | Script | Output |
|---|---|---|
| 1 | `fix_map_origin.py --data data` | `*_frame_fixed.cif` next to ModelAngelo's 26DK models |
| 2 | `check_references_in_map.py --data data` | `results/reference_check.tsv` |
| 3 | `grade_models.py --data data --usalign USalign` | `results/grades.tsv` |
| 4 | `grade_symmetry_copies.py --data data` | `results/symmetry_copy_grades.tsv` |
| 5 | `run_openstructure.py --data data` | `results/openstructure_scores.tsv` |
| 6 | `run_phenix.py --data data --phenix-env phenix_env.sh` | `results/phenix_scores.tsv` |
| 7 | `make_report_tables.py` | `results/report_tables.json`, `results/report_tables.md` |

Steps 3, 5 and 6 take the longest; most of the time goes to US-align on the many-copy models of
26DK and 9I8C. Steps 5 and 6 can be interrupted: on restart they skip the rows already written.

The output tables correspond to the files in `grades/` on Zenodo: `grades.tsv` to
`pilot15_grades.tsv`, `symmetry_copy_grades.tsv` to `pilot15_symmetry_copy_grades.tsv`,
`openstructure_scores.tsv` to `pilot15_ost.tsv`, `phenix_scores.tsv` to `pilot15_phenix.tsv`
and `report_tables.*` to `pilot15_report_tables.*`. The column names are the same in both.

## Expected results

Headline (`results/report_tables.md`): primary scores, mean over the 14 maps.

| Tool | TM-score | Completeness % | Backbone precision % | Backbone RMSD (Å) |
|---|---|---|---|---|
| ModelAngelo | 0.902 | 88.9 | 99.2 | 0.54 |
| CryoZeta | 0.984 | 96.1 | 96.4 | 0.68 |
| CryoAtom v1 | 0.954 | 91.8 | 98.1 | 0.53 |
| EModelX(+AF) | 0.955 | 92.1 | 97.9 | 0.87 |
| MICA | 0.956 | 92.4 | 98.0 | 0.81 |

## Repository layout

```
config/targets.tsv                     the 15 targets
evaluation/common.py                   shared definitions: targets, tools, file layout, residue helpers
evaluation/make_references.py          answer keys from the PDB (optional, included on Zenodo)
evaluation/fix_map_origin.py           map start-offset correction (26DK, ModelAngelo)
evaluation/check_references_in_map.py  answer keys lie in their maps' density
evaluation/grade_models.py             US-align and residue scores
evaluation/grade_symmetry_copies.py    per-copy grading of 26DK, 9UBC, 9I8C
evaluation/merge_fragments.py          fragment joining before OpenStructure
evaluation/run_openstructure.py        lDDT, QS-score, DockQ
evaluation/run_phenix.py               chain_comparison, map_model_cc, MolProbity
evaluation/make_report_tables.py       all report tables from the grade files
run_all.sh                             the whole pipeline
```

## License

The code is released under the MIT License (see `LICENSE`). The data on Zenodo carry their own
licences: maps and deposited models from EMDB and wwPDB under CC0, predicted models and grades
under CC BY 4.0. If you use the scores, please also cite US-align, OpenStructure, Phenix and gemmi.
