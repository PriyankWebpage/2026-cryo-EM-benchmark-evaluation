#!/usr/bin/env bash
# Runs the whole evaluation in order (see README). Paths are set with environment variables:
#   DATA        unpacked Zenodo record (default: data)
#   OUT         output directory (default: results)
#   USALIGN     US-align executable (default: USalign on the PATH)
#   OST         command that starts OpenStructure (default: conda run -n ost ost)
#   PHENIX_ENV  phenix_env.sh of the Phenix installation (required)
#   PYTHON      Python interpreter (default: python3)
set -euo pipefail

SRC="$(cd "$(dirname "$0")" && pwd)/evaluation"
DATA=${DATA:-data}
OUT=${OUT:-results}
USALIGN=${USALIGN:-USalign}
OST=${OST:-conda run -n ost ost}
PHENIX_ENV=${PHENIX_ENV:?set PHENIX_ENV to the phenix_env.sh of your Phenix installation}
PYTHON=${PYTHON:-python3}

"$PYTHON" "$SRC/fix_map_origin.py" --data "$DATA"
"$PYTHON" "$SRC/check_references_in_map.py" --data "$DATA" --out "$OUT/reference_check.tsv"
"$PYTHON" "$SRC/grade_models.py" --data "$DATA" --usalign "$USALIGN" --out "$OUT/grades.tsv"
"$PYTHON" "$SRC/grade_symmetry_copies.py" --data "$DATA" --out "$OUT/symmetry_copy_grades.tsv"
"$PYTHON" "$SRC/run_openstructure.py" --data "$DATA" --ost "$OST" \
    --work "$OUT/work/openstructure" --out "$OUT/openstructure_scores.tsv"
"$PYTHON" "$SRC/run_phenix.py" --data "$DATA" --phenix-env "$PHENIX_ENV" \
    --work "$OUT/work/phenix" --out "$OUT/phenix_scores.tsv"
"$PYTHON" "$SRC/make_report_tables.py" --grades "$OUT/grades.tsv" \
    --copies "$OUT/symmetry_copy_grades.tsv" --openstructure "$OUT/openstructure_scores.tsv" \
    --phenix "$OUT/phenix_scores.tsv" \
    --out-json "$OUT/report_tables.json" --out-md "$OUT/report_tables.md"
