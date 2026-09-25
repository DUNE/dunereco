#!/bin/bash
# n-nbar TransformerCVN inference chain on RecoEnergyS output files.
#
#   run_nnbar_inference.sh -i INPUT -o OUTDIR [-t TYPE] [-j NPAR] [-m MODEL_TAG]
#                          [-w WEIGHTS -O OPTIONS.json] [-s STOP_AFTER]
#
#   INPUT      a RecoEnergyS output file (*_cvnpreprocess.root, tree recoEnergy/WC), a directory
#              of them, or a text file listing them (one per line; /pnfs paths are fine)
#   OUTDIR     receives pixelmap/<name>.root (stage 5), pixelmap.h5 (stage 6),
#              sparse.h5 (stage 7) and predictions.h5 (stage 8)
#   TYPE       image geometry of the stage-5 macro: nnbar (default) or atmnu
#   NPAR       parallel stage-5 jobs (default 1)
#   MODEL_TAG  sample tag stored in the h5 "model" column (atm, ha_br, ...); default none
#   WEIGHTS    checkpoint of the trained network; OPTIONS.json its options file (stage 8)
#   STOP_AFTER pixelmap | h5 | sparse : stop after that stage
#
# Stages 7 and 8 need sparsify.py and evaluate.py from the TransformerCVN training toolkit
# next to this script; until they are present the chain stops after stage 6 and says so.
# Requires: source setup_nnbar_cvn.sh
set -uo pipefail
_here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MACRO="$_here/make_text_file_to_root_trks_shws.C"
INPUT=""; OUTDIR=""; TYPE=nnbar; NPAR=1; MODEL_TAG=""; WEIGHTS=""; OPTIONS=""; STOP=""
usage() { sed -n '2,20p' "$0"; exit 1; }
while getopts "i:o:t:j:m:w:O:s:h" opt; do
  case $opt in i) INPUT=$OPTARG;; o) OUTDIR=$OPTARG;; t) TYPE=$OPTARG;; j) NPAR=$OPTARG;; m) MODEL_TAG=$OPTARG;;
    w) WEIGHTS=$OPTARG;; O) OPTIONS=$OPTARG;; s) STOP=$OPTARG;; *) usage;; esac
done
[ -n "$INPUT" ] && [ -n "$OUTDIR" ] || usage
command -v root >/dev/null || { echo "root not in PATH; source setup_nnbar_cvn.sh"; exit 1; }
python -c "import tables, uproot" 2>/dev/null || { echo "python env missing; source setup_nnbar_cvn.sh"; exit 1; }

# ---- input list ---------------------------------------------------------------------------
mkdir -p "$OUTDIR/pixelmap"
LIST="$OUTDIR/inputs.txt"
if   [ -d "$INPUT" ]; then find "$INPUT" -maxdepth 1 -type f -name '*.root' -size +0c | LC_ALL=C sort > "$LIST"
elif [ -f "$INPUT" ] && [[ "$INPUT" != *.root ]]; then grep -v '^\s*$' "$INPUT" > "$LIST"
elif [ -f "$INPUT" ]; then printf '%s\n' "$INPUT" > "$LIST"
else echo "input $INPUT not found"; exit 1; fi
n=$(wc -l < "$LIST"); [ "$n" -gt 0 ] || { echo "no input files"; exit 1; }
echo "== stage 5: $n RecoEnergyS file(s) -> $OUTDIR/pixelmap (type=$TYPE, $NPAR parallel)"

# /pnfs inputs go through xrootd when a bearer token is available, else the NFS mount.
have_token() { [ -n "${BEARER_TOKEN:-}" ] || [ -r "${BEARER_TOKEN_FILE:-/run/user/$(id -u)/bt_u$(id -u)}" ]; }
pnfs_to_xrootd() {
  if have_token; then printf '%s\n' "${1/\/pnfs\/dune/root:\/\/fndcadoor.fnal.gov:1094\/pnfs\/fnal.gov\/usr\/dune}"; else printf '%s\n' "$1"; fi
}
process_one() {  # process_one <input> <output.root>   (outputs above 10 kB are kept: resumable)
  local f=$1 o=$2 base; base=$(basename "$f")
  if [ -f "$o" ] && (( $(wc -c <"$o") > 10240 )); then echo "skip $base (done)"; return 0; fi
  echo "process $base"
  root -l -b -q "${MACRO}(\"$(pnfs_to_xrootd "$f")\", \"$o\", \"$TYPE\")" > "${o%.root}.log" 2>&1 \
    || echo "FAILED $base (see ${o%.root}.log)"
}
export -f process_one pnfs_to_xrootd have_token; export MACRO TYPE
while IFS= read -r f; do [ -n "$f" ] && printf '%s\0%s\0' "$f" "$OUTDIR/pixelmap/$(basename "${f%.*}").root"; done < "$LIST" \
  | xargs -0 -n 2 -P "$NPAR" bash -c 'process_one "$0" "$1"'
[ "$STOP" = pixelmap ] && exit 0

# ---- stage 6 ------------------------------------------------------------------------------
PMLIST="$OUTDIR/pixelmap.txt"
find "$OUTDIR/pixelmap" -maxdepth 1 -type f -name '*.root' -size +10000c | LC_ALL=C sort > "$PMLIST"
echo "== stage 6: $(wc -l < "$PMLIST") pixelmap file(s) -> $OUTDIR/pixelmap.h5"
python "$_here/preprocess.py" "$OUTDIR/pixelmap.h5" --filelist "$PMLIST" --no-stage ${MODEL_TAG:+--model-tag "$MODEL_TAG"} || exit 1
[ "$STOP" = h5 ] && exit 0

# ---- stage 7: sparsify (TransformerCVN training toolkit) ------------------------------------
if [ ! -f "$_here/sparsify.py" ]; then
  echo "== stage 7: sparsify.py not available in $_here; chain stops after stage 6 ($OUTDIR/pixelmap.h5)"; exit 0
fi
echo "== stage 7: $OUTDIR/pixelmap.h5 -> $OUTDIR/sparse.h5"
python "$_here/sparsify.py" "$OUTDIR/pixelmap.h5" "$OUTDIR/sparse.h5" || exit 1
[ "$STOP" = sparse ] && exit 0

# ---- stage 8: evaluate (TransformerCVN training toolkit) ------------------------------------
if [ ! -f "$_here/evaluate.py" ]; then
  echo "== stage 8: evaluate.py not available in $_here; chain stops after stage 7 ($OUTDIR/sparse.h5)"; exit 0
fi
[ -n "$WEIGHTS" ] && [ -n "$OPTIONS" ] || { echo "stage 8 needs -w WEIGHTS and -O OPTIONS.json"; exit 1; }
echo "== stage 8: evaluating $OUTDIR/sparse.h5 with $WEIGHTS -> $OUTDIR/predictions.h5"
python "$_here/evaluate.py" --options "$OPTIONS" --checkpoint "$WEIGHTS" --training-file "$OUTDIR/sparse.h5" \
  --split testing --testing-file "$OUTDIR/sparse.h5" --output "$OUTDIR/predictions.h5" --device cpu
