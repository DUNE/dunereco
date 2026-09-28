#!/bin/bash
# n-nbar TransformerCVN inference chain on RecoEnergyS output files.
#
#   run_nnbar_inference.sh -i INPUT -o OUTDIR [-t TYPE] [-j NPAR] [-m MODEL_TAG]
#                          [-w WEIGHTS] [-O OPTIONS.json] [-s STOP_AFTER] [-P]
#   run_nnbar_inference.sh -6 STAGE6.h5 -o OUTDIR [...]      start from an existing stage-6 file
#
#   INPUT      a RecoEnergyS output file (*_cvnpreprocess.root, tree recoEnergy/WC), a directory
#              of them, or a text file listing them (one per line; /pnfs paths are fine)
#   OUTDIR     receives pixelmap/<name>.root (stage 5), pixelmap.h5 (stage 6),
#              sparse.h5 (stage 7) and predictions.h5 (stage 8)
#   TYPE       image geometry of the stage-5 macro: nnbar (default) or atmnu
#   NPAR       parallel stage-5 jobs (default 1)
#   MODEL_TAG  sample tag stored in the h5 "model" column (atm, ha_br, ...); default none
#   WEIGHTS    checkpoint of the trained network   [$NNBAR_CVN_WEIGHTS, fetched by setup_nnbar_cvn.sh]
#   OPTIONS    options file matching the checkpoint [training/option_files/example.json]
#   STOP_AFTER pixelmap | h5 | sparse : stop after that stage
#   -P         do not apply the analysis precut in stage 7
#
# Stage 7 (../sparsify) and stage 8 (../training/evaluate.py) run on the CPU; stage 8 accepts
# -d cuda:0. Requires: source ../setup/setup_nnbar_cvn.sh (locates the tree, sets NNBAR_CVN_DIR)
set -uo pipefail
_here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
[ -f "$_here/../training/evaluate.py" ] && NNBAR_CVN_DIR="$( cd "$_here/.." && pwd )"
[ -n "${NNBAR_CVN_DIR:-}" ] && [ -f "$NNBAR_CVN_DIR/training/evaluate.py" ] || { echo "nnbar tree not found; source setup_nnbar_cvn.sh"; exit 1; }
MACRO="$NNBAR_CVN_DIR/pixelmap/make_text_file_to_root_trks_shws.C"
INPUT=""; STAGE6=""; OUTDIR=""; TYPE=nnbar; NPAR=1; MODEL_TAG=""; WEIGHTS=""; OPTIONS=""; STOP=""; NOPRECUT=""; DEVICE=cpu
usage() { sed -n '2,24p' "$0"; exit 1; }
while getopts "i:6:o:t:j:m:w:O:s:d:Ph" opt; do
  case $opt in i) INPUT=$OPTARG;; 6) STAGE6=$OPTARG;; o) OUTDIR=$OPTARG;; t) TYPE=$OPTARG;; j) NPAR=$OPTARG;; m) MODEL_TAG=$OPTARG;;
    w) WEIGHTS=$OPTARG;; O) OPTIONS=$OPTARG;; s) STOP=$OPTARG;; d) DEVICE=$OPTARG;; P) NOPRECUT=--no-precut;; *) usage;; esac
done
{ [ -n "$INPUT" ] || [ -n "$STAGE6" ]; } && [ -n "$OUTDIR" ] || usage
python -c "import tables, uproot" 2>/dev/null || { echo "python env missing; source setup_nnbar_cvn.sh"; exit 1; }
mkdir -p "$OUTDIR"
H5="$OUTDIR/pixelmap.h5"

if [ -n "$STAGE6" ]; then
  [ -f "$STAGE6" ] || { echo "stage-6 file $STAGE6 not found"; exit 1; }
  echo "== starting from stage-6 file $STAGE6"; H5="$STAGE6"
else
command -v root >/dev/null || { echo "root not in PATH; source setup_nnbar_cvn.sh"; exit 1; }

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
python "$NNBAR_CVN_DIR/preprocess/preprocess.py" "$H5" --filelist "$PMLIST" --no-stage ${MODEL_TAG:+--model-tag "$MODEL_TAG"} || exit 1
[ "$STOP" = h5 ] && exit 0
fi

# ---- stage 7: sparsify (network input format; applies the analysis precut) -------------------
echo "== stage 7: $H5 -> $OUTDIR/sparse.h5"
python "$NNBAR_CVN_DIR/sparsify/sparsify.py" "$H5" "$OUTDIR/sparse.h5" --jobs "$NPAR" $NOPRECUT || exit 1
[ "$STOP" = sparse ] && exit 0

# ---- stage 8: evaluate with ../training/evaluate.py -------------------------------------------
TOOLKIT="$NNBAR_CVN_DIR/training"
WEIGHTS=${WEIGHTS:-${NNBAR_CVN_WEIGHTS:-}}; OPTIONS=${OPTIONS:-$TOOLKIT/option_files/example.json}
[ -n "$WEIGHTS" ] || { echo "stage 8: no checkpoint (-w WEIGHTS or NNBAR_CVN_WEIGHTS from setup_nnbar_cvn.sh)"; exit 1; }
[ -f "$WEIGHTS" ] && [ -f "$OPTIONS" ] || { echo "stage 8: weights $WEIGHTS or options $OPTIONS not found"; exit 1; }
rm -f "$OUTDIR/predictions.h5"
echo "== stage 8: evaluating $OUTDIR/sparse.h5 with $WEIGHTS on $DEVICE -> $OUTDIR/predictions.h5"
PYTHONPATH="$TOOLKIT${PYTHONPATH:+:$PYTHONPATH}" python "$TOOLKIT/evaluate.py" --options "$OPTIONS" --checkpoint "$WEIGHTS" \
  --training-file "$OUTDIR/sparse.h5" --testing-file "$OUTDIR/sparse.h5" --split testing \
  --output "$OUTDIR/predictions.h5" --device "$DEVICE" || exit 1
echo "== done: $OUTDIR/predictions.h5 (one row per real event of sparse.h5, in order; the padding event is not evaluated)"
