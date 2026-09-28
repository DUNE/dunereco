#!/bin/bash
# Environment for the n-nbar TransformerCVN chain: ROOT for the stage-5 macro, a python venv for
# stages 6-8 and for training, and the trained checkpoint.
#
#   source setup_nnbar_cvn.sh [--no-eval] [--train]
#
# ROOT:    taken from the current environment (e.g. after "setup dunesw ..." or "setup root ...");
#          if none is found, ROOT 6.28.12 is loaded from the larsoft spack-packages area on cvmfs.
# venv:    $NNBAR_CVN_VENV, default /exp/dune/app/users/$USER/nnbar-cvn-venv when that area exists
#          (home directories on the gpvms are small), else $HOME/nnbar-cvn-venv. Created on first use
#          from requirements.txt plus requirements-eval.txt (CPU torch; skipped with --no-eval) or,
#          with --train, requirements-train.txt after a CUDA torch (see that file).
# code:    the network code, train.py and evaluate.py live in ../training of this tree; it is put on
#          PYTHONPATH. The tree is located from this script, or from $NNBAR_CVN_DIR, or from the
#          installed sources of dunereco ($DUNERECO_DIR/source/dunereco/TransformerCVN/nnbar).
# weights: $NNBAR_CVN_WEIGHTS (default: nnbar_best.ckpt next to the venv); downloaded once from the
#          author's repository at the pinned commit $NNBAR_CVN_WEIGHTS_COMMIT if missing.
_here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
_with_eval=1; _train=0
for _a in "$@"; do case $_a in --no-eval) _with_eval=0;; --train) _train=1;; esac; done

if   [ -f "$_here/../training/evaluate.py" ]; then NNBAR_CVN_DIR="$( cd "$_here/.." && pwd )"
elif [ -n "${NNBAR_CVN_DIR:-}" ] && [ -f "$NNBAR_CVN_DIR/training/evaluate.py" ]; then :
elif [ -n "${DUNERECO_DIR:-}" ] && [ -f "$DUNERECO_DIR/source/dunereco/TransformerCVN/nnbar/training/evaluate.py" ]; then
  NNBAR_CVN_DIR="$DUNERECO_DIR/source/dunereco/TransformerCVN/nnbar"
else echo "ERROR: cannot locate the nnbar tree (set NNBAR_CVN_DIR)"; return 1 2>/dev/null || exit 1; fi
export NNBAR_CVN_DIR

if ! command -v root >/dev/null 2>&1; then
  # The spack-packages instance has a single root@6.28.12 (AL9, gcc 12); the other instances
  # on cvmfs carry several builds of the same version and need a hash-qualified spec.
  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh 2>/dev/null
  spack load root@6.28.12 2>/dev/null || echo "WARNING: no ROOT in PATH and spack load failed; set up ROOT by hand"
fi

if [ -z "${NNBAR_CVN_VENV:-}" ]; then
  if [ -d "/exp/dune/app/users/$USER" ]; then NNBAR_CVN_VENV=/exp/dune/app/users/$USER/nnbar-cvn-venv
  else NNBAR_CVN_VENV=$HOME/nnbar-cvn-venv; fi
fi
export NNBAR_CVN_VENV
if [ ! -x "$NNBAR_CVN_VENV/bin/python" ]; then
  echo "creating venv $NNBAR_CVN_VENV"
  python3 -m venv "$NNBAR_CVN_VENV" && "$NNBAR_CVN_VENV/bin/pip" install -q --upgrade pip \
    && "$NNBAR_CVN_VENV/bin/pip" install -q -r "$NNBAR_CVN_DIR/setup/requirements.txt" \
    && { if [ $_train -eq 1 ]; then
           "$NNBAR_CVN_VENV/bin/pip" install -q torch==2.0.1 torchvision==0.15.2 --index-url https://download.pytorch.org/whl/cu118 \
             && "$NNBAR_CVN_VENV/bin/pip" install -q -r "$NNBAR_CVN_DIR/setup/requirements-train.txt"
         elif [ $_with_eval -eq 1 ]; then
           "$NNBAR_CVN_VENV/bin/pip" install -q -r "$NNBAR_CVN_DIR/setup/requirements-eval.txt"
         fi; } \
    || echo "WARNING: venv creation failed"
fi
source "$NNBAR_CVN_VENV/bin/activate"

# Network code (transformercvn package), train.py and evaluate.py.
export NNBAR_CVN_TOOLKIT="$NNBAR_CVN_DIR/training"
case ":${PYTHONPATH:-}:" in *":$NNBAR_CVN_TOOLKIT:"*) ;; *) export PYTHONPATH="$NNBAR_CVN_TOOLKIT${PYTHONPATH:+:$PYTHONPATH}" ;; esac

# Trained checkpoint (70 MB, not kept in dunereco): best epoch of the 13 July 2026 training.
NNBAR_CVN_WEIGHTS_COMMIT=${NNBAR_CVN_WEIGHTS_COMMIT:-deb1014}
NNBAR_CVN_WEIGHTS_URL=${NNBAR_CVN_WEIGHTS_URL:-https://github.com/KaiwenYu2001/dune-nnbar-transformercvn_v2/raw/$NNBAR_CVN_WEIGHTS_COMMIT/nnbar_best.ckpt}
export NNBAR_CVN_WEIGHTS=${NNBAR_CVN_WEIGHTS:-$(dirname "$NNBAR_CVN_VENV")/nnbar-cvn-weights/nnbar_best.ckpt}
if [ ! -s "$NNBAR_CVN_WEIGHTS" ] && [ $_with_eval -eq 1 ]; then
  echo "downloading $NNBAR_CVN_WEIGHTS_URL -> $NNBAR_CVN_WEIGHTS"
  mkdir -p "$(dirname "$NNBAR_CVN_WEIGHTS")" && curl -fsSL -o "$NNBAR_CVN_WEIGHTS" "$NNBAR_CVN_WEIGHTS_URL" \
    || { rm -f "$NNBAR_CVN_WEIGHTS"; echo "WARNING: checkpoint download failed; stage 8 needs -w WEIGHTS"; }
fi
echo "nnbar CVN env: root $(root-config --version 2>/dev/null || echo none), python $(python --version 2>&1), venv $NNBAR_CVN_VENV, tree $NNBAR_CVN_DIR"
unset _here _with_eval _train _a
