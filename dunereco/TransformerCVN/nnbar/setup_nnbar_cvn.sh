#!/bin/bash
# Environment for the n-nbar TransformerCVN inference chain (stages 5-8): ROOT for the
# stage-5 macro and a python venv for stages 6-8.
#
#   source setup_nnbar_cvn.sh [--no-eval]
#
# ROOT: taken from the current environment (e.g. after "setup dunesw ..." or "setup root ...");
#       if none is found, ROOT 6.28.12 is loaded from the larsoft spack-packages area on cvmfs.
# venv: $NNBAR_CVN_VENV, default /exp/dune/app/users/$USER/nnbar-cvn-venv when that area
#       exists (home directories on the gpvms are small), else $HOME/nnbar-cvn-venv.
#       Created on first use from requirements.txt and, unless --no-eval, requirements-eval.txt.
_here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
_with_eval=1; [ "${1:-}" = "--no-eval" ] && _with_eval=0

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
    && "$NNBAR_CVN_VENV/bin/pip" install -q -r "$_here/requirements.txt" \
    && { [ $_with_eval -eq 0 ] || "$NNBAR_CVN_VENV/bin/pip" install -q -r "$_here/requirements-eval.txt"; } \
    || echo "WARNING: venv creation failed"
fi
source "$NNBAR_CVN_VENV/bin/activate"
export NNBAR_CVN_DIR="$_here"
echo "nnbar CVN env: root $(root-config --version 2>/dev/null || echo none), python $(python --version 2>&1), venv $NNBAR_CVN_VENV"
unset _here _with_eval
