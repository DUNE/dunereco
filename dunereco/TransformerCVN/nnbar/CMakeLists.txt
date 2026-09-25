# n-nbar TransformerCVN inference chain (stages 5-8): scripts only, nothing is compiled.
# They are installed into the product's bin directory so that, after "setup dunereco",
# run_nnbar_inference.sh and setup_nnbar_cvn.sh are on the PATH.
install_scripts(LIST
  run_nnbar_inference.sh
  setup_nnbar_cvn.sh
  preprocess.py
  make_text_file_to_root_trks_shws.C
  requirements.txt
  requirements-eval.txt
)
install_source()
