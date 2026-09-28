# n-nbar TransformerCVN: pixel maps, preprocessing, sparsify, training and inference.
# Python and ROOT scripts only, nothing is compiled. The two entry points go into the product's
# bin directory; the whole tree is installed with the sources, where the entry points find it
# ($DUNERECO_DIR/source/dunereco/TransformerCVN/nnbar, see setup/setup_nnbar_cvn.sh).
install_scripts(LIST setup/setup_nnbar_cvn.sh inference/run_nnbar_inference.sh inference/summarize_scores.py)
install(DIRECTORY setup pixelmap preprocess sparsify training inference
        DESTINATION source/dunereco/TransformerCVN/nnbar
        USE_SOURCE_PERMISSIONS
        PATTERN "__pycache__" EXCLUDE
        PATTERN "*.ckpt" EXCLUDE)
