# install files
file(GLOB jsonnet_file *.jsonnet )
message(STATUS "jsonnet files ${jsonnet_file}")
install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/experiment/protodunevd)

file(GLOB jsonbz2_file *.json.bz2)
install_wp(LIST ${jsonbz2_file})

# TorchScript DNN models, installed to the WIRECELL_PATH sub-paths the
# protodunevd jsonnet refers to (dnnroi_model = 'dnnroi/pdvd/...',
# l1sp_pd_dnn_model = 'l1sp/pdvd/...' in wcls-nf-sp-dnnroi.jsonnet).
# Source files are kept flat in this dir; install_wp's SUBDIRNAME places them.
#install_wp(LIST pipe_distill_nestedunet_6ch.ts SUBDIRNAME dnnroi/pdvd)
#install_wp(LIST l1sp_dnn_pdvd_v1.ts SUBDIRNAME l1sp/pdvd)


