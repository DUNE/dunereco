file(GLOB jsonnet_file *.jsonnet )

message(STATUS "jsonnet files ${jsonnet_file}")

install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/common)

