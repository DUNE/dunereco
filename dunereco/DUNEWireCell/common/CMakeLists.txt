file(GLOB jsonnet_file *.jsonnet )

message(STATUS "jsonnet files ${jsonnet_file}")

install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/common)

file(GLOB jsonnet_file sim/*.jsonnet )

message(STATUS "jsonnet files ${jsonnet_file}")

install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/common/sim)
