
# install txt files
file(GLOB jsonnet_file *.jsonnet )

message(STATUS "jsonnet files ${jsonnet_file}")

install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/experiment/dune10kt-1x2x2)


file(GLOB jsonbz2_file *.json.bz2)
install_wp(LIST ${jsonbz2_file})
