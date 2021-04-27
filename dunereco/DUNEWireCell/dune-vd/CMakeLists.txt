
# install txt files
file(GLOB jsonnet_file *.jsonnet )
message(STATUS "jsonnet files ${jsonnet_file}")
install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/experiment/dune-vd)

file(GLOB json_file *.json )
message(STATUS "json files ${json_file}")
install_wp(LIST ${json_file} SUBDIRNAME pgrapher/experiment/dune-vd)

file(GLOB jsonbz2_file *.json.bz2)
install_wp(LIST ${jsonbz2_file})
