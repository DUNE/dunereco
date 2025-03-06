
# install txt files
#set( ${product}_wirecell_dir "${product}/${version}/wire-cell-cfg/pdsp/" )
file(GLOB jsonnet_file *.jsonnet )
#install( FILES ${jsonnet_file} DESTINATION ${${product}_wirecell_dir} )

message(STATUS "jsonnet files ${jsonnet_file}")

install_wp(LIST ${jsonnet_file} SUBDIRNAME pdsp)


file(GLOB jsonbz2_file *.json.bz2)
install_wp(LIST ${jsonbz2_file})
