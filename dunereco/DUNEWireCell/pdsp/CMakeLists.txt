
# install txt files
set( ${product}_wirecell_dir "${product}/${version}/wire-cell-cfg/pgrapher/experiment/pdsp/" )
file(GLOB jsonnet_file *.jsonnet )
install( FILES ${jsonnet_file} DESTINATION ${${product}_wirecell_dir} )
