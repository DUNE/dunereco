if ( DEFINED ENV{LIBTORCH_DIR})
   add_subdirectory(art)
   add_subdirectory(func)
endif()
# nnbar/ is a standalone tree (python and ROOT scripts of the n-nbar study); it is not built or installed.
# RecoEnergyStudies: RecoEnergyS analyzer (n-nbar CVN preparation, first stage of nnbar/); does not need libtorch.
add_subdirectory(RecoEnergyStudies)
