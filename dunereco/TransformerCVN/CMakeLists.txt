if ( DEFINED ENV{LIBTORCH_DIR})
   add_subdirectory(art)
   add_subdirectory(func)
endif()
add_subdirectory(nnbar)   # scripts only, no libtorch needed
# RecoEnergyStudies: RecoEnergyS analyzer (n-nbar CVN preparation, first stage of nnbar/); does not need libtorch.
add_subdirectory(RecoEnergyStudies)
