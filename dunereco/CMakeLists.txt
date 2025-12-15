add_subdirectory(AnaUtils)
add_subdirectory(ClusterFinderDUNE)
add_subdirectory(CVN)
add_subdirectory(DUNEPandora)
add_subdirectory(DUNEWireCell)
add_subdirectory(FDSensOpt)
add_subdirectory(HitFinderDUNE)
add_subdirectory(InfillChannels)
add_subdirectory(LowEUtils)
add_subdirectory(FDSelections)
add_subdirectory(RecoAlgDUNE)
add_subdirectory(RegCNN)
add_subdirectory(TrackFinderDUNE)
add_subdirectory(TrackPID)
add_subdirectory(TransformerCVN)
add_subdirectory(VLNets)
add_subdirectory(PointResTree)
add_subdirectory(WireFilter)

if ( NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Supera/.git" )
  execute_process(COMMAND git submodule update --init --recursive
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# We need to run Supera/setup.sh first, otherwise it won't have a CMakeLists.txt.
execute_process(COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/Supera/setup.sh pdvd
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Supera )
add_subdirectory(Supera)

# these two were not built in dunetpc at the time of the split
# add_subdirectory(SNSlicer)
# add_subdirectory(SNUtils)
