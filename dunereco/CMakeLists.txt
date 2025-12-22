add_subdirectory(AnaUtils)
add_subdirectory(Cluster3D)
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

# Supera
#
# Why are we checking explicitly for git submodule initialization?
# ----------------------------------------------------------------
# Because it is stored as a git submodule, it should be initialized first before
# any build is attempted. `mrb gitCheckout` should run `--recurse-submodules` to
# this effect. However, it is implemented in a hacky way: the option for submodules
# is only added to `mrb gitCheckout` command if we are cloning from SBNSoftware/sbncode
# and nothing else. This means cloning from SBNSoftware/sbncode#somePRnumber will
# NOT carry this option. This type of cloning can happen in CI builds for example.
# In other words, we need to fix the hack with another hack: if Supera does not
# look initialized, we take the matter in our hands and run the git submodule
# initialization command.
#
# (This explanation above might be relevant for any other git submodules.)
if ( NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Supera/.git" )
  execute_process(COMMAND git submodule update --init --recursive
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# We need to run Supera/setup.sh first, otherwise it won't have a CMakeLists.txt.
execute_process(COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/Supera/setup.sh dune
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Supera )
                
add_subdirectory(Supera)
# these two were not built in dunetpc at the time of the split
# add_subdirectory(SNSlicer)
# add_subdirectory(SNUtils)
