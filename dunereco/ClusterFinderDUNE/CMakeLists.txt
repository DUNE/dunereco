
#art_make( 
#	  LIBRARY_NAME  dune_ClusterFinderDUNE          
#          LIB_LIBRARIES larcorealg_Geometry
#	  		dune_RecoAlgDUNE
#                        lardataobj_RecoBase
#                        larreco_RecoAlg_ClusterRecoUtil
#          MODULE_LIBRARIES dune_RecoAlgDUNE_Cluster3DAlgs
#                           larreco_RecoAlg_ClusterRecoUtil			   
#                           larreco_ClusterFinder
#                           larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
#                           larevt_Filters
#                           lardataobj_RecoBase
#                           larcorealg_Geometry
#                           larcore_Geometry_Geometry_service
#                           nurandom_RandomUtils_NuRandomService_service
#                           nusimdata_SimulationBase
#                           ${ART_FRAMEWORK_CORE}
#                           ${ART_FRAMEWORK_PRINCIPAL}
#                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                           ${ART_FRAMEWORK_SERVICES_BASIC}
#                           ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
#                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#                           art_Persistency_Common canvas
#                           art_Utilities canvas
#                           ${messagefacility::MF_MessageLogger}
#                           
#                           ${FHICLCPP}
#                           cetlib cetlib_except
#         )

install_headers()
install_fhicl()
install_source()


