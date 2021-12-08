
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
#                           nurandom::RandomUtils_NuRandomService_service
#                           nusimdata::SimulationBase
#                           art::Framework_Core
#                           art::Framework_Principal
#                           art::Framework_Services_Registry
#                           ${ART_FRAMEWORK_SERVICES_BASIC}
#                           art_root_io::tfile_support ROOT::Core
#                           art_root_io::TFileService_service
#                           art::Persistency_Common canvas
#                           art::Utilities canvas
#                           messagefacility::MF_MessageLogger
#                           fhiclcpp::fhiclcpp
#                           cetlib::cetlib cetlib_except
#         )

install_headers()
install_fhicl()
install_source()


