
art_make( BASENAME_ONLY
          LIBRARY_NAME     HitFinderDUNE
          LIB_LIBRARIES    lardataobj_RecoBase
                           lardataalg_DetectorInfo
	  		   larreco_RecoAlg
			   larreco_RecoAlg_PMAlg
			   larevt_Filters
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_RecoObjects
			   larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
			   nusimdata_SimulationBase
                           ART_FRAMEWORK_CORE
   			   ART_FRAMEWORK_PRINCIPAL
			   ART_FRAMEWORK_SERVICES_REGISTRY
			   ART_ROOT_IO_TFILE_SUPPORT
			   ROOT::Core
			   ART_ROOT_IO_TFILESERVICE_SERVICE
                           ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
			   art_Persistency_Common canvas
		       	   art_Persistency_Provenance canvas
			   art_Utilities canvas
			   messagefacility::MF_MessageLogger
			   
	                   FHICLCPP
			   cetlib cetlib_except
                           CLHEP
			   ${Boost_SYSTEM_LIBRARY}
			   ROOT_BASIC_LIB_LIST
                           ROOT_GEOM
         MODULE_LIBRARIES HitFinderDUNE
                           lardataobj_RecoBase
                           lardata_ArtDataHelper
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
			   larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
			   nurandom_RandomUtils_NuRandomService_service
			   nusimdata_SimulationBase
 			   ART_FRAMEWORK_CORE
			   ART_FRAMEWORK_PRINCIPAL
			   ART_FRAMEWORK_SERVICES_REGISTRY
			   ART_ROOT_IO_TFILE_SUPPORT
			   ROOT::Core
			   ART_ROOT_IO_TFILESERVICE_SERVICE
			   art_Persistency_Common canvas
			   art_Utilities canvas
			   messagefacility::MF_MessageLogger
			   
	                   FHICLCPP
			   cetlib cetlib_except
                           ROOT_BASIC_LIB_LIST
                           ROOT_GEOM
         )

install_headers()
install_fhicl()
install_source()
