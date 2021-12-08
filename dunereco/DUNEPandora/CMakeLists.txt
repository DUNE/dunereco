# where should the scripts/..xml file be installed?  Perhaps in bin?

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make( BASENAME_ONLY
          LIBRARY_NAME  DUNEPandora
          LIB_LIBRARIES larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			LARPANDORACONTENT
                        PANDORASDK
                        PANDORAMONITORING
                        nusimdata::SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
                        ROOT_GEOM
			
          MODULE_LIBRARIES DUNEPandora
          )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
