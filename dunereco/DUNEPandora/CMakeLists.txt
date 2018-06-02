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
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			${LARPANDORACONTENT}
                        ${PANDORASDK}
                        ${PANDORAMONITORING}
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
			${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES DUNEPandora
          )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
