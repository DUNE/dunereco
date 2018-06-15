# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

simple_plugin(SNSlicer "module"
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_Simulation
                        lardataobj_RawData
			lardataobj_RecoBase
			nusimdata_SimulationBase
                        nutools_ParticleNavigation
                        lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
)

cet_make( LIBRARIES lardataobj_RawData
                    canvas_Persistency_Common
                    canvas_Persistency_Provenance
                    canvas_Utilities
                    MF_MessageLogger
                    
                    cetlib
                    ${ROOT_BASIC_LIB_LIST}
                    NO_DICTIONARY )


art_dictionary( DICTIONARY_LIBRARIES ) # dune_SNSlicer )

install_headers()
install_fhicl()
install_source()
install_scripts()
