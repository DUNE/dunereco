include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  IniSegAlg
          LIB_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
			larreco_RecoAlg
			larreco_RecoAlg_PMAlg
                        lardata_Utilities
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        art_Persistency_Common canvas
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${FHICLCPP}
			${MF_MESSAGELOGGER}
			
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
          )


install_headers()
install_fhicl()
install_source()
