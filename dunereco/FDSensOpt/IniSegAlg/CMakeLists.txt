include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  IniSegAlg
          LIB_LIBRARIES lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
                        lardata_DetectorInfoServices_LArPropertiesServiceStandard_service
                        Geometry
                        larcore_Geometry_Geometry_service
                        Simulation
                        RawData
                        RecoBase
			RecoAlg
			PMAlg
                        Utilities
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${FHICLCPP}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
          )


install_headers()
install_fhicl()
install_source()
