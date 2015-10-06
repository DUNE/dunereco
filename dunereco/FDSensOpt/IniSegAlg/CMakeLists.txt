include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  	IniSegAlg
          LIB_LIBRARIES 	RecoBase
					RecoAlg
					PMAlg
					DetectorProperties_service
                        	LArProperties_service
                        	Geometry
                        	Geometry_service
                        	Simulation
                        	RawData
                        	Utilities
                        	${SIMULATIONBASE}
                        	${ART_FRAMEWORK_CORE}
                        	${ART_FRAMEWORK_PRINCIPAL}
                        	${ART_PERSISTENCY_COMMON}
                        	${ART_FRAMEWORK_SERVICES_REGISTRY}
                        	${ART_FRAMEWORK_SERVICES_BASIC}
                       	${ART_FRAMEWORK_SERVICES_OPTIONAL}
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
