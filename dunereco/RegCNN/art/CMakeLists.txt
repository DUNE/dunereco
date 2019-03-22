
art_make(BASENAME_ONLY
  LIBRARY_NAME      RegCNNArt
  EXCLUDE ${EXCLUDE_TF}
  LIB_LIBRARIES RegCNNFunc
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${ART_PERSISTENCY_COMMON}
  canvas
  art_Persistency_Provenance
  ${ART_UTILITIES}

  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  lardataobj_RecoBase
  lardata_Utilities
  nusimdata_SimulationBase
  larsim_MCCheater_BackTrackerService_service           
  larsim_MCCheater_ParticleInventoryService_service
  MVAAlg
  larreco_RecoAlg
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  ${LARRECO_RECOALG_IMAGEPATTERNALGS_TF}
  MODULE_LIBRARIES  RegCNNFunc
  RegCNNArt
  )


install_headers()
install_fhicl()
install_source()
