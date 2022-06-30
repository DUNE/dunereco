function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
        endif()
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
endfunction()

message(STATUS "Torch include: ${TORCH_INCLUDE_DIRS}")
dump_cmake_variables("TORCH")

art_make(BASENAME_ONLY
  LIBRARY_NAME      RegCNNArt
  EXCLUDE ${EXCLUDE_TF}
  LIB_LIBRARIES RegCNNFunc
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  ART_PERSISTENCY_COMMON
  canvas::canvas
  art::Persistency_Provenance
  art::Utilities
  torch
  torch_cpu
  c10
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  lardataobj_RecoBase
  lardata_Utilities
  larpandora_LArPandoraInterface
  nusimdata::SimulationBase
  larsim_MCCheater_BackTrackerService_service           
  larsim_MCCheater_ParticleInventoryService_service
  MVAAlg
  larreco_RecoAlg
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  MODULE_LIBRARIES  RegCNNFunc
  RegCNNArt
  )

# use add or remove
##cet_set_compiler_flags(DIAGS CAUTIOUS
##    NO_UNDEFINED
##    EXTRA_FLAGS
##    -pedantic)

install_headers()
install_fhicl()
install_source()
