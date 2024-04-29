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
  LIBRARY_NAME      TransformerCVNArt
  EXCLUDE ${EXCLUDE_TF}
  LIB_LIBRARIES
      RegCNNFunc
      TransformerCVNFunc
  ROOT::Core
  canvas::canvas
  art::Persistency_Provenance
  art::Utilities
  torch
  torch_cpu
  c10
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except::cetlib_except
  lardataobj::RecoBase
  lardata::Utilities
  larpandora::LArPandoraInterface
  nusimdata::SimulationBase
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  MVAAlg
  larreco::RecoAlg
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
 
  MODULE_LIBRARIES
      RegCNNFunc
      TransformerCVNFunc
  TransformerCVNArt
  )

# use add or remove
##cet_set_compiler_flags(DIAGS CAUTIOUS
##    NO_UNDEFINED
##    EXTRA_FLAGS
##    -pedantic)

install_headers()
install_fhicl()
install_source()
