art_make(
    LIBRARY_NAME VLNEvaluators
    MODULE_LIBRARIES
        VLNData
        VLNModels
        VLNVarExtractors
        VLNDataGenerators
    BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

