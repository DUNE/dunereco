art_make(
    LIBRARY_NAME VLNDataGenerators
    MODULE_LIBRARIES
        VLNData
        VLNVarExtractors
        VLNDataGenerators
    BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
