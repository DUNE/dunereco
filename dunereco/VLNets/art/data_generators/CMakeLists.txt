art_make(
    LIBRARY_NAME VLNDataGenerators
    ALLOW_UNDERSCORES
    MODULE_LIBRARIES
        VLNData
        VLNVarExtractors
        VLNDataGenerators
    BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
