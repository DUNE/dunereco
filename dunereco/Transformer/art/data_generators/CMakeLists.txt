art_make(
    LIBRARY_NAME TransformerDataGenerators
    ALLOW_UNDERSCORES
    MODULE_LIBRARIES
        TransformerData
        TransformerVarExtractors
        TransformerDataGenerators
    BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
