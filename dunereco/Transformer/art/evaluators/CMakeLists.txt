art_make(
    LIBRARY_NAME TransformerEvaluators
    MODULE_LIBRARIES
        TransformerData
        TransformerModels
        TransformerVarExtractors
        TransformerDataGenerators
    BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

