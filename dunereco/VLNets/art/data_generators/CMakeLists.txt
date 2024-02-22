art_make(
    LIBRARY_NAME VLNDataGenerators
    ALLOW_UNDERSCORES
    MODULE_LIBRARIES
        VLNData
        VLNVarExtractors
        VLNDataGenerators
    BASENAME_ONLY
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
