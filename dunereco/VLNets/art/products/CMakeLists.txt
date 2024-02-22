art_make(
    LIBRARY_NAME VLNProducts
    LIB_LIBRARIES
        VLNData
    BASENAME_ONLY
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()

