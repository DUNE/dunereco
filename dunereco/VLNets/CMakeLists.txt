add_subdirectory(art)
add_subdirectory(data)

if ((DEFINED ENV{TENSORFLOW_DIR}) AND (DEFINED ENV{TRTIS_CLIENTS_DIR}))
    add_subdirectory(models)
endif ()

