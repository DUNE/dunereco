add_subdirectory(products)
add_subdirectory(var_extractors)
add_subdirectory(data_generators)

if ((DEFINED ENV{TENSORFLOW_DIR}) AND (DEFINED ENV{TRITON_DIR}))
    add_subdirectory(evaluators)
endif ()

