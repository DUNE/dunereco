include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  add_subdirectory(algorithms)
  add_subdirectory(fcls)
  add_subdirectory(modules)
  add_subdirectory(tf)
endif ()

add_subdirectory(products)

install_headers()
install_fhicl()
install_source()

