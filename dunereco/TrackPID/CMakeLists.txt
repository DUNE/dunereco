include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

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

