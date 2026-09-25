if ( DEFINED ENV{LIBTORCH_DIR})
   add_subdirectory(art)
   add_subdirectory(func)
endif()
add_subdirectory(nnbar)   # scripts only, no libtorch needed
