if ( DEFINED ENV{LIBTORCH_DIR})
   add_subdirectory(art)
   add_subdirectory(func)
endif()
# nnbar/ is a standalone tree (python and ROOT scripts of the n-nbar study); it is not built or installed.
