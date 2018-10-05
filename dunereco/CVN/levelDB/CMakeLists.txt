include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{CAFFE_INC} )
include_directories( $ENV{GSL_INC}/gsl )
include_directories( $ENV{HDF5_INC} )
include_directories( $ENV{PROTOBUF_INC} )
include_directories( $ENV{LEVELDB_INC} )
include_directories( $ENV{LMDB_INC} )
include_directories( $ENV{BOOST_INC} )
include_directories( $ENV{GFLAGS_INC} )
include_directories( $ENV{GLOG_INC} )

cet_find_library( CAFFELIB NAMES caffe-d caffe PATHS ENV CAFFE_LIB NO_DEFAULT_PATH )
cet_find_library( HDF5_CPPLIB NAMES hdf5_cpp PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
cet_find_library( HDF5LIB NAMES hdf5 PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
cet_find_library( PROTOBUFLIB NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
cet_find_library( LEVELDBLIB NAMES leveldb PATHS ENV LEVELDB_LIB NO_DEFAULT_PATH )
cet_find_library( LMDBLIB NAMES lmdb PATHS ENV LMDB_LIB NO_DEFAULT_PATH )
cet_find_library( BOOSTLIB NAMES boost_program_options PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( GFLAGSLIB NAMES gflags PATHS ENV GFLAGS_LIB NO_DEFAULT_PATH )
cet_find_library( GLOGLIB NAMES glog PATHS ENV GLOG_LIB NO_DEFAULT_PATH )

message("HDF5: ${HDF5_CPPLIB}")
message("Boost manual lib      : ${BOOSTLIB}")

art_make_exec( cvnCreateDB         
               SOURCE    cvnCreateDB.cc
               LIBRARIES ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_BASIC}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                         ${ART_PERSISTENCY_COMMON}
                         canvas
                         art_Persistency_Provenance
                         ${ART_UTILITIES}
                         ${BOOSTLIB}
                         ${CAFFELIB}
                         ${HDF5_CPPLIB} ${HDF5LIB}
                         ${PROTOBUFLIB}
                         ${LEVELDBLIB}
                         ${LMDBLIB}
                         ${GFLAGSLIB}
                         ${GLOGLIB}
                         ${FHICLCPP}
                         ${CETLIB}
                         ${ROOT_BASIC_LIB_LIST}
                         CVNFunc
               )

art_make_exec( cvnCreateZlibImages       
               SOURCE    cvnCreateZlibImages.cc
               LIBRARIES ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_BASIC}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                         ${ART_PERSISTENCY_COMMON}
                         canvas
                         art_Persistency_Provenance
                         z
                         ${ART_UTILITIES}
                         ${BOOSTLIB}
                         ${GLOGLIB}
                         ${FHICLCPP}
                         ${CETLIB}
                         ${ROOT_BASIC_LIB_LIST}
                         CVNFunc
               )


install_source()
install_fhicl()
