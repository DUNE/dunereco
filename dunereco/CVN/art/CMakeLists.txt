include_directories( $ENV{CAFFE_INC} )
include_directories( $ENV{CAFFE_SRC} )
include_directories( $ENV{GSL_INC}/gsl )
include_directories( $ENV{HDF5_INC} )
include_directories( $ENV{PROTOBUF_INC} )
include_directories( $ENV{LEVELDB_INC} )
include_directories( $ENV{GFLAGS_INC} )
include_directories( $ENV{GLOG_INC} )

cet_find_library( CAFFELIB NAMES caffe PATHS ENV CAFFE_LIB NO_DEFAULT_PATH )
cet_find_library( GSLLIB NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )
cet_find_library( HDF5LIB NAMES hdf5 PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
cet_find_library( PROTOBUFLIB NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
cet_find_library( LEVELDBLIB NAMES leveldb PATHS ENV LEVELDB_LIB NO_DEFAULT_PATH )
cet_find_library( GFLAGSLIB NAMES gflags PATHS ENV GFLAGS_LIB NO_DEFAULT_PATH )
cet_find_library( GLOGLIB NAMES glog PATHS ENV GLOG_LIB NO_DEFAULT_PATH )

  art_make(BASENAME_ONLY
           LIBRARY_NAME      CVNArt
           LIB_LIBRARIES     ${ART_FRAMEWORK_CORE}
                             ${ART_FRAMEWORK_PRINCIPAL}
                             ${ART_FRAMEWORK_SERVICES_REGISTRY}
                             ${ART_FRAMEWORK_SERVICES_BASIC}
                             ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                             ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                             ${ART_PERSISTENCY_COMMON}
                             canvas
                             art_Persistency_Provenance
                             canvas 
                             ${ART_UTILITIES}
                             canvas
                             ${MF_MESSAGELOGGER}
                             ${MF_UTILITIES}
                             ${FHICLCPP}
                             cetlib cetlib_except
                             ${ROOT_BASIC_LIB_LIST}
                             CVNFunc
                             lardataobj_RecoBase
                             lardata_Utilities
                             nusimdata_SimulationBase
                             larsim_MCCheater_BackTrackerService_service           
                             ${CAFFELIB}
                             ${GSLLIB}                  
                             ${HDF5LIB}
                             ${PROTOBUFLIB}
                             ${LEVELDBLIB}
                             ${GFLAGSLIB}
                             ${GLOGLIB}                             

           MODULE_LIBRARIES  CVNFunc
                             CVNArt
)


install_headers()
install_fhicl()
install_source()
