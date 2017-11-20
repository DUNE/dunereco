#
# rocco@physics.umn.edu 
#

#include SoftRelTools/arch_spec_root.mk
#include SoftRelTools/arch_spec_novadaq.mk

#override CPPFLAGS := $(CPPFLAGS)

BINS = $(basename $(wildcard *.cc))  # This finds all the .cc files, drops .cc 
SIMPLEBINS = $(BINS)


BINLIBS := $(BINLIBS) 

#override CPPFLAGS := $(CPPFLAGS) -I$(FANNHOME)/include -I$(BOOST_INC) -I$(NUTOOLS_INC) -I$(NOVADAQ_INC) -I$(CAFFE_INC) -I$(CAFFE_SRC) -I$(GSL_INC) -I$(GSL_INC)/gsl  -I$(PROTOBUF_INC) -I$(GLOG_INC) -I$(LEVELDB_INC) -I$(LMDB_INC) -I$(GFLAGS_INC) -I$(HDF5_INC)
override CPPFLAGS := $(CPPFLAGS) -I$(BOOST_INC) -I$(NUTOOLS_INC) -I$(CAFFE_INC) -I$(CAFFE_SRC) -I$(GSL_INC) -I$(GSL_INC)/gsl  -I$(PROTOBUF_INC) -I$(GLOG_INC) -I$(LEVELDB_INC) -I$(LMDB_INC) -I$(GFLAGS_INC) -I$(HDF5_INC)

########################################################################

#include SoftRelTools/standard.mk
include standard.mk
#include SoftRelTools/arch_spec_art.mk
#include SoftRelTools/arch_spec_wda.mk

override BINLIBS += $(LOADLIBES)  -L$(BOOST_LIB) -lboost_program_options -L$(FHICLCPP_LIB)/ -L$(CETLIB_LIB)/ -lfhiclcpp -lcetlib -L$(ART_LIB) -lart_Utilities -lCintex -lCVNFunc -lCalibratorFunc -lCalibrator -lCalibrator_service -lCMap_service -lUtilitiesFunc -lLiveGeometry_service -lBackTracker_service -L$(CAFFE_LIB) -lcaffe -lCVN -L$(HDF5_LIB) -lhdf5_cpp -lhdf5


