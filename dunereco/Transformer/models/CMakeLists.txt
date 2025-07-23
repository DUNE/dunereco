
cet_add_compiler_flags(CXX -Wno-pedantic)

art_make_library(
    LIBRARY_NAME TransformerModels
    SOURCE
        zoo/TransformerEnergyModel.cxx
        pt_model/ModelConfig.cxx
        pt_model/PTModel.cxx
        pt_model/TorchHandler.cxx
    LIBRARIES
        torch
        torch_cpu
        c10
	    Boost::boost
)

install_headers(SUBDIRS pt_model zoo)
install_source(SUBDIRS pt_model zoo)

