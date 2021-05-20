art_make_library(
    LIBRARY_NAME VLNData
    SOURCE
        exporters/CSVExporter.cxx
        structs/VarDict.h
        structs/VLNEnergy.h
)

install_headers(SUBDIRS exporters structs)
install_source(SUBDIRS exporters structs)

