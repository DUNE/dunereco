art_make_library(
    LIBRARY_NAME TransformerData
    SOURCE
        exporters/CSVExporter.cxx
        structs/VarDict.h
        structs/TransformerEnergy.h
)

install_headers(SUBDIRS exporters structs)
install_source(SUBDIRS exporters structs)

