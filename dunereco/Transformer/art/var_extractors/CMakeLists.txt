art_make_library(
    LIBRARY_NAME TransformerVarExtractors
    SOURCE
        DefaultInputVarExtractor.cxx
        EventAddrVarExtractor.cxx
        EventMCVarExtractor.cxx
        EventRecoEVarExtractor.cxx
        EventRecoVarExtractor.cxx
        FiducialCutVarExtractor.cxx
        PFParticleVarExtractor.cxx
        TransformerEnergyVarExtractor.cxx
        VarExtractorBase.cxx
        utils.cxx
    LIBRARIES
        dunereco_AnaUtils
        larreco::Calorimetry
)

install_headers()
install_source()
install_fhicl()

