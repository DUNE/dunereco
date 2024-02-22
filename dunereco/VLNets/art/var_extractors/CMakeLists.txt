art_make_library(
    LIBRARY_NAME VLNVarExtractors
    SOURCE
        DefaultInputVarExtractor.cxx
        EventAddrVarExtractor.cxx
        EventMCVarExtractor.cxx
        EventRecoEVarExtractor.cxx
        EventRecoVarExtractor.cxx
        FiducialCutVarExtractor.cxx
        PFParticleVarExtractor.cxx
        VLNEnergyVarExtractor.cxx
        VarExtractorBase.cxx
        utils.cxx
    LIBRARIES
        dunereco_AnaUtils
        larreco::Calorimetry
)

install_headers()
install_source()
install_fhicl()

