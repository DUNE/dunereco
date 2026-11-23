cet_build_plugin(FlashMatchBDT art::EDProducer
  SOURCE
  FlashMatchBDT_module.cc
  FlashMatchBDTAlg.cxx
  BDTModel.cxx
  LIBRARIES PRIVATE
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  canvas::canvas
  cetlib::cetlib
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_fhicl(LIST bdt_augmented.bin bdt_baseline.bin SUBDIRNAME FlashMatch)
install_source()
install_scripts(LIST scripts/export_bin_bdt.py)
