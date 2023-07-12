cet_build_plugin(
  SampledInputInfo art::EDProducer
  LIBRARIES
  PRIVATE
  larcoreobj::SummaryData
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_source()
install_fhicl()
