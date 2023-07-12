cet_build_plugin(
  PMTDoubleReadoutFilter art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RawData
)

install_headers()
install_fhicl()
install_source()
