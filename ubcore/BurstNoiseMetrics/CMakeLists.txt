cet_build_plugin(
  BurstNoiseMetrics art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::AssociationUtil
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  BurstNoiseMetricsFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardata::AssociationUtil
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

install_headers()
install_fhicl()
install_source()

