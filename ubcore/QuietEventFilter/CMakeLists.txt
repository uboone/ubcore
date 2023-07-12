cet_build_plugin(
  QuietEventAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  QuietEventFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
)

install_headers()
install_source()
install_fhicl()
