cet_build_plugin(
  DLSignalSample art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::MCBase
  nusimdata::SimulationBase
)

install_headers()
install_source()
install_fhicl()
