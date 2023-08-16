cet_make_library(
  SOURCE
  DQMFlashAlg.cxx
  DQMHitAlg.cxx
  LIBRARIES
  PUBLIC
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  DQMChannelNoiseModule art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
)

cet_build_plugin(
  DQMFlashModule art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
)

cet_build_plugin(
  DQMHitModule art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
)

cet_build_plugin(
  DQMTPCBurstModule art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
)

cet_build_plugin(
  DQMTrackModule art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
)

cet_build_plugin(
  GoodRunSelectionAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
  ubobj::Trigger
)

cet_build_plugin(
  HVBurstFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubcore::DQMTools
)

install_headers()
install_fhicl()
install_source()
