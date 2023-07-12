cet_make_library(
  SOURCE
  UBTriggerAlgo.cxx
  LIBRARIES
  PUBLIC
  lardata::DetectorClocksService
  larcore::ServiceUtil
  ROOT::Core
)

cet_build_plugin(
  TriggerTiming art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::TriggerSim
  art_root_io::TFileService_service
)

cet_build_plugin(
  UBTriggerSim art::EDProducer
  LIBRARIES
  PRIVATE
  ubcore::TriggerSim
  art_root_io::TFileService_service
)

install_headers()
install_fhicl()
install_source()
