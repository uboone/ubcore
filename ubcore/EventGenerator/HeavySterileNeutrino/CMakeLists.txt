art_make(
      MODULE_LIBRARIES
            HSNGenHelpers
            HSNGenDataObjects
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom::RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support ROOT::Core
            art_root_io::TFileService_service
            art::Framework_Services_Optional_RandomNumberGenerator_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            messagefacility::MF_MessageLogger
            
            fhiclcpp::fhiclcpp
            CLHEP::CLHEP
            SQLite::SQLite3
            ROOT::EG
            ROOT::Tree
            )

add_subdirectory(Helpers)
add_subdirectory(DataObjects)

install_headers()
install_fhicl()
install_source()

