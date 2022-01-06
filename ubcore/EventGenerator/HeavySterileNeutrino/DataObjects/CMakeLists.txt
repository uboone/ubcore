art_make(
      LIBRARY_NAME HSNGenDataObjects
      LIB_LIBRARIES
            nusimdata::SimulationBase
            messagefacility::MF_MessageLogger
            
            fhiclcpp::fhiclcpp
            cetlib::cetlib
            CLHEP::CLHEP
      MODULE_LIBRARIES
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
            canvas::canvas
            messagefacility::MF_MessageLogger
            
            fhiclcpp::fhiclcpp
            cetlib::cetlib
            cetlib_except::cetlib_except
            CLHEP::CLHEP
            SQLite::SQLite3
            ROOT::EG
            ifdh::ifdh
            )

install_headers()
install_fhicl()
install_source()
