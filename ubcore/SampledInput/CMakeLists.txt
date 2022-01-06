art_make(
          LIBRARY_NAME ubcore_SampledInput
          # Needed for standalone c++ helper files
          LIB_LIBRARIES

          # Needed for the main module
          MODULE_LIBRARIES
                    larcoreobj_SummaryData
                    nusimdata::SimulationBase
                    art::Framework_Core
                    art::Framework_Principal
                    art::Framework_Services_Registry
                    #
		    art_root_io::tfile_support ROOT::Core
		    art_root_io::TFileService_service
	            art_root_io::TFileService_service
                    art::Framework_Services_Registry
                    art::Persistency_Common
                    art::Persistency_Provenance
                    art::Utilities
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    cetlib_except::cetlib_except
                    ROOT::Tree
          )

install_headers()
install_source()
install_fhicl()

