
art_make( BASENAME_ONLY
	  LIBRARY_NAME  uboonecode_TriggerSim
	  LIB_LIBRARIES lardataobj_RawData
			lardata_Utilities
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Persistency_Common canvas::canvas
			art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

			cetlib::cetlib
                        cetlib_except::cetlib_except

	  MODULE_LIBRARIES
			nug4::ParticleNavigation lardataobj_Simulation
			lardataobj_RawData
			lardataobj_OpticalDetectorData
			lardata_Utilities
			uboonecode_TriggerSim
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except

          SERVICE_LIBRARIES
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except

        )

install_headers()
install_fhicl()
install_source()
