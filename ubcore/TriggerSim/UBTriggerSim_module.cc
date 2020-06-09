/**
 * \file UBTriggerSim_module.cc
 *
 * \ingroup UBTriggerSim
 * 
 * \brief raw::Trigger producer module for microboone
 *
 * @author kazuhiro
 */

/** \addtogroup UBTriggerSim

@{*/

#ifndef UBTriggerSim_H
#define UBTriggerSim_H

// LArSoft includes
//#include "Geometry/Geometry.h"

// framework
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

/// LArSoft
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "UBTrigException.h"
#include "UBTriggerAlgo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/OpticalDetectorData/PMTTrigger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

/// nutools
#include "lardataobj/Simulation/BeamGateInfo.h"

namespace trigger {
  /**
     \class UBTriggerSim
     MicroBooNE's producer module for raw::Trigger data product.
     Current implementation takes (1) trigger module configuration (UBTriggerAlgo configuration) and
     (2) G4-event-wise Calibration, External, PC, BNB or NuMI trigger input with G4-time as time-stamp.
     Note (2) is treated in the same way as pulse input into trigger module in real electronics.
     It is possible to also set (2) per run or per sub-run, but such scheme is probably not useful for
     the current implementation of event generation/simulation scheme. Talk to the author if such option
     becomes useful/possible and needs to be implemented.
  */ 
  class UBTriggerSim : public art::EDProducer{
  public:

    /// Default ctor
    UBTriggerSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~UBTriggerSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&); 

  private:

    /// Algorithm
    UBTriggerAlgo fAlg;

    /// BNB delay
    double fBNBFireTime;
    /// NuMI delay
    double fNuMIFireTime;

    /// Module labels for event generator(s) that produced sim::BeamGateInfo data product
    std::vector<std::string> fBeamModName;

    /// Module label for OpticalFEM
    std::string fOpticalFEMMod;

    //-- ElecClock for user-defined trigger times --//
    std::vector<detinfo::ElecClock> fTriggerCalib; ///< user-defined calibration trigger (per-event)
    std::vector<detinfo::ElecClock> fTriggerPC;    ///< user-defined PC trigger (per-event)
    std::vector<detinfo::ElecClock> fTriggerExt;   ///< user-defined Ext trigger (per-event)
    std::vector<detinfo::ElecClock> fTriggerBNB;   ///< user-defined BNB trigger (per-event)
    std::vector<detinfo::ElecClock> fTriggerNuMI;  ///< user-defined NuMI trigger (per-event)

  };

} 

#endif//  UBTriggerSim_H

// UBTriggerSim.cc

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
namespace trigger {
  DEFINE_ART_MODULE(UBTriggerSim)
}

namespace trigger {

  //#########################################################
  UBTriggerSim::UBTriggerSim(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    fAlg{art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob()}
  //#########################################################
  {
    /// Get beam event generator label
    fBeamModName = pset.get< std::vector<std::string> > ("BeamModName", std::vector<std::string>());

    /// Get OpticalFEM module label
    fOpticalFEMMod = pset.get< std::string > ("OpticalFEMMod","");

    // Get trigger module configuration parameters
    fAlg.SetDebugMode  ( pset.get< bool                  > ("DebugMode"      ) );
    fAlg.SetMask       ( pset.get< std::vector<uint32_t> > ("Mask"           ) );
    fAlg.SetPrescale   ( pset.get< std::vector<bool>     > ("Prescale"       ) );
    fAlg.SetDeadtime   ( pset.get< unsigned short        > ("DeadTime"       ) );

    fAlg.SetBNBParams  ( pset.get< unsigned short        > ("BNBGateWidth"   ),
			 pset.get< unsigned short        > ("BNBGateDelay"   ),
			 pset.get< unsigned short        > ("BNBCosmicStart" ),
			 pset.get< unsigned short        > ("BNBCosmicEnd"   ) );

    fAlg.SetNuMIParams ( pset.get< unsigned short        > ("NuMIGateWidth"   ),
			 pset.get< unsigned short        > ("NuMIGateDelay"   ),
			 pset.get< unsigned short        > ("NuMICosmicStart" ),
			 pset.get< unsigned short        > ("NuMICosmicEnd"   ) );

    // Get user-defined trigger timings
    std::vector<double> trig_calib ( pset.get< std::vector<double> > ("CalibTrigger",
								      std::vector<double>()) 
				     );
    std::vector<double> trig_ext   ( pset.get< std::vector<double> > ("ExtTrigger",
								      std::vector<double>()) 
				     );
    std::vector<double> trig_pc    ( pset.get< std::vector<double> > ("PCTrigger",
								      std::vector<double>()) 
				     );
    std::vector<double> trig_bnb   ( pset.get< std::vector<double> > ("BNBTrigger",
								      std::vector<double>()) 
				     );
    std::vector<double> trig_numi  ( pset.get< std::vector<double> > ("NuMITrigger",
								      std::vector<double>()) 
				     );

    fBNBFireTime  = pset.get<double>("BNBFireTime");
    fNuMIFireTime = pset.get<double>("NuMIFireTime");

    // Store user-defined trigger timings to the attributes
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto clock = clockData.OpticalClock();

    fTriggerCalib.clear();
    fTriggerExt.clear();
    fTriggerPC.clear();
    fTriggerBNB.clear();
    fTriggerNuMI.clear();

    auto clock_with_time = [&clockData, &clock](auto const g4_time) {
                             auto const elec_time = clockData.G4ToElecTime(g4_time);
                             auto const clock_time = clock.Time(clock.Sample(elec_time), clock.Frame(elec_time));
                             return clock.WithTime(clock_time);
                           };

    fTriggerCalib.reserve(trig_calib.size());
    for(auto const t : trig_calib) {
      fTriggerCalib.push_back(clock_with_time(t));
    }

    fTriggerExt.reserve(trig_ext.size());
    for(auto const t : trig_ext) {
      fTriggerExt.push_back(clock_with_time(t));
    }

    fTriggerPC.reserve(trig_pc.size());
    for(auto const t : trig_pc) {
      fTriggerPC.push_back(clock_with_time(t));
    }

    fTriggerBNB.reserve(trig_bnb.size());
    for(auto const t : trig_bnb) {
      fTriggerBNB.push_back(clock_with_time(t));
    }

    fTriggerNuMI.reserve(trig_numi.size());
    for(auto const t : trig_numi) {
      fTriggerNuMI.push_back(clock_with_time(t));
    }

    // Produces raw::Trigger data product
    produces< std::vector<raw::Trigger> >();
  }

  //###########################
  UBTriggerSim::~UBTriggerSim()
  //###########################
  {}
   
  //#########################################
  void UBTriggerSim::produce(art::Event& evt) 
  //#########################################
  {
    // Initialize
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    std::unique_ptr< std::vector<raw::Trigger>   >  triggers(new std::vector<raw::Trigger>);
    fAlg.ClearInputTriggers(clockData);
    auto clock = clockData.OpticalClock();

    // Register user-defined triggers
    for( auto const &t : fTriggerCalib ) fAlg.AddTriggerCalib (clockData, t);
    for( auto const &t : fTriggerExt   ) fAlg.AddTriggerExt   (clockData, t);
    for( auto const &t : fTriggerPC    ) fAlg.AddTriggerPC    (clockData, t);
    for( auto const &t : fTriggerBNB   ) fAlg.AddTriggerBNB   (clockData, t);
    for( auto const &t : fTriggerNuMI  ) fAlg.AddTriggerNuMI  (clockData, t);

    auto clock_with_time = [&clockData, &clock](auto const g4_time) {
                             auto const elec_time = clockData.G4ToElecTime(g4_time);
                             auto const clock_time = clock.Time(clock.Sample(elec_time), clock.Frame(elec_time));
                             return clock.WithTime(clock_time);
                           };

    // Register input BNB/NuMI beam trigger
    for(auto const &name : fBeamModName) {
      art::Handle<std::vector<sim::BeamGateInfo> > beamArray;
      evt.getByLabel(name,beamArray);
      if(!(beamArray.isValid())) continue;
      for(size_t i=0; i<beamArray->size(); ++i) {
	art::Ptr<sim::BeamGateInfo> beam_ptr (beamArray,i);
	if(beam_ptr->BeamType() == sim::kBNB) {
          fAlg.AddTriggerBNB(clockData, clock_with_time(fBNBFireTime));
	}
	else if(beam_ptr->BeamType() == sim::kNuMI) {
          fAlg.AddTriggerNuMI(clockData, clock_with_time(fNuMIFireTime));
	}
	else
	  throw UBTrigException(Form("Beam type %d not supported!",beam_ptr->BeamType()));
      }
    }

    // Register input PMTTrigger
    if(!(fOpticalFEMMod.empty())) {

      art::Handle<std::vector<optdata::PMTTrigger> > pmtArray;
      evt.getByLabel(fOpticalFEMMod,pmtArray);
      if(pmtArray.isValid()) {

        for(auto const& pmt_trigger : *pmtArray) {
          auto const clock_time = clock.Time(pmt_trigger.TimeSlice(), pmt_trigger.Frame());
          auto const new_clock = clock.WithTime(clock_time);
          switch(pmt_trigger.Category()) {
	  case optdata::kCosmicPMTTrigger:
            fAlg.AddTriggerPMTCosmic(clockData, new_clock);
	    break;
	  case optdata::kBeamPMTTrigger:
            fAlg.AddTriggerPMTBeam(clockData, new_clock);
	    break;
	  default:
            throw UBTrigException(Form("Unexpected OpticalCategory (%d) found from PMTTrigger",pmt_trigger.Category()));
	  }
	}
      }
    }

    // Run trigger simulation
    fAlg.ProcessTrigger(*triggers);
    
    // Store
    //if(triggers->size()) // now unnecessary
      evt.put(std::move(triggers));
  }
} 
/** @} */ // end of doxygen group 
