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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

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
    std::vector<double> fTriggerCalib; ///< user-defined calibration trigger in G4 ns (per-event)
    std::vector<double> fTriggerPC;    ///< user-defined PC trigger in G4 ns (per-event)
    std::vector<double> fTriggerExt;   ///< user-defined Ext trigger in G4 ns (per-event)
    std::vector<double> fTriggerBNB;   ///< user-defined BNB trigger in G4 ns (per-event)
    std::vector<double> fTriggerNuMI;  ///< user-defined NuMI trigger in G4 ns (per-event)

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
  UBTriggerSim::UBTriggerSim(fhicl::ParameterSet const& pset)
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
    fTriggerCalib = pset.get< std::vector<double> > ("CalibTrigger",std::vector<double>());
    fTriggerExt   = pset.get< std::vector<double> > ("ExtTrigger",std::vector<double>());
    fTriggerPC    = pset.get< std::vector<double> > ("PCTrigger",std::vector<double>());
    fTriggerBNB   = pset.get< std::vector<double> > ("BNBTrigger",std::vector<double>());
    fTriggerNuMI  = pset.get< std::vector<double> > ("NuMITrigger",std::vector<double>());

    fBNBFireTime  = pset.get<double>("BNBFireTime");
    fNuMIFireTime = pset.get<double>("NuMIFireTime");

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
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    std::unique_ptr< std::vector<raw::Trigger>   >  triggers(new std::vector<raw::Trigger>);
    fAlg.ClearInputTriggers();
    auto clock = ts->OpticalClock();

    // Register triggers
    for(auto const t : fTriggerCalib) {
      auto const elec_time = ts->G4ToElecTime(t);
      clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
      fAlg.AddTriggerCalib(clock);
    }

    for(auto const t : fTriggerExt) {
      auto const elec_time = ts->G4ToElecTime(t);
      clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
      fAlg.AddTriggerExt(clock);
    }

    for(auto const t : fTriggerPC) {
      auto const elec_time = ts->G4ToElecTime(t);
      clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
      fAlg.AddTriggerPC(clock);
    }

    for(auto const t : fTriggerBNB) {
      auto const elec_time = ts->G4ToElecTime(t);
      clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
      fAlg.AddTriggerBNB(clock);
    }

    for(auto const t : fTriggerNuMI) {
      auto const elec_time = ts->G4ToElecTime(t);
      clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
      fAlg.AddTriggerNuMI(clock);
    }

    // Register input BNB/NuMI beam trigger
    for(auto const &name : fBeamModName) {
      art::Handle<std::vector<sim::BeamGateInfo> > beamArray;
      evt.getByLabel(name,beamArray);
      if(!(beamArray.isValid())) continue;
      for(size_t i=0; i<beamArray->size(); ++i) {
	art::Ptr<sim::BeamGateInfo> beam_ptr (beamArray,i);
	if(beam_ptr->BeamType() == sim::kBNB) {
	  //auto const elec_time = ts->G4ToElecTime(beam_ptr->Start());
	  auto const elec_time = ts->G4ToElecTime(fBNBFireTime);
	  clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
	  fAlg.AddTriggerBNB(clock);
	}
	else if(beam_ptr->BeamType() == sim::kNuMI) {
	  //auto const elec_time = ts->G4ToElecTime(beam_ptr->Start());
	  auto const elec_time = ts->G4ToElecTime(fNuMIFireTime);
	  clock.SetTime(clock.Sample(elec_time),clock.Frame(elec_time));
	  fAlg.AddTriggerNuMI(clock);
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

	for(size_t i=0; i<pmtArray->size(); ++i) {

	  art::Ptr<optdata::PMTTrigger> pmt_ptr (pmtArray,i);
	  clock.SetTime( pmt_ptr->TimeSlice(), pmt_ptr->Frame() );
	  switch(pmt_ptr->Category()) {
	  case optdata::kCosmicPMTTrigger:
	    fAlg.AddTriggerPMTCosmic(clock);
	    break;
	  case optdata::kBeamPMTTrigger:
	    fAlg.AddTriggerPMTBeam(clock);
	    break;
	  default:
	    throw UBTrigException(Form("Unexpected OpticalCategory (%d) found from PMTTrigger",pmt_ptr->Category()));
	    return;
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
