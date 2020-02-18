////////////////////////////////////////////////////////////////////////
// Class:       PMTDoubleReadoutFilter
// Plugin Type: producer (art v3_04_00)
// File:        PMTDoubleReadoutFilter_module.cc
//
// Generated at Wed Jan 29 21:07:05 2020 by Pawel Guzowski using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

#include <memory>

namespace pdrf {
  class PMTDoubleReadoutFilter;



  class PMTDoubleReadoutFilter : public art::EDProducer {
    public:
      explicit PMTDoubleReadoutFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      PMTDoubleReadoutFilter(PMTDoubleReadoutFilter const&) = delete;
      PMTDoubleReadoutFilter(PMTDoubleReadoutFilter&&) = delete;
      PMTDoubleReadoutFilter& operator=(PMTDoubleReadoutFilter const&) = delete;
      PMTDoubleReadoutFilter& operator=(PMTDoubleReadoutFilter&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      art::InputTag fTriggerLabel, fWaveformLabel;
      unsigned int fNumOfChannelsExpected;
      double fTimestampDiffTolerance;
      // Declare member data here.

  };

}

pdrf::PMTDoubleReadoutFilter::PMTDoubleReadoutFilter(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTriggerLabel(p.get<art::InputTag>("triggerLabel")),
  fWaveformLabel(p.get<art::InputTag>("waveformLabel")),
  fNumOfChannelsExpected(p.get<unsigned int>("numberOfChannelsToExpect")),
  fTimestampDiffTolerance(p.get<double>("timestampDiffTolerance"))// ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  produces<std::vector<raw::OpDetWaveform>>(fWaveformLabel.instance());

  consumes<std::vector<raw::OpDetWaveform>>(fWaveformLabel);
  consumes<std::vector<raw::Trigger>>(fTriggerLabel);
}

void pdrf::PMTDoubleReadoutFilter::produce(art::Event& e)
{
  // Implementation of required member function here.
  double TriggerTime = -1.;
  art::Handle<std::vector<raw::Trigger>> trigHandle;
  if(e.getByLabel(fTriggerLabel, trigHandle)) {
    if(trigHandle->size() != 1) {
      throw art::Exception(art::errors::DataCorruption) << "Could not find unique trigger data in event.\n";
    }
    TriggerTime = trigHandle->at(0).TriggerTime();
  }
  else throw art::Exception(art::errors::ProductNotFound) << "Could not find trigger data with label '"<<fTriggerLabel<<"' in event.\n";

  art::Handle<std::vector<raw::OpDetWaveform>> wfHandle;
  if(e.getByLabel(fWaveformLabel, wfHandle)) {
    std::unique_ptr< std::vector<raw::OpDetWaveform>   > new_waveforms(new std::vector<raw::OpDetWaveform>);
    new_waveforms->reserve(fNumOfChannelsExpected);

    // 'seen' is used to ensure there is only one single waveform for each channel, and all fNumOfChannelsExpected channels have one.
    std::vector<int> seen(fNumOfChannelsExpected,0);
    for(auto const& wf : *wfHandle) {
      // we want to select only the waveforms with their timestamp at the trigger time
      if(std::abs(TriggerTime - wf.TimeStamp()) < fTimestampDiffTolerance) {
        if(wf.ChannelNumber() < fNumOfChannelsExpected) {
          seen[wf.ChannelNumber()] ++;
        }
        else throw art::Exception(art::errors::DataCorruption) << "PMT channel seen with bad channel ID "<<wf.ChannelNumber()<<".\n";
        new_waveforms->push_back(wf);
      }
    }
    if(std::count(seen.begin(), seen.end(), 1) != fNumOfChannelsExpected) {
      throw art::Exception(art::errors::DataCorruption) << "Could not find "<<fNumOfChannelsExpected<<" PMT waveforms at the trigger time.\n";
    }

    e.put(std::move(new_waveforms),fWaveformLabel.instance());
  }
  else throw art::Exception(art::errors::ProductNotFound) << "Could not find waveform data with label '"<<fWaveformLabel<<"' in event.\n";
}

DEFINE_ART_MODULE(pdrf::PMTDoubleReadoutFilter)
