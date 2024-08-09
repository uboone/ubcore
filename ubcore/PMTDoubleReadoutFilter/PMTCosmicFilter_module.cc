////////////////////////////////////////////////////////////////////////
// Class:       PMTCosmicFilter
// Plugin Type: producer (art v3_01_02)
// File:        PMTCosmicFilter_module.cc
//
// Generated at Thu Aug  8 14:51:03 2024 by Herbert Greenlee using cetskelgen
// from cetlib version v3_05_01.
//
// Description:
//
// This producer module reads an input OpDetWaveform containing a mixture of
// beam and cosmic waveforms, and produces a new output OpDetWaveform containing
// only the cosmic waveforms from the input data product.
//
// Beam and cosmic waveforms are distringuished by their length.  Beam waveforms
// normally have a length of 1501, and cosmic waveforms normally have a length
// of 41.  This module uses a length cut set by fcl parameter.  There is no
// selection based on timing.
//
// Fcl parameters:
//
// inputLabel: Module label and instance name of input waveform
//             (<module>:<instance>).
// outputInstance: Instance name of output waveform
// maxLength: Maximum length of cosmic waveforms.
// channelIncrement: Increment for channel numbers in output waveform.
//                
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include <memory>

class PMTCosmicFilter;


class PMTCosmicFilter : public art::EDProducer {
public:
  explicit PMTCosmicFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTCosmicFilter(PMTCosmicFilter const&) = delete;
  PMTCosmicFilter(PMTCosmicFilter&&) = delete;
  PMTCosmicFilter& operator=(PMTCosmicFilter const&) = delete;
  PMTCosmicFilter& operator=(PMTCosmicFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  // Fcl parameters.

  art::InputTag fInputLabel;    // Module label + instance name of input waveform.
  std::string fOutputInstance;  // Instance name of output waveform.
  unsigned int fMaxLength;      // Maximum length cut.
  int fChannelIncrement;        // Channel increment.
};


PMTCosmicFilter::PMTCosmicFilter(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fInputLabel(p.get<art::InputTag>("inputLabel")),
  fOutputInstance(p.get<std::string>("outputInstance")),
  fMaxLength(p.get<unsigned int>("maxLength")),
  fChannelIncrement(p.get<int>("channelIncrement"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  produces<std::vector<raw::OpDetWaveform> >(fOutputInstance);
  consumes<std::vector<raw::OpDetWaveform> >(fInputLabel);
}

void PMTCosmicFilter::produce(art::Event& e)
{
  // Implementation of required member function here.

  // Get the input waveform.

  art::Handle<std::vector<raw::OpDetWaveform> > wfh;
  e.getByLabel(fInputLabel, wfh);

  // We consider it a fatal error if the input waveform is not found.

  if(!wfh.isValid()) {
    throw art::Exception(art::errors::ProductNotFound)
      << "Could not find input waveform with label " << fInputLabel;
  }

  // Construct the new data product that will hold the selected waveforms.

  std::unique_ptr<std::vector<raw::OpDetWaveform> > new_wfvec(new std::vector<raw::OpDetWaveform>);
  new_wfvec->reserve(wfh->size());

  // Loop over OpDetWaveforms.

  for(auto const& wf : *wfh) {

    // Maybe insert waveform into new data product.

    if(wf.size() <= fMaxLength) {
      new_wfvec->push_back(wf);

      // Increment channel number.

      new_wfvec->back().SetChannelNumber(wf.ChannelNumber() + fChannelIncrement);
    }
  }

  // Insert data product into event.

  e.put(std::move(new_wfvec), fOutputInstance);
}

DEFINE_ART_MODULE(PMTCosmicFilter)
