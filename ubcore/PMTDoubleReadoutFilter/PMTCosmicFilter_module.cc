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
// This producer module reads input OpDetWaveforms containing a mixture of
// beam and cosmic waveforms, and produces a new output OpDetWaveform containing
// only the cosmic waveforms from the input data products.
//
// Beam and cosmic waveforms are distringuished by their length.  Beam waveforms
// normally have a length of 1501, and cosmic waveforms normally have a length
// of 41.  This module uses a length cut set by fcl parameter.  There is no
// selection based on timing.
//
// Fcl parameters:
//
// inputLabel: Input waveform module label (just the label without the instance name).
// inputInstance: Input waveform instance name.
// outputInstance: Output waveform instance name.
// maxLength: Maximum length of cosmic waveforms.
//
// Usage notes:
//
// 1.  The default mode in which this module is intended to be used will have
//     the following fcl parameters.
//
//     inputLabel = "pmtreadout"   (read OpDetWaveforms produced by swizzler).
//     inputInstance = "OpdetBeamHighGain"  (where HG cosmic waveforms are found in run 1a data).
//     ouputInstance = "OpdetCosmicHighGain"  (normal instance for HG cosmic waveforms).
//
// 2.  This module will attempt to read cosmic waveforms from two input tags:
//     a) <inputLabel>:<inputInstance>
//     b) <inputLabel>:<outputInstance>
//
// Input tag (a) is where cosmic waveforms would be found in run 1a data.
// Input tag (b) retains existing cosmic waveforms that have the correct instance already.
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

  std::string fInputLabel;      // Module label of input waveform.
  std::string fInputInstance;   // Instance name of input waveform.
  std::string fOutputInstance;  // Instance name of output waveform.
  unsigned int fMaxLength;      // Maximum length cut.
};


PMTCosmicFilter::PMTCosmicFilter(fhicl::ParameterSet const& p) :
  EDProducer{p},
  fInputLabel(p.get<std::string>("inputLabel")),
  fInputInstance(p.get<std::string>("inputInstance")),
  fOutputInstance(p.get<std::string>("outputInstance")),
  fMaxLength(p.get<unsigned int>("maxLength"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  produces<std::vector<raw::OpDetWaveform> >(fOutputInstance);
  consumes<std::vector<raw::OpDetWaveform> >(art::InputTag(fInputLabel, fInputInstance));
  consumes<std::vector<raw::OpDetWaveform> >(art::InputTag(fInputLabel, fOutputInstance));
}

void PMTCosmicFilter::produce(art::Event& e)
{
  // Implementation of required member function here.

  // Get the input waveforms.
  // We try to read the following two waveforms.
  // 1.  <inputLabel>:<inputInstance>
  // 2.  <inputLabel>:<outputInstance>

  std::vector<art::Handle<std::vector<raw::OpDetWaveform> > > wfhvec(2);
  e.getByLabel(fInputLabel, fInputInstance, wfhvec[0]);
  e.getByLabel(fInputLabel, fOutputInstance, wfhvec[1]);

  // Construct the new data product that will hold the selected waveforms.

  std::unique_ptr<std::vector<raw::OpDetWaveform> > new_wfvec(new std::vector<raw::OpDetWaveform>);

  // Loop over handles.

  for(auto const& wfh : wfhvec) {

    // It is not an error if handle is not valid.

    if(wfh.isValid()) {

      if(new_wfvec->size() == 0)
        new_wfvec->reserve(wfh->size());

      // Loop over OpDetWaveforms from this handle.

      for(auto const& wf : *wfh) {

        // Maybe insert waveform into new data product.

        if(wf.size() <= fMaxLength) {
          new_wfvec->push_back(wf);

          // Adjust channel number.

          new_wfvec->back().SetChannelNumber(wf.ChannelNumber() % 100 + 200);
        }
      }
    }
  }

  // Insert data product into event.

  e.put(std::move(new_wfvec), fOutputInstance);
}

DEFINE_ART_MODULE(PMTCosmicFilter)
