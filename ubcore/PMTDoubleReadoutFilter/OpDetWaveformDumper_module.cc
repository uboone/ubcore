////////////////////////////////////////////////////////////////////////
// Class:       OpDetWaveformDumper
// Plugin Type: analyzer (art v3_01_02)
// File:        OpDetWaveformDumper_module.cc
//
// Generated at Thu Aug  8 11:37:55 2024 by Herbert Greenlee using cetskelgen
// from cetlib version v3_05_01.
//
// Fcl parameters:
//
// triggerLabel: Module label of Trigger data product.
// waveformLabel: Module label of OpDetWaveform data product.
//                This can a simple module label ("<module>") or a
//                label and instance name separated by a colon 
//                ("<module>:<instance>"). or it can be an empty string.
//                If nonempty, only matching waveforms are dumped.
// processName: Process name (if empty, dump all processes).
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

class OpDetWaveformDumper;


class OpDetWaveformDumper : public art::EDAnalyzer {
public:
  explicit OpDetWaveformDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpDetWaveformDumper(OpDetWaveformDumper const&) = delete;
  OpDetWaveformDumper(OpDetWaveformDumper&&) = delete;
  OpDetWaveformDumper& operator=(OpDetWaveformDumper const&) = delete;
  OpDetWaveformDumper& operator=(OpDetWaveformDumper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

  // Fcl parameters.

  std::string fTriggerLabel;
  std::string fWaveformLabel;
  std::string fInstanceName;
  std::string fProcessName;
};


OpDetWaveformDumper::OpDetWaveformDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fTriggerLabel(p.get<std::string>("triggerLabel")),
  fWaveformLabel(p.get<std::string>("waveformLabel")),
  fProcessName(p.get<std::string>("processName"))
{
  // Parse module label parameter into module + instance.

  size_t n = fWaveformLabel.find(':');
  if(n < std::string::npos) {
    fInstanceName = fWaveformLabel.substr(n+1);
    fWaveformLabel = fWaveformLabel.substr(0,n);
  }
  std::cout << "OpDetWaveformDumper module configured with the following parameters:" << std::endl;
  std::cout << "Trigger label = " << fTriggerLabel << std::endl;
  std::cout << "Waveform label = " << fWaveformLabel << std::endl;
  std::cout << "Waveform instance = " << fInstanceName << std::endl;
  std::cout << "Process name = " << fProcessName << std::endl;

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<raw::OpDetWaveform>>(fWaveformLabel);
  consumes<std::vector<raw::Trigger>>(fTriggerLabel);
}

void OpDetWaveformDumper::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  std::cout << "Dumping OpDetWaveform information." << std::endl;

  // Get raw::Trigger data product.

  art::Handle<std::vector<raw::Trigger> > trigh;
  e.getByLabel(fTriggerLabel, trigh);
  if(!trigh.isValid()) {
    std::cout << "No Trigger data product with module label " << fTriggerLabel << std::endl;
    return;
  }

  // Don't know what to do if trigger data product length is different than one.

  if(trigh->size() != 1) {
    std::cout << "Trigger data product bad length = " << trigh->size() << std::endl;
    return;
  }

  // Get trigger time.

  double trigger_time = trigh->front().TriggerTime();
  //std::cout << "Trigger time = " << trigger_time << std::endl;

  // Get all raw::OpDetWaveform data products.

  std::vector<art::Handle<std::vector<raw::OpDetWaveform> > > wfhs;
  e.getManyByType(wfhs);
  std::cout << wfhs.size() << " waveform types." << std::endl;

  // Filter handles to retain only matching labels and instances.

  std::vector<art::Handle<std::vector<raw::OpDetWaveform> > > filtered_wfhs;
  for(const auto& wfh : wfhs) {
    const art::Provenance* prov = wfh.provenance();
    if((fWaveformLabel.empty() || fWaveformLabel == prov->moduleLabel()) &&
       (fInstanceName.empty() || fInstanceName == prov->productInstanceName()) &&
       (fProcessName.empty() || fProcessName == prov->processName()))
      filtered_wfhs.push_back(wfh);
  }
  std::cout << filtered_wfhs.size() << " filtered waveform types." << std::endl;

  // Loop over handles.

  for(const auto& wfh : filtered_wfhs) {

    // Count the number of channels.

    std::set<unsigned int> chset;
    for(auto const& wf : *wfh) {
      unsigned int ch = wf.ChannelNumber();
      if(chset.count(ch) == 0)
        chset.insert(ch);
    }
    const art::Provenance* prov = wfh.provenance();
    std::cout << "\nModule label = " << prov->moduleLabel() << std::endl;
    std::cout << "Instance name = " << prov->productInstanceName() << std::endl;
    std::cout << "Process name = " << prov->processName() << std::endl;
    std::cout << wfh->size() << " waveforms" << std::endl;
    std::cout << chset.size() << " channels" << std::endl;

    // Loop over OpDetWaveforms.

    for(auto const& wf : *wfh) {
      unsigned int ch = wf.ChannelNumber();
      size_t len = wf.size();
      double wftime = wf.TimeStamp() - trigger_time;

      short adc_max = 0;
      double adc_tot = 0.;
      int adc_num = 0;
      double adc_avg = 0.;
      for(auto adc : wf) {
        ++adc_num;
        if(adc > adc_max)
          adc_max = adc;
        adc_tot += adc;
      }
      if(adc_num > 0)
        adc_avg = adc_tot / adc_num;
      std::cout << "\nChennel " << ch << std::endl;
      std::cout << "Length = " << len << std::endl;
      std::cout << "Start time = " << wftime << std::endl;
      std::cout << "Average ADC = " << adc_avg << std::endl;
      std::cout << "Maximum ADC = " << adc_max << std::endl;
    }
  }
}

DEFINE_ART_MODULE(OpDetWaveformDumper)
