////////////////////////////////////////////////////////////////////////
//
//  UBTriggerAlgo source
//
////////////////////////////////////////////////////////////////////////

#ifndef UBTRIGGERALGO_CXX
#define UBTRIGGERALGO_CXX

#include <iostream>
#include <sstream>
#include <TString.h>
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "ubobj/Trigger/UBTriggerTypes.h"
#include "UBTriggerAlgo.h"

namespace trigger{

  //##############################################################
  UBTriggerAlgo::UBTriggerAlgo(detinfo::DetectorClocksData const& clockData) :
    _tpc_clock{clockData.TPCClock()},
    _pmt_clock{clockData.OpticalClock()},
    _trig_clock{clockData.TriggerClock()},
    _mask(9,0),
				   _prescale(9,false)
  //##############################################################
  {
    // LArLight version constructor
    SetDeadtime(4);
    SetBNBParams(102,256,0,0);
    SetNuMIParams(102,256,0,0);
  }

  //######################################################
  void UBTriggerAlgo::Report(const std::string &msg) const
  //######################################################
  {
    std::cout << msg.c_str() << std::endl;
  }

  //#########################################################
  void UBTriggerAlgo::SetBNBParams(unsigned short width,
				   unsigned short delay,
				   unsigned short cosmic_min,
				   unsigned short cosmic_max) 
  //#########################################################
  {
    if(cosmic_max < cosmic_min)
      throw UBTrigException(Form("BNB Cosmic allow window max (%d) is smaller than min (%d)!",
				 cosmic_min,
				 cosmic_max)
			    );
    _bnb_gate_width       = width;
    _bnb_delay            = delay;
    _bnb_cosmic_allow_min = cosmic_min;
    _bnb_cosmic_allow_max = cosmic_max;
  }

  //##########################################################
  void UBTriggerAlgo::SetNuMIParams(unsigned short width,
				    unsigned short delay,
				    unsigned short cosmic_min,
				    unsigned short cosmic_max) 
  //##########################################################
  {
    if( cosmic_max < cosmic_min)
      throw UBTrigException(Form("NuMI Cosmic allow window max (%d) is smaller than min (%d)!",
				 cosmic_min,
				 cosmic_max)
			    );
    _numi_gate_width       = width;
    _numi_delay            = delay;
    _numi_cosmic_allow_min = cosmic_min;
    _numi_cosmic_allow_max = cosmic_max;
  }


  //#############################################################
  void UBTriggerAlgo::SetMask(unsigned char index, uint32_t mask)
  //#############################################################
  {
    if(index>9) throw UBTrigException("Index >9 cannot be set!");
    _mask.at(index)=mask;
  }

  //#################################################################
  void UBTriggerAlgo::SetPrescale(unsigned char index, bool prescale)
  //#################################################################
  {
    if(index>9) throw UBTrigException("Index >9 cannot be set!");
    _prescale.at(index)=prescale;
  }

  //##################################################################
  void UBTriggerAlgo::SetMask(const std::vector<uint32_t> &mask) 
  //##################################################################
  {
    if(_mask.size()!=9) {
      throw UBTrigException(Form("Length of masks (=%zu) is invalid! Initializing to right length...",
				 mask.size()));
      _mask.resize(9,0);
    }else
      _mask = mask;

  }

  //##################################################################
  void UBTriggerAlgo::SetPrescale(const std::vector<bool> &prescale) 
  //##################################################################
  {
    if(_prescale.size()!=9) {
      throw UBTrigException(Form("Length of prescales (=%zu) is invalid! Initializing to right length...",
				 prescale.size()));
      _prescale.resize(9,0);
    }else
      _prescale = prescale;

  }

  //######################################
  void UBTriggerAlgo::ClearInputTriggers(detinfo::DetectorClocksData const& clockData)
  //######################################
  {
    _candidates.clear();

    for( auto const t : _bnb_timings ) AddTriggerBNB(clockData, t);

    for( auto const t : _numi_timings ) AddTriggerNuMI(clockData, t);

    for( auto const t : _calib_timings ) AddTriggerCalib(clockData, t);

    for( auto const t : _ext_timings ) AddTriggerExt(clockData, t);

    for( auto const t : _pc_timings ) AddTriggerPC(clockData, t);

  }

  //############################################################################
  detinfo::ElecClock UBTriggerAlgo::BNBStartTime(const detinfo::DetectorClocksData& clockData,
                                                 const detinfo::ElecClock& time) const
  //############################################################################
  {
    return clockData.OpticalClock().WithTime(_pmt_clock.Time(time.Time())).AdvanceTimeBy(_bnb_delay);
  }

  //#############################################################################
  detinfo::ElecClock UBTriggerAlgo::NuMIStartTime(const detinfo::DetectorClocksData& clockData,
                                                  const detinfo::ElecClock& time) const
  //#############################################################################
  {
    return clockData.OpticalClock().WithTime(_pmt_clock.Time(time.Time())).AdvanceTimeBy(_numi_delay);
  }

  //##############################################################
  void UBTriggerAlgo::AddTriggerCalib(const detinfo::DetectorClocksData& clockData,
                                      const detinfo::ElecClock &time)
  //##############################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
    // Calibration triggers
    AddTrigger(clockData, raw::Trigger(0,
			    _pmt_clock.Time(),
			    _pmt_clock.Time(),
			    ( (0x1) << trigger::kTriggerCalib )) 
	       );
    _pmt_clock = _pmt_clock.WithTime(0);
  }

  //############################################################
  void UBTriggerAlgo::AddTriggerExt(const detinfo::DetectorClocksData& clockData,
                                    const detinfo::ElecClock& time)
  //############################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
    // EXT triggers
    AddTrigger(clockData, raw::Trigger(0,
			    _pmt_clock.Time(),
			    _pmt_clock.Time(),
			    ( (0x1) << trigger::kTriggerEXT ))
	       );
    _pmt_clock = _pmt_clock.WithTime(0);
  }

  //###########################################################
  void UBTriggerAlgo::AddTriggerPC(const detinfo::DetectorClocksData& clockData,
                                   const detinfo::ElecClock& time)
  //###########################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
    // PC triggers
    AddTrigger(clockData, raw::Trigger(0,
			    _pmt_clock.Time(),
			    _pmt_clock.Time(),
			    ( (0x1) << trigger::kTriggerPC ))
	       );
    _pmt_clock = _pmt_clock.WithTime(0);
  }

  //#######################################
  void UBTriggerAlgo::ReportConfig() const
  //#######################################
  {
    std::ostringstream msg;

    msg
      << std::endl
      << " UBTriggerAlgo Configuration:              " << std::endl
      << "---------------------------------------------" << std::endl;

    msg << " Debug Mode ... " << (_debug_mode ? "enabled!" : "disabled!") << std::endl
	<< std::endl;

    for(size_t i=0; i<_mask.size(); ++i)

      msg << Form("  Mask %zu     : %d", i, _mask[i]) << std::endl;
    
    msg << std::endl;

    for(size_t i=0; i<_prescale.size(); ++i)

      msg << Form("  Prescale %zu : %d", i, int(_prescale[i])) << std::endl;

    msg << std::endl;

    msg
      << Form(" Trigger Deadtime    : %d frames",_deadtime) << std::endl
      << std::endl
      << Form("  NuMI Trigger Delay : %d samples",_numi_delay) << std::endl
      << Form("  NuMI Cosmic Start  : %d samples",_numi_cosmic_allow_min) << std::endl
      << Form("  NuMI Cosmic End    : %d samples",_numi_cosmic_allow_max) << std::endl
      << std::endl
      << Form("  BNB Trigger Delay  : %d samples",_bnb_delay) << std::endl
      << Form("  BNB Cosmic Start   : %d samples",_bnb_cosmic_allow_min) << std::endl
      << Form("  BNB Cosmic End     : %d samples",_bnb_cosmic_allow_max) << std::endl
      << "---------------------------------------------" << std::endl
      << std::endl;

    //mf::LogInfo(__FUNCTION__) << msg.str();
    Report(msg.str());
  }

  //#######################################################################################
  const raw::Trigger UBTriggerAlgo::CombineTriggers(const detinfo::DetectorClocksData& clockData,
                                                    const raw::Trigger &trigger1,
						    const raw::Trigger &trigger2)
  //#######################################################################################
  {

    auto trig1_time = trigger1.TriggerTime();
    auto trig2_time = trigger2.TriggerTime();

    auto trig1_number = trigger1.TriggerNumber();
    auto trig2_number = trigger2.TriggerNumber();

    if(_debug_mode)
      
      Report(Form("  Attempting to combine two triggers: (%d,%d) and (%d,%d)",
		  _trig_clock.Sample(trig1_time),
		  _trig_clock.Frame(trig1_time),
		  _trig_clock.Sample(trig2_time),
		  _trig_clock.Frame(trig2_time))
	     );

    if( trig1_number != trig2_number ) {
      throw UBTrigException("Cannot combine triggers with different trigger counters!");
      return raw::Trigger();
    }

    if( _trig_clock.Frame(trig1_time) != _trig_clock.Frame(trig2_time) ) {
      throw UBTrigException("Cannot combine triggers in different frames!");
      return raw::Trigger();
    }

    if( _trig_clock.Sample(trig1_time) != _trig_clock.Sample(trig2_time) ) {
      throw UBTrigException("Cannot combine triggers in different trigger clock sample number!");
      return raw::Trigger();
    }

    auto const res_number = trig1_number;
    auto const res_frame  = _pmt_clock.Frame(trig1_time);
    auto const trig1_sample = _pmt_clock.Sample(trig1_time);
    auto const trig2_sample = _pmt_clock.Sample(trig2_time);
    auto const res_sample = (trig1_sample < trig2_sample ? trig1_sample : trig2_sample);

    // Construct beam gate region. Note that this is relevant only if either of
    // two candidates is beam. Else do nothing.
    unsigned int beam_frame  = 0;
    unsigned int beam_sample = 0;

    // Case1: both are beam ... not supported for now
    if( (trigger1.Triggered(trigger::kTriggerBNB) || trigger1.Triggered(trigger::kTriggerNuMI)) &&
	(trigger2.Triggered(trigger::kTriggerBNB) || trigger2.Triggered(trigger::kTriggerNuMI)) ) {
      throw UBTrigException("Combining two beam gates not supported for now!");
      return raw::Trigger();
    }
    // Case2: only trigger 1 is beam
    else if( (trigger1.Triggered(trigger::kTriggerBNB) || trigger1.Triggered(trigger::kTriggerNuMI)) ) {
      
      beam_frame  = _pmt_clock.Frame(trigger1.BeamGateTime());
      beam_sample = _pmt_clock.Sample(trigger1.BeamGateTime());

    }
    // Case3: only trigger 2 is beam
    else if( (trigger2.Triggered(trigger::kTriggerBNB) || trigger2.Triggered(trigger::kTriggerNuMI)) ) {

      beam_frame  = _pmt_clock.Frame(trigger2.BeamGateTime());
      beam_sample = _pmt_clock.Sample(trigger2.BeamGateTime());

    }
    // Case4: neither is beam ... use trigger timing
    else{

      beam_frame  = res_frame;
      beam_sample = res_sample;

    }

    if(_debug_mode) Report(Form("    Combined beam timing @ (%d,%d)", beam_sample, beam_frame));

    // Combine trigger bits
    uint32_t res_bits = 0x0;
    for(unsigned char i=0; i < 4*(sizeof(res_bits)); ++i) {
      
      if( trigger1.Triggered(i) || trigger2.Triggered(i) )
	{
	  res_bits += (uint32_t)( (0x1) << i );
	  if(_debug_mode) Report(Form("    Combined bit %d ... now %d",i,res_bits));
	}
    }
    auto const& clock = clockData.OpticalClock();
    auto trig_time = clock.WithTime(clock.Time(res_sample,  res_frame));
    auto beam_time = clock.WithTime(clock.Time(beam_sample, beam_frame));
    return raw::Trigger(res_number,
			trig_time.Time(),
			beam_time.Time(),
			res_bits);
    
  }

  //##################################################################
  void UBTriggerAlgo::AddTrigger(detinfo::DetectorClocksData const& clockData,
                                 const raw::Trigger &new_trigger)
  //##################################################################
  {
    auto frame  = _trig_clock.Frame(new_trigger.TriggerTime());
    auto sample = _trig_clock.Sample(new_trigger.TriggerTime());

    if(_debug_mode)

      Report(Form("Requested to add a new trigger @ (%d, %d)",sample,frame));

    auto frame_iter = _candidates.find(frame);

    if( frame_iter == _candidates.end() ) {

      _candidates[frame].insert(std::pair<unsigned int,raw::Trigger>(sample,new_trigger));
      
    }else{

      auto sample_iter = (*frame_iter).second.find(sample);

      if( sample_iter == (*frame_iter).second.end() )

	(*frame_iter).second.insert(std::pair<unsigned int,raw::Trigger>(sample,new_trigger));
      
      else {

        raw::Trigger combined_trigger = CombineTriggers(clockData, new_trigger, (*sample_iter).second);

	_candidates[frame][sample]=combined_trigger;

      }
    }
  }

  //##################################################################
  void UBTriggerAlgo::AddTriggerPMTCosmic(const detinfo::DetectorClocksData& clockData,
                                          const detinfo::ElecClock& time)
  //##################################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
    // Trigger bits
    unsigned int trig_bits=0x0;
    trig_bits += ( (0x1) << kPMTTriggerCosmic );
    trig_bits += ( (0x1) << kPMTTrigger );
    
    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
				_pmt_clock.Time(),
				trig_bits);
    
    // Add this trigger candidate
    AddTrigger(clockData, trig_candidate);
    _pmt_clock = _pmt_clock.WithTime(0);
  }

  //################################################################
  void UBTriggerAlgo::AddTriggerPMTBeam(const detinfo::DetectorClocksData& clockData,
                                        const detinfo::ElecClock& time)
  //################################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
   // Trigger bits
    unsigned int trig_bits=0x0;
    trig_bits += ( (0x1) << kPMTTriggerBeam );
    trig_bits += ( (0x1) << kPMTTrigger );
    
    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
				_pmt_clock.Time(),
				trig_bits);
    
    // Add this trigger candidate
    AddTrigger(clockData, trig_candidate);
    _pmt_clock = _pmt_clock.WithTime(0);
  }

  //############################################################
  void UBTriggerAlgo::AddTriggerBNB(const detinfo::DetectorClocksData& clockData,
                                    const detinfo::ElecClock& time)
  //############################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
    uint32_t trig_bits = ( (0x1) << trigger::kTriggerBNB );
 
    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
                                BNBStartTime(clockData, _pmt_clock).Time(),
				trig_bits);
    
    // Add this trigger candidate
    AddTrigger(clockData, trig_candidate);
    _pmt_clock = _pmt_clock.WithTime(0);
  }

  //#############################################################
  void UBTriggerAlgo::AddTriggerNuMI(const detinfo::DetectorClocksData& clockData,
                                     const detinfo::ElecClock& time)
  //#############################################################
  {
    _pmt_clock = _pmt_clock.WithTime(_pmt_clock.Time(time.Time()));
    uint32_t trig_bits = ( (0x1) << trigger::kTriggerNuMI );

    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
                                NuMIStartTime(clockData, _pmt_clock).Time(),
				trig_bits);

    // Add this trigger candidate
    AddTrigger(clockData, trig_candidate);
    _pmt_clock = _pmt_clock.WithTime(0);
  }
  
  //###############################################
  void UBTriggerAlgo::ShowCandidateTriggers() const
  //###############################################
  {

    auto frame_iter = _candidates.begin();
    std::ostringstream msg;

    msg
      << "######################################################"
      << "Trigger Inputs"
      << "######################################################" << std::endl
      << std::endl
      << Form("  %-18s : %-18s : %-18s : %-18s : %-18s : %-18s : %-32s",
	      "Trig Sample(16MHz)",
	      "Trig Frame (16MHz)",
	      "Trig Sample(64MHz)",
	      "Trig Frame (64MHz)",
	      "Beam Sample(64Mhz)",
	      "Beam Frame (64Mhz)",
	      "Bits") << std::endl
      << std::endl;

    while(frame_iter!=_candidates.end()) {
      
      auto sample_iter = (*frame_iter).second.begin();

      while(sample_iter != (*frame_iter).second.end()) {

	unsigned int decision_sample = (*sample_iter).first;
	unsigned int decision_frame  = (*frame_iter).first;

	unsigned int trigger_sample  = _pmt_clock.Sample(((*sample_iter).second.TriggerTime()));
	unsigned int trigger_frame   = _pmt_clock.Frame(((*sample_iter).second.TriggerTime()));

	unsigned int beam_sample     = _pmt_clock.Sample(((*sample_iter).second.BeamGateTime()));
	unsigned int beam_frame      = _pmt_clock.Sample(((*sample_iter).second.BeamGateTime()));
	
	msg
	  << Form("  %-18d : %-18d : %-18d : %-18d : %-18d : %-18d : ",
		  decision_sample,
		  decision_frame,
		  trigger_sample,
		  trigger_frame,
		  beam_sample,
		  beam_frame);

	for(unsigned char k=0; k<32; ++k)

	  msg << ((*sample_iter).second.Triggered(k) ? 1 : 0);

	msg << std::endl;

	sample_iter++;
      }
      frame_iter++;
    }
    msg
      << std::endl
      << "######################################################"
      << "##############"
      << "######################################################" << std::endl
      << std::endl;

    Report(msg.str());
  }

  //#########################################################################
  void UBTriggerAlgo::ProcessTrigger(std::vector<raw::Trigger> &triggers)
  //#########################################################################
  {

    triggers.clear();

    if(_debug_mode)
      
      ShowCandidateTriggers();

    // Just make sure trigger clock is only used for sample/frame getter (i.e. no time used)
    _trig_clock = _trig_clock.WithTime(::detinfo::kTIME_MAX);
    detinfo::ElecClock trig_time = _trig_clock;

    time_window_t bnb_active  (_trig_clock,_trig_clock);
    time_window_t numi_active (_trig_clock,_trig_clock);
    time_window_t bnb_gate    (_trig_clock,_trig_clock);
    time_window_t numi_gate   (_trig_clock,_trig_clock);
    time_window_t deadtime    (_trig_clock,_trig_clock);

    _trig_clock = _trig_clock.WithTime(0);

    auto const mask0 = _mask.at(0);
    auto const mask1 = _mask.at(1);
    auto const mask8 = _mask.at(8);
    auto const scale0 = _prescale.at(0);
    auto const scale1 = _prescale.at(1);
    auto const scale8 = _prescale.at(8);

    // Loop over candidate frames
    auto frame_iter = _candidates.begin();
    while( frame_iter != _candidates.end() ) {

      // Make sure sample exists if frame is found
      auto sample_iter = (*frame_iter).second.begin();
      
      if( sample_iter == (*frame_iter).second.end() ) {
	
	throw UBTrigException("Logic error: found candidate frame but no associated sample & trigger!");
	triggers.clear();
	return;
      }    

      // Loop over samples
      while( sample_iter != (*frame_iter).second.end() ) {

	if(_debug_mode) Report(Form("\n  Processing candidate @ %d, %d",
				    (*sample_iter).first,
				    (*frame_iter).first));

        trig_time = trig_time.WithTime(_trig_clock.Time((*sample_iter).first,(*frame_iter).first));

	// If in deadtime, continue
	if( InWindow(trig_time, deadtime) ) {

	  if(_debug_mode) Report(Form("    Skipping as deadtime is (%d,%d) => (%d,%d)",
				      deadtime.first.Sample(),
				      deadtime.first.Frame(),
				      deadtime.second.Sample(),
				      deadtime.second.Frame())
				 );
	  sample_iter++;
	  continue;
	}else if(_debug_mode){
	  
	  Report(Form("    Not in deadtime: (%d,%d) => (%d,%d)",
		      deadtime.first.Sample(),
		      deadtime.first.Frame(),
		      deadtime.second.Sample(),
		      deadtime.second.Frame())
		 );
	}

	auto const pmt0   = (*sample_iter).second.Triggered(trigger::kPMTTriggerBeam);
	auto const pmt1   = (*sample_iter).second.Triggered(trigger::kPMTTriggerCosmic);
	auto const numi   = (*sample_iter).second.Triggered(trigger::kTriggerNuMI);
	auto const bnb    = (*sample_iter).second.Triggered(trigger::kTriggerBNB);
	auto const calib  = (*sample_iter).second.Triggered(trigger::kTriggerCalib);
	auto const ext    = (*sample_iter).second.Triggered(trigger::kTriggerEXT);
	auto const pc     = (*sample_iter).second.Triggered(trigger::kTriggerPC);

	// Evaluate "active" = this trigger is within bnb/numi cosmic allow window
	auto const active = ( InWindow(trig_time, bnb_active) || InWindow(trig_time, numi_active) );

	// If BNB or NuMI, open new beam trigger gate
	if(bnb) {
          auto& [min, max] = bnb_gate;
          auto new_min_time = min.Time((*sample_iter).first, (*frame_iter).first);
          min = min.WithTime(new_min_time);
          max = max.WithTime((bnb_gate.first.Time() + _trig_clock.Time(_bnb_gate_width)));
	}
	if(numi) {
          auto& [min, max] = numi_gate;
          auto new_min_time = min.Time((*sample_iter).first, (*frame_iter).first);
          min = min.WithTime(new_min_time);
          max = max.WithTime(numi_gate.first.Time() + _trig_clock.Time(_numi_gate_width));
	}
	
	auto const bnb_gate_active  = InWindow(trig_time, bnb_gate);  // true = this trigger is within bnb gate
	auto const numi_gate_active = InWindow(trig_time, numi_gate); // true = this trigger is within numi gate

	if(_debug_mode)
	  {
	    std::ostringstream msg;
	    msg 
	      << Form("    PMT Beam?   %s", (pmt0  ? "yes" : "no")) << std::endl
	      << Form("    PMT Cosmic? %s", (pmt1  ? "yes" : "no")) << std::endl
	      << Form("    NuMI Gate?  %s", (numi  ? "yes" : "no")) << std::endl
	      << Form("    BNB Gate?   %s", (bnb   ? "yes" : "no")) << std::endl
	      << Form("    Calib?      %s", (calib ? "yes" : "no")) << std::endl
	      << Form("    External?   %s", (ext   ? "yes" : "no")) << std::endl
	      << Form("    PC?         %s", (pc    ? "yes" : "no")) << std::endl
	      << std::endl
	      << Form("    NuMI   Window? %s", (numi_gate_active ? "yes" : "no")) << std::endl
	      << Form("    BNB    Window? %s", (bnb_gate_active  ? "yes" : "no")) << std::endl
	      << Form("    Cosmic Window? %s", (active           ? "yes" : "no")) << std::endl;
	    Report(msg.str());
	  }



	// Evaludate trigger condition p0 ... this is for PMT Cosmic trigger
	// Equation follows Nevis FPGA code including similar notations
	bool p0 = ( (numi_gate_active && pmt0 && ((0x1) & (mask0 >> 0))           ) ||
		    (bnb_gate_active  && pmt0 && ((0x1) & (mask0 >> 1))           ) ||
		    (active           && pmt0 && ((0x1) & (mask0 >> 2)) && !scale0) ||
		    (                    pmt0 && ((0x1) & (mask0 >> 3)) && !scale0) );

	// Evaluate trigger condition p1 ... this is for PMT Beam trigger
	// Equation follows Nevis FPGA code including similar notations
	bool p1 = ( (numi_gate_active && pmt1 && ((0x1) & (mask1 >> 0))           ) || 
		    (bnb_gate_active  && pmt1 && ((0x1) & (mask1 >> 1))           ) ||
		    (active           && pmt1 && ((0x1) & (mask1 >> 2)) && !scale1) ||
		    (                    pmt1 && ((0x1) & (mask1 >> 3)) && !scale1) );

	// Evaluate trigger condition p8 ... this is for Beam trigger input
	// Equation follows Nevis FPGA code including similar notations
	bool p8 = ( (active && ext    && ((0x1) & (mask8 >> 0)) && !scale8) ||
		    (          ext    && ((0x1) & (mask8 >> 1)) && !scale8) ||
		    (active && pc     && ((0x1) & (mask8 >> 2))           ) ||
		    (          pc     && ((0x1) & (mask8 >> 3))           ) ||
		    (bnb_gate_active  && ((0x1) & (mask8 >> 4))           ) ||
		    (numi_gate_active && ((0x1) & (mask8 >> 5))           ) ||
		    (           calib && ((0x1) & (mask8 >> 6))           ) );

	if(_debug_mode) {
	  
	  std::ostringstream msg;

	  msg << Form("    p0 condition ... %s", (p0 ? "satisfied" : "not met")) << std::endl
	      << Form("    p1 condition ... %s", (p1 ? "satisfied" : "not met")) << std::endl
	      << Form("    p8 condition ... %s", (p8 ? "satisfied" : "not met")) << std::endl;
	  Report(msg.str());
	}


	// Ask OR condition of p0, p1, and p8
	if( p0 || p1 || p8 ) {
	  // New trigger found. 
          auto& [min, max] = deadtime;
          min = min.WithTime(trig_time.Time());
          max = max.WithTime(trig_time.Time()).AdvanceTimeBy((int)(_deadtime * _trig_clock.FrameTicks() - 1));
	  if(_debug_mode) {
	    std::cout<<Form("    Dead time set: (%d,%d) => (%d,%d) with trigger @ (%d,%d)",
                            min.Sample(),
                            min.Frame(),
                            max.Sample(),
                            max.Frame(),
			    trig_time.Sample(),
			    trig_time.Frame()
			    )
		     << std::endl;
	  }

	  // Assign active window for bnb/numi
	  if(bnb) {
            auto& [min, max] = bnb_active;
            min = min.WithTime(trig_time.Time()).AdvanceTimeBy((int)(_bnb_cosmic_allow_min));
            max = max.WithTime(trig_time.Time()).AdvanceTimeBy((int)(_bnb_cosmic_allow_max));
	  }
	  if(numi) {
            auto& [min, max] = numi_active;
            min = min.WithTime(trig_time.Time()).AdvanceTimeBy((int)(_numi_cosmic_allow_min));
            max = max.WithTime(trig_time.Time()).AdvanceTimeBy((int)(_numi_cosmic_allow_max));
	  }
	  
	  // Store trigger object
	  triggers.push_back( raw::Trigger(_trigger_counter,
					   (*sample_iter).second.TriggerTime(),
					   (*sample_iter).second.TriggerTime(),//(*sample_iter).second.BeamGateTime(),
					   (*sample_iter).second.TriggerBits()) 
			      );
	  _trigger_counter++;
	} // end if trigger found
	sample_iter++;
      } // end looping over samples
      frame_iter++;
    } // end looping over frames
    
  }

}

#endif
