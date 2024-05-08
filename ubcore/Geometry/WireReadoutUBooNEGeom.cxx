////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapUBooNEAlg.cxx
/// \brief MicroBooNE optical detector channel mapping
///   
/// Terminology:
///  OpDet: the physical PMT or light guide paddle
///  HardwareChannel: The gain copy number. 0=Low A, 1 = Low B, 3=High A, 3 = High B  
///  OpChannel: The readout channel number
///
/// \version $Id:  $
/// \author  taritree@mit.edu
////////////////////////////////////////////////////////////////////////

#include "ubcore/Geometry/WireReadoutUBooNEGeom.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/WireReadoutSorterStandard.h"

namespace geo {

  //----------------------------------------------------------------------------
  WireReadoutUBooNEGeom::WireReadoutUBooNEGeom( geo::GeometryCore const* geom,
                                            fhicl::ParameterSet const& pvals)
    : WireReadoutStandardGeom( pvals, geom , std::make_unique<WireReadoutSorterStandard>() )
  {
    // parameter set will come from UBooNEGomeotryHelper service
    LoadOpticalMapData( pvals );
  }

  //----------------------------------------------------------------------------
  WireReadoutUBooNEGeom::~WireReadoutUBooNEGeom()
  {
  }

  //----------------------------------------------------------------------------
  // OPTICAL CHANNELS
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  unsigned int WireReadoutUBooNEGeom::NOpChannels(unsigned int NOpDets) const
  {
    return fNReadoutChannels;
  }

  unsigned int WireReadoutUBooNEGeom::MaxOpChannel(unsigned int NOpDets) const
  {
    return fMaxOpChannel;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutUBooNEGeom::NOpHardwareChannels(unsigned int opDet) const
  {
    auto it = fPMT2channels.find( opDet );
    if ( it!=fPMT2channels.end() )
      return (*it).second.size();
    else
      throw std::runtime_error( "WireReadoutUBooNEGeom::NOpHardwareChannels : Invalid opdet value" );
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutUBooNEGeom::OpChannel(unsigned int pmtID, unsigned int copynum) const
  {
    // converts OpDet and Gain Channel
    unsigned int uniqueChannel = fPMT2channels.at(pmtID).at(copynum);
    return uniqueChannel;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutUBooNEGeom::OpDetFromOpChannel(unsigned int opChannel) const
  {
    unsigned int pmtID = fChannel2pmt.at( opChannel );
    return pmtID;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutUBooNEGeom::HardwareChannelFromOpChannel(unsigned int opChannel) const
  {
    auto it = fPMT2channels.find( OpDetFromOpChannel(opChannel)  );
    if ( it!=fPMT2channels.end() ) {
      auto readoutChList = it->second;
      for ( unsigned int hwch=0; hwch<readoutChList.size(); hwch++ )
	if ( opChannel==readoutChList.at(hwch) )
	  return hwch;
    }

    throw std::runtime_error( "Could not find copy index of readout channel" );
  }

  //----------------------------------------------------------------------------
  bool WireReadoutUBooNEGeom::IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const {
    auto it=fChannel2pmt.find( opChannel );
    if ( it!=fChannel2pmt.end() ) {
      return true;
    }
    return false;
  }
  
  //----------------------------------------------------------------------------
  void WireReadoutUBooNEGeom::LoadOpticalMapData( fhicl::ParameterSet const& pset ) {
    fNOpDets = pset.get< unsigned int >("numberOfDetectors");
    fNReadoutChannels = 0;    
    fMaxOpChannel = 0;
    // ----------------------------------------------------------------------
    // map between opdet ID and Readout Channel Number
    for (unsigned int iop=0; iop<fNOpDets; iop++) {
      char entryname[50];
      sprintf(entryname,"OpDet%d_channels",iop);
      std::vector< unsigned int > chinput =  pset.get< std::vector<unsigned int> >( entryname );
      fPMT2channels[ iop ] = chinput;

      //std::cout << entryname << ": [";
      for (std::vector<unsigned int>::iterator it_ch=chinput.begin(); it_ch!=chinput.end(); it_ch++) {
	fChannel2pmt[ *it_ch ] = iop;
	//std::cout << *it_ch << ",";
	fNReadoutChannels++;
	if ( *it_ch > fMaxOpChannel )
	  fMaxOpChannel = *it_ch;
      }
      //std::cout << "]" << std::endl;
    }
    
  }
  
} // namespace
