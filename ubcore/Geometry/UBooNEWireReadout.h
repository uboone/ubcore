////////////////////////////////////////////////////////////////////////////////
/// \file UBooNEGeometryHelper.h
/// \brief Geometry helper service for UBooNE geometries.
///
/// Handles UBooNE-specific information for the generic Geometry service
/// within LArSoft. Derived from the WireReadout class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef UBooNE_WireReadout_h
#define UBooNE_WireReadout_h

#include "ubcore/Geometry/WireReadoutUBooNEGeom.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/fwd.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace uboone
{
  class UBooNEWireReadout : public geo::WireReadout
  {
  public:
    explicit UBooNEWireReadout(fhicl::ParameterSet const& pset);

  private:
    geo::WireReadoutGeom const& wireReadoutGeom() const override { return fWireReadout; }

    geo::WireReadoutUBooNEGeom fWireReadout;
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(uboone::UBooNEWireReadout, geo::WireReadout, SHARED)

#endif // UBooNE_WireReadout_h
