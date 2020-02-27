////////////////////////////////////////////////////////////////////////////////
/// \file UBooNEGeometryHelper.h
/// \brief Geometry helper service for UBooNE geometries.
///
/// Handles UBooNE-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef UBooNE_ExptGeoHelperInterface_h
#define UBooNE_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

#include <memory>

namespace geo
{
  class ChannelMapAlg;
}

namespace uboone
{
  class UBooNEGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
    explicit UBooNEGeometryHelper(fhicl::ParameterSet const& pset);

  private:
    ChannelMapAlgPtr_t
    doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                             std::string const& detectorName) const override;

    fhicl::ParameterSet fPset; ///< copy of configuration parameter set
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(uboone::UBooNEGeometryHelper,
                                   geo::ExptGeoHelperInterface,
                                   SHARED)

#endif // UBooNE_ExptGeoHelperInterface_h
