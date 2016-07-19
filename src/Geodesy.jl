__precompile__()

module Geodesy

using CoordinateTransformations
using FixedSizeArrays
using Compat

import CoordinateTransformations.transform_deriv,
       CoordinateTransformations.transform_deriv_params,
       CoordinateTransformations.compose,
       CoordinateTransformations.∘

export
    # Points
    ECEF,
    ENU,
    LLA,
    LatLon,
    UTM,
    UTMZ,

    # Other types
    Ellipsoid, ellipsoid,

    # Constants
    WGS84, wgs84,
    OSGB36, osgb36,
    NAD27, nad27,
    GRS80, grs80,

    # Methods
    distance,

    # transformation methods
    transform_deriv, transform_deriv_params, compose, ∘,
    ECEFfromLLA, LLAfromECEF, ENUfromECEF, ECEFfromENU, ENUfromLLA, LLAfromENU,
    UTMfromLLA, LLAfromUTM, UTMfromECEF, ECEFfromUTM, ENUfromUTM, UTMfromENU,
    UTMZfromLLA, LLAfromUTMZ, UTMZfromECEF, ECEFfromUTMZ, ENUfromUTMZ, UTMZfromENU,
    UTMZfromUTM, UTMfromUTMZ,

    # Datum transformations
    GDA94_from_ITRF, ITRF_from_GDA94,
    GDA94_from_ITRF2008, GDA94_from_ITRF2005, GDA94_from_ITRF2000, GDA94_from_ITRF1997, GDA94_from_ITRF1996,
    ITRF2008_from_GDA94, ITRF2005_from_GDA94, ITRF2000_from_GDA94, ITRF1997_from_GDA94, ITRF1996_from_GDA94,

    # UTM helpers
    utm_zone

include("points.jl")
include("datums.jl")
include("transverse_mercator.jl")
include("polar_stereographic.jl")
include("transformations.jl")
include("conversion.jl")
include("distances.jl")
include("utm.jl")
include("datum_transformations.jl")

end # module Geodesy
