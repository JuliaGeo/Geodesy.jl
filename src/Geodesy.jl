module Geodesy

using Dates
using CoordinateTransformations
using StaticArrays
using HTTP
using GeoArrays

import CoordinateTransformations.transform_deriv,
       CoordinateTransformations.transform_deriv_params,
       CoordinateTransformations.compose,
       CoordinateTransformations.∘

include("ellipsoids.jl")
include("datums.jl")
include("points.jl")
include("transverse_mercator.jl")
include("polar_stereographic.jl")
include("web_mercator.jl")
include("transformations.jl")
include("conversion.jl")
include("distances.jl")
include("utm.jl")
include("datum_transformations.jl")
include("geoid.jl")

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

    # Ellipsoids
    wgs84_ellipsoid, airy1830, clarke1866, grs80,

    # Datums
    Datum,
    WGS84, wgs84,
    OSGB36, osgb36,
    NAD27, nad27,
    NAD83, nad83,
    GRS80, grs80,
    GDA94, gda94,

    # Methods
    euclidean_distance,

    # transformation methods
    transform_deriv, transform_deriv_params, compose, ∘,
    ECEFfromLLA, LLAfromECEF, ENUfromECEF, ECEFfromENU, ENUfromLLA, LLAfromENU,
    UTMfromLLA, LLAfromUTM, UTMfromECEF, ECEFfromUTM, ENUfromUTM, UTMfromENU,
    UTMZfromLLA, LLAfromUTMZ, UTMZfromECEF, ECEFfromUTMZ, ENUfromUTMZ, UTMZfromENU,
    UTMZfromUTM, UTMfromUTMZ,
    WebMercatorfromLLA, LLAfromWebMercator, geoid

    # Datum transformations
    datum_shift_ECEF,
    ITRF, GDA94,

    # UTM helpers
    utm_zone

end # module Geodesy
