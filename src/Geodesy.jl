module Geodesy

using Compat

export
    # Points
    ECEF,
    ENU,
    LLA,

    # Other types
    Ellipsoid,

    # Constants
    WGS84,
    OSGB36,
    NAD27,

    # Methods
    distance,
    getX,
    getY,
    getZ,

include("points.jl")
include("conversion.jl")

end # module Geodesy
