module Geodesy

export
    # Points
    ECEF,
    ENU,
    LL,
    LLA,

    # Other types
    Bounds,
    Datum,
    Ellipsoid,

    # Constants
    WGS84,
    OSGB36,
    NAD27,

    # Methods
    center,
    distance,
    getX,
    getY,
    getZ,
    inBounds

    #= Unexported / Experimental
    decimal2dms
    dms2decimal

    haversine_distance

    boundaryPoint
    onBounds
    =#

for f in ["datum", "point", "bounds", "transform", "vicenty", "distance"]
    include("$f.jl")
end

end # module Geodesy
