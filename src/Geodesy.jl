module Geodesy

export
    # Points
    ECEF,
    ENU,
    LLA,

    # Other types
    Bounds,
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

for f in ["datum", "point", "bounds", "conversion"]
    include("$f.jl")
end

end # module Geodesy
