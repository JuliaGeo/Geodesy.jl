
###################
### Point Types ###
###################

### Point in Latitude-Longitude-Altitude (LLA) coordinates
# Used to store node data in OpenStreetMap XML files
"""
Latitude, longitude, and alititude co-ordinate system.
(Note: assumes degrees not radians)
"""
immutable LLA
    lat::Float64
    lon::Float64
    alt::Float64
end
LLA(lat, lon) = LLA(lat, lon, 0.0)

"""
Latitude and longitude co-ordinates.
(Note: assumes degrees not radians)
"""
immutable LatLon
    lat::Float64
    lon::Float64
end
LatLon(lla::LLA) = LatLon(lla.lat, lla.lon)
function LatLon(x, datum)
    lla = LLA(x, datum)
    return LatLon(lla.lat, lla.lon)
end

"""
Point in Earth-Centered-Earth-Fixed (ECEF) coordinates
Global cartesian coordinate system rotating with the Earth
"""
immutable ECEF # <: FixedVectorNoTuple{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

"""
Point in East-North-Up (ENU) coordinates
Local cartesian coordinate system
Linearized about a reference point
"""
immutable ENU # <: FixedVectorNoTuple{3,Float64}
    e::Float64
    n::Float64
    u::Float64
end
ENU(x, y) = ENU(x, y, 0.0)


### distance
# Point translators
distance(a::ENU, b::ENU) = distance(a.e, a.n, a.u,
                                    b.e, b.n, b.u)

distance(a::ECEF, b::ECEF) = distance(a.x, a.y, a.z,
                                      b.x, b.y, b.z)

function distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
end
