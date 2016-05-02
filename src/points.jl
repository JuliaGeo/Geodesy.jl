###################
### Point Types ###
###################

"""
LLA(lat,lon,alt): Latitude, longitude, and alititude co-ordinates.

Note: assumes degrees not radians
"""
immutable LLA{T}
    lat::T
    lon::T
    alt::T
end
LLA(;lat=NaN,lon=NaN,alt=0.0) = LLA(lat,lon,alt) # Constructor that is idependent of storage order
LLA(lat, lon) = LLA(lat, lon, 0.0)
Base.show(io::IO, lla::LLA) = print(io, "LLA(lat=$(lla.lat)°, lon=$(lla.lon)°, alt=$(lla.alt))") # Maybe show to nearest mm, or 6 decimel places?
Base.isapprox(lla1::LLA, lla2::LLA; atol = 1e-6, kwargs...) = isapprox(lla1.lat, lla2.lat; atol = 180*atol/6.371e6, kwargs...) & isapprox(lla1.lon, lla2.lon; atol = 180*atol/6.371e6, kwargs...) & isapprox(lla1.alt, lla2.alt; atol = atol, kwargs...) # atol in metres (1μm)

"""
Latitude and longitude co-ordinates.
(Note: assumes degrees not radians)
"""
immutable LatLon{T}
    lat::T
    lon::T
end
LatLon(;lat=NaN,lon=NaN) = LatLon(lat,lon) # Constructor that is idependent of storage order
LatLon(lla::LLA) = LatLon(lla.lat, lla.lon)
function LatLon(x, datum) #?
    lla = LLA(x, datum)
    return LatLon(lla.lat, lla.lon)
end
LLA(ll::LatLon) = LLA(ll.lat, ll.lon)
Base.show(io::IO, lla::LatLon) = print(io, "LatLon(lat=$(lla.lat)°, lon=$(lla.lon)°)") # Maybe show to nearest mm, or 6 decimel places?
Base.isapprox(ll1::LatLon, ll2::LatLon; atol = 1e-6, kwargs...) = isapprox(ll1.lat, ll2.lat; atol = 180*atol/6.371e6, kwargs...) & isapprox(ll1.lon, ll2.lon; atol = 180*atol/6.371e6, kwargs...) # atol in metres (1μm)

"""
ECEF(x,y,z): Earth-Centered-Earth-Fixed (ECEF) coordinates

A global Cartesian coordinate system rotating with the Earth.
"""
immutable ECEF{T} <: FixedVectorNoTuple{3,Float64}
    x::T
    y::T
    z::T
end

"""
ENU(e,n,u): East-North-Up (ENU) coordinates

A local Cartesian coordinate system, linearized about a reference point.
"""
immutable ENU{T} <: FixedVectorNoTuple{3,Float64}
    e::T
    n::T
    u::T
end
ENU(x, y) = ENU(x, y, 0.0)
