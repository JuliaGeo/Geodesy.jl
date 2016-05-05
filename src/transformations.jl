##################
## LLA <-> ECEF ##
##################

immutable LLAfromECEF{Datum} <: AbstractTransformation{LLA, ECEF}
    datum::Datum
end

function transform(trans::LLAfromECEF, ecef::ECEF)
    x, y, z = ecef.x, ecef.y, ecef.z
    d = ellipsoid(trans.datum)

    p = hypot(x, y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e′² * d.b * sin(θ)^3, p - d.e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    return LLA(rad2deg(ϕ), rad2deg(λ), h)
end

immutable ECEFfromLLA{Datum} <: AbstractTransformation{ECEF, LLA}
    datum::Datum
end

function transform(trans::ECEFfromLLA, lla::LLA)
    ϕdeg, λdeg, h = lla.lat, lla.lon, lla.alt
    d = ellipsoid(trans.datum)

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = d.a / sqrt(1 - d.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - d.e²) + h) * sinϕ

    return ECEF(x, y, z)
end

Base.inv(trans::LLAfromECEF) = ECEFfromLLA(trans.datum)
Base.inv(trans::ECEFfromLLA) = LLAfromECEF(trans.datum)
compose(trans1::ECEFfromLLA, trans2::LLAfromECEF) = t1.datum === t2.datum ? IdentityTransform{ECEF}() : ComposedTransformation{outtype(trans1), intype(trans2), typeof(trans1), typeof(trans2)}(trans1, trans2)
compose(trans1::LLAfromECEF, trans2::ECEFfromLLA) = t1.datum === t2.datum ? IdentityTransform{LLA}() : ComposedTransformation{outtype(trans1), intype(trans2), typeof(trans1), typeof(trans2)}(trans1, trans2)

##################
## ECEF <-> ENU ##
##################

immutable ENUfromECEF{T} <: AbstractTransformation{ENU, ECEF}
    origin::ECEF{T}
    lat::T
    lon::T
end
ENUfromECEF(origin::LLA, datum) = ENUfromECEF(transform(ECEFfromLLA(datum),origin), origin.lat, origin.lon)
function ENUfromECEF(origin::ECEF, datum)
    origin_lla = transform(LLAfromECEF(datum), origin)
    ENUfromECEF(origin, origin_lla.lat, origin_lla.lon)
end

function transform(trans::ENUfromECEF, ecef::ECEF)
    ϕdeg, λdeg = trans.lat, trans.lon

    ∂x = ecef.x - trans.origin.x
    ∂y = ecef.y - trans.origin.y
    ∂z = ecef.z - trans.origin.z

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # R = [     -sinλ       cosλ  0.0
    #      -cosλ*sinϕ -sinλ*sinϕ cosϕ
    #       cosλ*cosϕ  sinλ*cosϕ sinϕ]
    #
    # east, north, up = R * [∂x, ∂y, ∂z]
    east  = ∂x * -sinλ      + ∂y * cosλ       + ∂z * 0.0
    north = ∂x * -cosλ*sinϕ + ∂y * -sinλ*sinϕ + ∂z * cosϕ
    up    = ∂x * cosλ*cosϕ  + ∂y * sinλ*cosϕ  + ∂z * sinϕ

    return ENU(east, north, up)
end

immutable ECEFfromENU{T} <: AbstractTransformation{ECEF, ENU}
    origin::ECEF{T}
    lat::T
    lon::T
end
ECEFfromENU(origin::LLA, datum) = ECEFfromENU(transform(ECEFfromLLA(datum),origin), origin.lat, origin.lon)
function ECEFfromENU(origin::ECEF, datum)
    origin_lla = transform(LLAfromECEF(datum), origin)
    ECEFfromENU(origin, origin_lla.lat, origin_lla.lon)
end


function transform(trans::ECEFfromENU, enu::ENU)
    ϕdeg, λdeg = trans.lat, trans.lon

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Rᵀ = [-sinλ -cosλ*sinϕ cosλ*cosϕ
    #        cosλ -sinλ*sinϕ sinλ*cosϕ
    #         0.0       cosϕ      sinϕ]
    # Δx, Δy, Δz = Rᵀ * [east, north, up]
    Δx = -sinλ * enu.e + -cosλ*sinϕ * enu.n + cosλ*cosϕ * enu.u
    Δy =  cosλ * enu.e + -sinλ*sinϕ * enu.n + sinλ*cosϕ * enu.u
    Δz =   0.0 * enu.e +       cosϕ * enu.n +      sinϕ * enu.u

    X = trans.origin.x + Δx
    Y = trans.origin.y + Δy
    Z = trans.origin.z + Δz

    return ECEF(X,Y,Z)
end

#################
## LLA <-> ENU ##
#################

ENUfromLLA(origin, datum) = ENUfromECEF(origin, datum) ∘ ECEFfromLLA(datum)

LLAfromENU(origin, datum) = LLAfromECEF(datum) ∘ ECEFfromENU(origin, datum)

#################
## LLA <-> UTM ##
#################

immutable LLAfromUTM{Datum} <: AbstractTransformation{LLA, UTM}
    zone::UInt8
    hemisphere::Bool # true = north, false = south
    datum::Datum
end
LLAfromUTM(zone::Integer, h, d) = LLAfromUTM(UInt8(zone), h, d)

function transform(trans::LLAfromUTM, utm::UTM)
    if trans.zone == 0
        # Do polar steriographic projection
        error("Zone 0 not implemented")
    else
        lat_ref = 0.0
        lon_ref = utm_meridian(trans.zone)
        k0 = 0.9996 # Horizontal scaling factor
        x = (utm.x - 5e5) / k0 # Convention has 500km offset for easting
        y = (utm.y - (trans.hemisphere ? 0.0 : 1e7)) / k0 # Northing offset for southern hemisphere
        (lat,lon,k,γ) = transverse_mercator_inv(lat_ref, lon_ref, x, y, ellipsoid(trans.datum))
        return LLA(lat, lon, utm.z) # Note: scaling not applied to vertical dimension
    end
end

immutable UTMfromLLA{Datum} <: AbstractTransformation{UTM, LLA}
    zone::UInt8
    hemisphere::Bool # true = north, false = south
    datum::Datum
end
UTMfromLLA(zone::Integer, h, d) = UTMfromLLA(UInt8(zone), h, d)

function transform(trans::UTMfromLLA, lla::LLA)
    if trans.zone == 0
        # Do polar steriographic projection
        error("Zone 0 not implemented")
    else
        lat_ref = 0.0
        lon_ref = utm_meridian(trans.zone)
        (x,y,k,γ) = transverse_mercator(lat_ref, lon_ref, lla.lat, lla.lon, ellipsoid(trans.datum))
        k0 = 0.9996 # Horizontal scaling factor
        x = 5e5 + k0*x # Convention has 500km offset for easting
        y = (trans.hemisphere ? 0.0 : 1e7) + k0*y # Northing offset for southern hemisphere
        # also, k = k * k0
        return UTM(x, y, lla.alt) # Note: scaling not applied to vertical dimension
    end
end

Base.inv(trans::LLAfromUTM) = UTMfromLLA(trans.zone, trans.hemisphere, trans.datum)
Base.inv(trans::UTMfromLLA) = LLAfromUTM(trans.zone, trans.hemisphere, trans.datum)


##################
## ECEF <-> UTM ##
##################

UTMfromECEF(zone, hemisphere, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromECEF(datum)
ECEFfromUTM(zone, hemisphere, datum) = ECEFfromLLA(datum) ∘ LLAfromUTM(zone, hemisphere, datum)

##################
## LLA <-> UTMZ ##
##################

immutable LLAfromUTMZ{Datum} <: AbstractTransformation{LLA, UTMZ}
    datum::Datum
end

function transform(trans::LLAfromUTMZ, utm::UTMZ)
    if utm.zone == 0
        # Do polar steriographic projection
        error("Zone 0 not implemented")
    else
        lat_ref = 0.0
        lon_ref = utm_meridian(utm.zone)
        k0 = 0.9996 # Horizontal scaling factor
        x = (utm.x - 5e5) / k0 # Convention has 500km offset for easting
        y = (utm.y - (utm.hemisphere ? 0.0 : 1e7)) / k0 # Northing offset for southern hemisphere
        (lat,lon,k,γ) = transverse_mercator_inv(lat_ref, lon_ref, x, y, ellipsoid(trans.datum))
        return LLA(lat, lon, utm.z) # Note: scaling not applied to vertical dimension
    end
end

immutable UTMZfromLLA{Datum} <: AbstractTransformation{UTMZ, LLA}
    datum::Datum
end

function transform(trans::UTMZfromLLA, lla::LLA)
    (zone, hemisphere) = utm_zone(lla)
    if zone == 0
        # Do polar steriographic projection
        error("Zone 0 not implemented")
    else
        lat_ref = 0.0
        lon_ref = utm_meridian(zone)
        (x,y,k,γ) = transverse_mercator(lat_ref, lon_ref, lla.lat, lla.lon, ellipsoid(trans.datum))
        k0 = 0.9996 # Horizontal scaling factor
        x = 5e5 + k0*x # Convention has 500km offset for easting
        y = (hemisphere ? 0.0 : 1e7) + k0*y # Northing offset for southern hemisphere
        # also, k = k * k0
        return UTMZ(x, y, lla.alt, zone, hemisphere) # Note: scaling not applied to vertical dimension
    end
end

Base.inv(trans::LLAfromUTMZ) = UTMZfromLLA(trans.datum)
Base.inv(trans::UTMZfromLLA) = LLAfromUTMZ(trans.datum)

###################
## ECEF <-> UTMZ ##
###################

UTMZfromECEF(datum) = UTMZfromLLA(datum) ∘ LLAfromECEF(datum)
ECEFfromUTMZ(datum) = ECEFfromLLA(datum) ∘ LLAfromUTMZ(datum)

##################
## ENU <-> UTMZ ##
##################

ENUfromECEF(origin::UTMZ, datum) = ENUfromECEF(transform(ECEFfromUTMZ(datum),origin), origin.lat, origin.lon)
ECEFfromENU(origin::UTMZ, datum) = ECEFfromENU(transform(ECEFfromUTMZ(datum),origin), origin.lat, origin.lon)

ENUfromUTMZ(origin, datum)  = ENUfromLLA(origin, datum) ∘ LLAfromUTMZ(datum)

UTMZfromENU(origin, datum)  = UTMZfromLLA(datum) ∘ LLAfromENU(origin, datum)
