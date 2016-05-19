##################
## LLA <-> ECEF ##
##################

immutable LLAfromECEF{Datum} <: AbstractTransformation{LLA, ECEF}
    a::Float64      # major axis
    f::Float64      # flattening
    e2::Float64     # Eccentricity squared
    e2m::Float64    # 1 - e2
    e2a::Float64    # |e2|
    e4a::Float64    # e2^2

    datum::Datum

    function LLAfromECEF(a,f,e2,e2m,e2a,e4a,datum)
        if !(isfinite(a) && a > 0)
            error("Major radius is not positive")
        end
        if !(isfinite(f) && f < 1)
            error("Minor radius is not positive")
        end
        return new(a,f,e2,e2m,e2a,e4a,datum)
    end
end
Base.show(io::IO, trans::LLAfromECEF) = print(io, "LLAfromECEF($(trans.datum))")

"""
    LLAfromECEF(datum)

Construct a `AbstractTransformation` object to convert from ECEF coordinates
to LLA coordinates. Pre-caches ellipsoidal parameters for efficiency.
"""
function LLAfromECEF{Datum}(datum::Datum)
    el = ellipsoid(datum)

    a = el.a
    b = el.b
    f = 1 - b/a
    e2 = f*(2-f) # or el.e²
    e2m = (1-f)*(1-f)  #1 - e2
    e2a = abs(e2)
    e4a = e2*e2

    return LLAfromECEF{Datum}(a, f, e2, e2m, e2a, e4a, datum)
end


function transform(trans::LLAfromECEF, ecef::ECEF)
    # Ported to Julia by Andy Ferris, 2016 and re-released under MIT license.
    #/**
    # * \file Geocentric.cpp
    # * \brief Implementation for GeographicLib::Geocentric class
    # *
    # * Copyright (c) Charles Karney (2008-2015) <charles@karney.com> and licensed
    # * under the MIT/X11 License.  For more information, see
    # * http://geographiclib.sourceforge.net/
    # **********************************************************************/
    R = hypot(ecef.x, ecef.y)
    if R == 0
        slam = 0.0
        clam = 1.0
    else
        slam = ecef.y / R
        clam = ecef.x / R
    end
    h = hypot(R, ecef.z)    # Distance to center of earth

    if (trans.e4a == 0)
        # Treat the spherical case.  Dealing with underflow in the general case
        # with _e2 = 0 is difficult.  Origin maps to north pole, same as the
        # ellipsoidal case below
        if h == 0
            sphi = 1.0
            cphi = 0.0
        else
            sphi = ecef.z / h
            cphi = R / h
        end
        h -= trans.a
    else # Ellipsoidal
        # Treat prolate spheroids by swapping r and z here and by switching
        # the arguments to phi = atan2(...) at the end.
        p = (R / trans.a) * (R / trans.a)
        q = trans.e2m * (ecef.z / trans.a) * (ecef.z / trans.a)
        r = (p + q - trans.e4a) / 6
        if (trans.f < 0)
            tmp = p
            p = q
            q = tmp
        end


        if ( !(trans.e4a * q == 0 && r <= 0) )

            # Avoid possible division by zero when r = 0 by multiplying
            # equations for s and t by r^3 and r, resp.
            S = trans.e4a * p * q / 4 # S = r^3 * s
            r2 = r * r
            r3 = r * r2
            disc = S * (2 * r3 + S)
            u = r
            if (disc >= 0)
                T3 = S + r3
                # Pick the sign on the sqrt to maximize abs(T3).  This minimizes
                # loss of precision due to cancellation.  The result is unchanged
                # because of the way the T is used in definition of u.
                T3 += (T3 < 0 ? -sqrt(disc) : sqrt(disc)) # T3 = (r * t)^3
                # N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                T = cbrt(T3) # T = r * t
                # T can be zero; but then r2 / T -> 0.
                u += T + (T != 0 ? r2 / T : 0.0)
            else
                # T is complex, but the way u is defined the result is real.
                ang = atan2(sqrt(-disc), -(S + r3))
                # There are three possible cube roots.  We choose the root which
                # avoids cancellation.  Note that disc < 0 implies that r < 0.
                u += 2 * r * cos(ang / 3)
            end

            v = sqrt(u*u + trans.e4a * q) # guaranteed positive
            # Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
            # e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
            uv = (u < 0 ? trans.e4a * q / (v - u) : u + v) # u+v, guaranteed positive
            # Need to guard against w going negative due to roundoff in uv - q.
            w = max(0.0, trans.e2a * (uv - q) / (2 * v))
            # Rearrange expression for k to avoid loss of accuracy due to
            # subtraction.  Division by 0 not possible because uv > 0, w >= 0.
            k = uv / (sqrt(uv + w*w) + w)
            k1 = (trans.f >= 0 ? k : k - trans.e2)
            k2 = (trans.f >= 0 ? k + trans.e2 : k)
            d = k1 * R / k2
            H = hypot(ecef.z/k1, R/k2)
            sphi = (ecef.z/k1) / H
            cphi = (R/k2) / H
            h = (1 - trans.e2m/k1) * hypot(d, ecef.z)
        else  # e4 * q == 0 && r <= 0
            # This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
            # (prolate, rotation axis) and the generation of 0/0 in the general
            # formulas for phi and h.  using the general formula and division by 0
            # in formula for h.  So handle this case by taking the limits:
            # f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e4 - p)
            # f < 0: R -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e4 - p)
            zz = sqrt((trans.f >= 0 ? trans.e4a - p : p) / trans.e2m)
            xx = sqrt( trans.f <  0 ? trans.e4a - p : p)
            H = hypot(zz, xx)
            sphi = zz / H
            cphi = xx / H
            if (Z < 0)
                sphi = -sphi # for tiny negative Z (not for prolate)
            end
            h = - trans.a * (trans.f >= 0 ? trans.e2m : 1.0) * H / trans.e2a
        end
    end
    lat = 180 * atan2(sphi, cphi) / pi
    lon = 180 * atan2(slam, clam) / pi

    return LLA(lat, lon, h)
end


"""
    ECEFfromLLA(datum)

Construct a `AbstractTransformation` object to convert from LLA coordinates
to ECEF coordinates.
"""
immutable ECEFfromLLA{Datum} <: AbstractTransformation{ECEF, LLA}
    a::Float64   # Ellipsoidal major axis
    e²::Float64  # Ellipsoidal square-eccentricity = 1 - b^2/a^2

    datum::Datum
end
Base.show(io::IO, trans::ECEFfromLLA) = print(io, "ECEFfromLLA($(trans.datum))")

function ECEFfromLLA{Datum}(datum::Datum)
    el = ellipsoid(datum)
    return ECEFfromLLA{Datum}(el.a, el.e², datum)
end


function transform(trans::ECEFfromLLA, lla::LLA)
    ϕdeg, λdeg, h = lla.lat, lla.lon, lla.alt

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = trans.a / sqrt(1 - trans.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - trans.e²) + h) * sinϕ

    return ECEF(x, y, z)
end

Base.inv(trans::LLAfromECEF) = ECEFfromLLA(trans.datum)
Base.inv(trans::ECEFfromLLA) = LLAfromECEF(trans.datum)

# It's not clear if this is worthwhile or not... (return is not type-certain)
compose(trans1::ECEFfromLLA, trans2::LLAfromECEF) = t1.datum === t2.datum ? IdentityTransform{ECEF}() : ComposedTransformation{outtype(trans1), intype(trans2), typeof(trans1), typeof(trans2)}(trans1, trans2)
compose(trans1::LLAfromECEF, trans2::ECEFfromLLA) = t1.datum === t2.datum ? IdentityTransform{LLA}() : ComposedTransformation{outtype(trans1), intype(trans2), typeof(trans1), typeof(trans2)}(trans1, trans2)

##################
## ECEF <-> ENU ##
##################

"""
    ENUfromECEF(origin, datum)
    ENUfromECEF(origin::ECEF, lat, lon)

Construct a `AbstractTransformation` object to convert from global `ECEF` coordinates
to local `ENU` coordinates centered at the `origin`. This object pre-caches both the
ECEF coordinates and latitude and longitude of the origin for maximal efficiency.
"""
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
Base.show(io::IO, trans::ENUfromECEF) = print(io, "ENUfromECEF($(trans.origin), lat=$(trans.lat)°, lon=$(trans.lon)°)")
Base.isapprox(t1::ENUfromECEF, t2::ENUfromECEF; kwargs...) = isapprox(t1.origin, t2.origin; kwargs...) && isapprox(t1.lat, t2.lat; kwargs...) && isapprox(t1.lon, t2.lon; kwargs...)


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

"""
    ECEFfromENU(origin, datum)
    ECEFfromENU(origin::ECEF, lat, lon)

Construct a `AbstractTransformation` object to convert from local `ENU` coordinates
centred at `origin` to global `ECEF` coodinates. This object pre-caches both the
ECEF coordinates and latitude and longitude of the origin for maximal efficiency.
"""
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
Base.show(io::IO, trans::ECEFfromENU) = print(io, "ECEFfromENU($(trans.origin), lat=$(trans.lat)°, lon=$(trans.lon)°)")
Base.isapprox(t1::ECEFfromENU, t2::ECEFfromENU; kwargs...) = isapprox(t1.origin, t2.origin; kwargs...) && isapprox(t1.lat, t2.lat; kwargs...) && isapprox(t1.lon, t2.lon; kwargs...)

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

Base.inv(trans::ECEFfromENU) = ENUfromECEF(trans.origin, trans.lat, trans.lon)
Base.inv(trans::ENUfromECEF) = ECEFfromENU(trans.origin, trans.lat, trans.lon)

#################
## LLA <-> ENU ##
#################

ENUfromLLA(origin, datum) = ENUfromECEF(origin, datum) ∘ ECEFfromLLA(datum)

LLAfromENU(origin, datum) = LLAfromECEF(datum) ∘ ECEFfromENU(origin, datum)

#################
## LLA <-> UTM ##
#################

"""
    LLAfromUTM(zone, northern_hemisphere::Bool, datum)

Construct a `AbstractTransformation` object to convert from `UTM` coordinates in
the specified zone and hemisphere to global `LLA` coordinates. Pre-caches
ellipsoidal parameters for efficiency and performs Charles Karney's accurate
6th-order series expansion algorithm.

(See also `LLAfromUTMZ`)
"""
immutable LLAfromUTM{Datum,Order} <: AbstractTransformation{LLA, UTM}
    zone::UInt8
    hemisphere::Bool # true = north, false = south
    tm::TransverseMercator{Order}
    datum::Datum
end
LLAfromUTM(zone::UInt8, h, d) = LLAfromUTM(UInt8(zone), h, TransverseMercator(d), d)
LLAfromUTM(zone::Integer, h, d) = LLAfromUTM(UInt8(zone), h, d)
Base.show(io::IO, trans::LLAfromUTM) = print(io, "LLAfromUTM(zone=$(trans.zone == 0 ? "polar" : trans.zone) ($(trans.hemisphere ? "north" : "south")), $(trans.datum))")

function transform(trans::LLAfromUTM, utm::UTM)
    if trans.zone == 0
        # Do inverse steriographic projection
        k0 = 0.994
        x = (utm.x - 2e6)
        y = (utm.y - 2e6)
        (lat,lon,γ,k) = polarst_inv(trans.hemisphere, k0, trans.tm, x, y)

        return LLA(lat, lon, utm.z)
    else
        lon_ref = Float64(utm_meridian(trans.zone))
        k0 = 0.9996 # Horizontal scaling factor
        x = (utm.x - 5e5) # Convention has 500km offset for easting
        y = (utm.y - (trans.hemisphere ? 0.0 : 1e7)) # Northing offset for southern hemisphere
        #(lat,lon,k,γ) = transverse_mercator_inv(lat_ref, lon_ref, x, y, ellipsoid(trans.datum))
        (lat,lon,γ,k) = transverse_mercator_reverse(lon_ref, x, y, k0, trans.tm)
        return LLA(lat, lon, utm.z) # Note: scaling not applied to vertical dimension
    end
end

"""
    UTMfromLLA(zone, northern_hemisphere::Bool, datum)

Construct a `AbstractTransformation` object to convert from global `LLA` coordinates
to `UTM` coordinates in the specified zone and hemisphere. Pre-caches
ellipsoidal parameters for efficiency and performs Charles Karney's accurate
6th-order series expansion algorithm.

(See also `UTMZfromLLA`)
"""
immutable UTMfromLLA{Datum,Order} <: AbstractTransformation{UTM, LLA}
    zone::UInt8
    hemisphere::Bool # true = north, false = south
    tm::TransverseMercator{Order}
    datum::Datum
end
UTMfromLLA(zone::UInt8, h, d) = UTMfromLLA(zone, h, TransverseMercator(d), d)
UTMfromLLA(zone::Integer, h, d) = UTMfromLLA(UInt8(zone), h, d)
Base.show(io::IO, trans::UTMfromLLA) = print(io, "UTMfromLLA(zone=$(trans.zone == 0 ? "polar" : trans.zone) ($(trans.hemisphere ? "north" : "south")), $(trans.datum))")

function transform(trans::UTMfromLLA, lla::LLA)
    if trans.zone == 0
        # Do polar steriographic projection
        k0 = 0.994
        (x,y,γ,k) = polarst_fwd(trans.hemisphere, k0, trans.tm, lla.lat, lla.lon)
        x = x + 2e6
        y = y + 2e6

        return UTM(x, y, lla.alt)
    else
        lon_ref = Float64(utm_meridian(trans.zone))
        k0 = 0.9996 # Horizontal scaling factor
        #(x,y,k,γ) = transverse_mercator(lat_ref, lon_ref, lla.lat, lla.lon, ellipsoid(trans.datum))
        (x,y,γ,k) = transverse_mercator_forward(lon_ref, lla.lat, lla.lon, k0, trans.tm)
        x = 5e5 + x # Convention has 500km offset for easting
        y = (trans.hemisphere ? 0.0 : 1e7) + y # Northing offset for southern hemisphere
        # also, k = k * k0
        return UTM(x, y, lla.alt) # Note: scaling not applied to vertical dimension
    end
end

Base.inv(trans::LLAfromUTM) = UTMfromLLA(trans.zone, trans.hemisphere, trans.tm, trans.datum)
Base.inv(trans::UTMfromLLA) = LLAfromUTM(trans.zone, trans.hemisphere, trans.tm, trans.datum)


##################
## ECEF <-> UTM ##
##################

UTMfromECEF(zone, hemisphere, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromECEF(datum)
ECEFfromUTM(zone, hemisphere, datum) = ECEFfromLLA(datum) ∘ LLAfromUTM(zone, hemisphere, datum)

##################
## LLA <-> UTMZ ##
##################

"""
    LLAfromUTMZ(datum)

Construct a `AbstractTransformation` object to convert from global `UTMZ`
coordinates to global `LLA` coordinates. Pre-caches
ellipsoidal parameters for efficiency and performs Charles Karney's accurate
6th-order series expansion algorithm.

(See also `LLAfromUTM`)
"""
immutable LLAfromUTMZ{Datum,Order} <: AbstractTransformation{LLA, UTMZ}
    tm::TransverseMercator{Order}
    datum::Datum
end
LLAfromUTMZ(datum) = LLAfromUTMZ(TransverseMercator(datum), datum)
Base.show(io::IO, trans::LLAfromUTMZ) = print(io, "LLAfromUTMZ($(trans.datum))")


function transform(trans::LLAfromUTMZ, utm::UTMZ)
    if utm.zone == 0
        # Do inverse steriographic projection
        k0 = 0.994
        x = (utm.x - 2e6)
        y = (utm.y - 2e6)
        (lat,lon,γ,k) = polarst_inv(utm.hemisphere, k0, trans.tm, x, y)

        return LLA(lat, lon, utm.z)
    else
        lon_ref = Float64(utm_meridian(utm.zone))
        k0 = 0.9996 # Horizontal scaling factor
        x = (utm.x - 5e5) # Convention has 500km offset for easting
        y = (utm.y - (utm.hemisphere ? 0.0 : 1e7)) # Northing offset for southern hemisphere
        #(lat,lon,k,γ) = transverse_mercator_inv(lat_ref, lon_ref, x, y, ellipsoid(trans.datum))
        (lat,lon,γ,k) = transverse_mercator_reverse(lon_ref, x, y, k0, trans.tm)
        return LLA(lat, lon, utm.z) # Note: scaling not applied to vertical dimension
    end
end

"""
    UTMZfromLLA(zone, northern_hemisphere::Bool, datum)

Construct a `AbstractTransformation` object to convert from global `LLA` coordinates
to global `UTMZ` coordinates. The zone and hemisphere is automatically calculated
following the standard definitions (including exceptions in Norway). Pre-caches
ellipsoidal parameters for efficiency and performs Charles Karney's accurate
6th-order series expansion algorithm.

(See also `UTMfromLLA`)
"""
immutable UTMZfromLLA{Datum,Order} <: AbstractTransformation{UTMZ, LLA}
    tm::TransverseMercator{Order}
    datum::Datum
end
UTMZfromLLA(datum) = UTMZfromLLA(TransverseMercator(datum), datum)
Base.show(io::IO, trans::UTMZfromLLA) = print(io, "UTMZfromLLA($(trans.datum))")


function transform(trans::UTMZfromLLA, lla::LLA)
    (zone, hemisphere) = utm_zone(lla)
    zone::Int64
    hemisphere::Bool
    if zone == 0
        # Do polar steriographic projection
        k0 = 0.994
        (x,y,γ,k) = polarst_fwd(hemisphere, k0, trans.tm, lla.lat, lla.lon)
        x = x + 2e6
        y = y + 2e6

        return UTMZ(x, y, lla.alt, 0, hemisphere)
    else
        lon_ref = Float64(utm_meridian(zone))
        k0 = 0.9996 # Horizontal scaling factor
        #(x,y,k,γ) = transverse_mercator(lat_ref, lon_ref, lla.lat, lla.lon, ellipsoid(trans.datum))
        (x,y,γ,k) = transverse_mercator_forward(lon_ref, lla.lat, lla.lon, k0, trans.tm)
        x = 5e5 + x # Convention has 500km offset for easting
        y = (hemisphere ? 0.0 : 1e7) + y # Northing offset for southern hemisphere
        # also, k = k * k0
        return UTMZ(x, y, lla.alt, zone, hemisphere) # Note: scaling not applied to vertical dimension
    end
end

Base.inv(trans::LLAfromUTMZ) = UTMZfromLLA(trans.tm, trans.datum)
Base.inv(trans::UTMZfromLLA) = LLAfromUTMZ(trans.tm, trans.datum)

###################
## UTM <-> UTMZ ##
###################

immutable UTMZfromUTM{Datum} <: AbstractTransformation{UTMZ, UTM}
    zone::Int
    hemisphere::Bool
    datum::Datum
end
UTMZfromUTM{D}(zone::Integer, h, d::D) = UTMZfromUTM{D}(UInt8(zone), h, d)
Base.show(io::IO, trans::UTMZfromUTM) = print(io, "UTMZfromUTM(zone=$(trans.zone == 0 ? "polar" : trans.zone) ($(trans.hemisphere ? "north" : "south")), $(trans.datum))")
#Base.show(io::IO, trans::UTMZfromUTM) = print(io, "UTMZfromUTM(zone=$(trans.zone == 0 ? "polar" : trans.zone) ($(trans.hemisphere ? "north" : "south")))")

immutable UTMfromUTMZ{Datum} <: AbstractTransformation{UTM, UTMZ}
    zone::Int
    hemisphere::Bool
    datum::Datum
end
UTMfromUTMZ{D}(zone::Integer, h, d::D) = UTMfromUTMZ{D}(UInt8(zone), h, d)
Base.show(io::IO, trans::UTMfromUTMZ) = print(io, "UTMfromUTMZ(zone=$(trans.zone == 0 ? "polar" : trans.zone) ($(trans.hemisphere ? "north" : "south")), $(trans.datum))")
#Base.show(io::IO, trans::UTMfromUTMZ) = print(io, "UTMfromUTMZ(zone=$(trans.zone == 0 ? "polar" : trans.zone) ($(trans.hemisphere ? "north" : "south")))")

transform(trans::UTMZfromUTM, utm::UTM) = UTMZ(utm.x, utm.y, utm.z, trans.zone, trans.hemisphere)
function transform(trans::UTMfromUTMZ, utm::UTMZ)
    if trans.zone == utm.zone && trans.hemisphere == utm.hemisphere
        UTM(utm.x, utm.y, utm.z)
    else
        # Should this be an error or an automatic transformation to the correct zone?
        #error("Incorrect UTM zone")
        transform(UTMfromLLA(trans.zone, trans.hemisphere, trans.datum), transform(LLAfromUTMZ(trans.datum), utm))
    end
end

Base.inv(trans::UTMfromUTMZ) = UTMZfromUTM(trans.zone, trans.hemisphere, trans.datum)
Base.inv(trans::UTMZfromUTM) = UTMfromUTMZ(trans.zone, trans.hemisphere, trans.datum)

###################
## ECEF <-> UTMZ ##
###################

UTMZfromECEF(datum) = UTMZfromLLA(datum) ∘ LLAfromECEF(datum)
ECEFfromUTMZ(datum) = ECEFfromLLA(datum) ∘ LLAfromUTMZ(datum)

##################
## ENU <-> UTMZ ##
##################

ENUfromECEF(origin::UTMZ, datum) = ENUfromECEF(transform(LLAfromUTMZ(datum), origin), datum)
ECEFfromENU(origin::UTMZ, datum) = ECEFfromENU(transform(LLAfromUTMZ(datum), origin), datum)

ENUfromUTMZ(origin, datum)  = ENUfromLLA(origin, datum) ∘ LLAfromUTMZ(datum)
UTMZfromENU(origin, datum)  = UTMZfromLLA(datum) ∘ LLAfromENU(origin, datum)

#################
## ENU <-> UTM ##
#################
ENUfromECEF(origin::UTM, zone::Integer, hemisphere::Bool, datum) = ENUfromECEF(transform(LLAfromUTM(zone, hemisphere, datum), origin), datum)
ECEFfromENU(origin::UTM, zone::Integer, hemisphere::Bool, datum) = ECEFfromENU(transform(LLAfromUTM(zone, hemisphere, datum), origin), datum)

# Assume origin and utm point share the same zone and hemisphere
UTMfromENU(origin::UTM, zone::Integer, hemisphere::Bool, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromENU(UTMZ(origin, zone, hemisphere), datum)
ENUfromUTM(origin::UTM, zone::Integer, hemisphere::Bool, datum) = ENUfromLLA(UTMZ(origin, zone, hemisphere), datum) ∘ LLAfromUTM(zone, hemisphere, datum)

UTMfromENU(origin, zone::Integer, hemisphere::Bool, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromENU(origin, datum)
ENUfromUTM(origin, zone::Integer, hemisphere::Bool, datum) = ENUfromLLA(origin, datum) ∘ LLAfromUTM(zone, hemisphere, datum)
