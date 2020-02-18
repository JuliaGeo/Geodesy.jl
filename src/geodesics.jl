# Great circle distances

include("GeographicLib/GeographicLib.jl")
import .GeographicLib

"""
    GreatCircle

A `GreatCircle` object stores a cache of values needed to quickly compute
great circles.
"""
struct GreatCircle
    geod::GeographicLib.Geodesic
    datum::Datum
end

"""
    GreatCircle(datum) -> gc

Construct a `GreatCircle` `gc` from a datum.

The object `gc` can then be used as a function to calculate great circle properties.
"""
GreatCircle(datum::Datum) =
    (el = ellipsoid(datum); GreatCircle(GeographicLib.Geodesic(el.a, el.f), datum))

"""
    (::GreatCircle)(a, b) -> (azi, baz, dist, angle)

Compute the great circle between points `a` and `b`.  The function returns a named
tuple giving the forward azimuth (`azi`, °) from `a` to `b`,
the backazimuth (`baz`, °) from `b` to `a`, the great circle distance (`dist`, m)
and the spherical equivalent angular distance (`angle`, °).

#### Example

Compute the distance `dist` (m) from Nelson's Column, Trafalgar Square, London to
the Empire State Building, New York, USA, as well as the forward azimuth (`azi`°)
from London, the backazimuth (`baz`°) from New York and the spherical equivalent angular
distance (`angle`°).

```
julia> using Geodesy

julia> gc = GreatCircle(wgs84);

julia> nelson = LatLon(51.5077, -0.1277); esb = LatLon(40.7484, -73.9857);

julia> gc(nelson, esb)
(azi = -71.60766384249631, baz = 51.2691499649614, dist = 5.581417932416286e6, angle = 50.20783561306145)
```

Distances are typically accurate to the nanometre.

---

    (::GreatCircle)(a, b, zone, isnorth)
    (::GreatCircle)(a, zone_a, isnorth_a, b, zone_b, isnorth_b)

If one of `a` or `b` are `UTM`s, then provide the UTM `zone` integer and whether
the zone is in the northern hemisphere (`isnorth=true`) or not.  If both `a` and
`b` are `UTM`s, then both will be taken to be in the same zone and hemisphere.

If both `a` and `b` are `UTM`s but are not in the same zone and hemisphere, provide
the zone numbers `zone_a` and `zone_b`, respectively for points `a` and `b`, and
likewise `isnorth_a` and `isnorth_b` for the hemisphere.

---

    (::GreatCircle)(a, b, origin)
    (::GreatCircle)(a, origin_a, b, origin_b)

If one or both of `a` and `b` are `ENU`s, then provide the `origin` of the local
coordinate system.  If both `a` and `b` are in `ENU` format but are related to
`LLA` with different origins, provide `origin_a` and `origin_b` respectively for
points `a` and `b`.

---

    (::GreatCircle)(utm, zone, isnorth, enu, origin)
    (::GreatCircle)(enu, origin, utm, zone, isnorth)

If one of the points is a `UTM` and the other an `ENU`, supply the `zone`
and `isnorth` of the former, and the `origin` of the latter.
"""
(g::GreatCircle)(a::Union{LatLon,LLA}, b::Union{LatLon,LLA}) =
    GeographicLib.inverse(g.geod, a.lon, a.lat, b.lon, b.lat)

# UTM conversions
(g::GreatCircle)(a::UTM, b, zone::Integer, isnorth::Bool) =
    g(LLA(a, zone, isnorth, g.datum), b)
(g::GreatCircle)(a, b::UTM, zone::Integer, isnorth::Bool) =
    g(a, LLA(b, zone, isnorth, g.datum))
(g::GreatCircle)(a::UTM, b::UTM, zone::Integer, isnorth::Bool) =
    g(a, zone, isnorth, b, zone, isnorth)
(g::GreatCircle)(a::UTM, zone_a::Integer, isnorth_a::Bool,
                 b::UTM, zone_b::Integer, isnorth_b::Bool) =
    g(LLA(a, zone_a, isnorth_a, g.datum), LLA(b, zone_b, isnorth_b, g.datum))

# ENU conversions
(g::GreatCircle)(a::ENU, origin_a, b::ENU, origin_b) =
    g(LLA(a, origin_a, g.datum), LLA(b, origin_b, g.datum))
(g::GreatCircle)(a::ENU, b::ENU, origin) = g(a, origin, b, origin)
(g::GreatCircle)(a::ENU, b, origin) = g(LLA(a, origin, g.datum), b)
(g::GreatCircle)(a, b::ENU, origin) = g(a, LLA(b, origin, g.datum))

# ENU and UTM conversions combined
(g::GreatCircle)(a::ENU, origin, b::UTM, zone::Integer, isnorth::Bool) =
    g(LLA(a, origin, g.datum), LLA(b, zone, isnorth, g.datum))
(g::GreatCircle)(a::UTM, zone::Integer, isnorth::Bool, b::ENU, origin) =
    g(LLA(a, zone, isnorth, g.datum), LLA(b, origin, g.datum))

# Other conversions
(g::GreatCircle)(a::Union{LatLon,LLA}, b) = g(a, LLA(b, g.datum))
(g::GreatCircle)(a, b::Union{LatLon,LLA}) = g(LLA(a, g.datum), b)
(g::GreatCircle)(a, b) = g(LLA(a, g.datum), LLA(b, g.datum))

"""
    (::GreatCircle)(p; azi, dist, angle) -> (lon, lat, baz, dist, angle)

Compute the endpoint when travelling either `dist` m or an equivalent spherical angular
distance `angle`° along an azimuth `azi`° from starting point `p`.

One and only one of `dist` and `angle` should be provided.  `azi` must always be given.

The function returns a named tuple giving the longitude (`lon`, °), latitude (`lat`, °),
backazimuth from the end point (`baz`, °), distance (`dist`, m) and spherical equivalent
angular distance (`angle`, °).

#### Example

Calculate the end point from travelling northeast from Nelson's Column, London for
1 km:

```
julia> using Geodesy

julia> gc = GreatCircle(wgs84);

julia> nelson = LatLon(51.5077, -0.1277);

juila> gc(nelson, azi=45, dist=1000)
(lon = -0.11751395294238295, lat = 51.51405511443968, baz = -134.99202711278312, dist = 1000, angle = 0.008994870339334382)
```

---

    (::GreatCircle)(p, zone, isnorth; azi, dist, angle)

If `p` is a `UTM`, then provide the UTM `zone` number and whether the zone is in the
northern hemisphere (`isnorth=true`) or not.

---

    (::GreatCircle)(p, origin; azi, dist, angle)

If `p` is a `ENU`, then provide the `origin` of the local ENU coordinate system.
"""
function (g::GreatCircle)(p::Union{LatLon,LLA}; azi, dist=nothing, angle=nothing)
    sum(x -> x === nothing, (dist, angle)) == 1 ||
        throw(ArgumentError("must specify one and only one of dist or angle"))
    if dist !== nothing
        GeographicLib.forward(g.geod, p.lon, p.lat, azi, dist)
    elseif angle !== nothing
        GeographicLib.forward_deg(g.geod, p.lon, p.lat, azi, angle)
    end
end

# Conversions
(g::GreatCircle)(p::UTM, zone::Integer, isnorth::Bool; kwargs...) =
    g(LLA(p, zone, isnorth, g.datum); kwargs...)
(g::GreatCircle)(p::ENU, origin; kwargs...) = g(LLA(p, origin, g.datum); kwargs...)
(g::GreatCircle)(p; kwargs...) = g(LLA(p, g.datum); kwargs...)
