immutable PolarStereographic
    a::Float64
    f::Float64
    e2::Float64
    es::Float64
    e2m::Float64
    c::Float64
end

function PolarStereographic(a, f)
    e2 = f * (2-f)
    es = sign(f) * sqrt(abs(e2))
    e2m = 1 - e2
    c = (1 - f) * exp(eatanhe(1.0, es))

    if !(isfinite(a) && a > 0)
        error("Major radius is not positive")
    end
    if !(isfinite(f) && f < 1)
        error("Minor radius is not positive")
    end

    return PolarStereographic(a, f, e2, es, e2m, c)
end
function PolarStereographic(datum)
    el = ellipsoid(datum)
    f = 1 - el.b/el.a
    PolarStereographic(el.a, f)
end

function polarst_fwd(northp::Bool, k0::Float64, ps::PolarStereographic, lat, lon) # k0 is scale factor...
    lat = LatFix(lat) * (northp ? 1 : -1)

    tau = tand(lat)
    secphi = hypot(1.0, tau)
    taup = taupf(tau, ps.es)
    rho = hypot(1.0, taup) + abs(taup)
    rho = (taup >= 0 ? (lat != 90 ? 1/rho : 0.0) : rho)
    rho *= 2 * k0 * ps.a / ps.c
    k = (lat != 90 ? (rho / ps.a) * secphi * sqrt(ps.e2m + ps.e2 / (secphi*secphi)) : k0)
    x = sind(lon)
    y = cosd(lon)
    x *= rho
    y *= (northp ? -rho : rho)
    gamma = AngNormalize(northp ? lon : -lon)

    return (x, y, gamma, k)
end

function polarst_inv(northp::Bool, k0::Float64, ps::PolarStereographic, x, y) # k0 is scale factor...
    rho = hypot(x, y)
    t = (rho > 0 ? rho / (2 * k0 * ps.a / ps.c) : eps(Float64)*eps(Float64))
    taup = (1 / t - t) / 2
    tau = tauf(taup, ps.es)
    secphi = hypot(1.0, tau)
    k = (rho > 0 ? (rho / ps.a) * secphi * sqrt(ps.e2m + ps.e2 / (secphi*secphi)) : k0)
    lat = (northp ? 1 : -1) * atand(tau)
    lon = atan2(x, northp ? -y : y ) * 180/pi
    gamma = AngNormalize(northp ? lon : -lon)

    return (lat, lon, gamma, k)
end
