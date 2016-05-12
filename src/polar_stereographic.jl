# Ported to Julia by Andy Ferris, 2016, and re-released under an MIT license.
#/**
# * Copyright (c) Charles Karney (2008-2015) <charles@karney.com> and licensed
# * under the MIT/X11 License.  For more information, see
# * http://geographiclib.sourceforge.net/
# **********************************************************************/

"""
    (x, y, gamma, k) = polarst_fwd(northpole::Bool, k0::Float64, tm::TransverseMercator, lat, lon)

Perform polar-stereographic projection of `lat` and `lon` with respect to north
or south pole `northpole` and horizontal scaling `k0` (`= 0.994` for UPS).
`γ` and `k` are the local convergence and scaling factors, respectively.
"""
function polarst_fwd(northp::Bool, k0::Float64, tm::TransverseMercator, lat, lon) # k0 is scale factor...
    lat = LatFix(lat) * (northp ? 1 : -1)

    tau = tand(lat)
    secphi = hypot(1.0, tau)
    taup = taupf(tau, tm.e2) # TODO revert to C++ es here?
    rho = hypot(1.0, taup) + abs(taup)
    rho = (taup >= 0 ? (lat != 90 ? 1/rho : 0.0) : rho)
    rho *= 2 * k0 * tm.a / tm.c
    k = (lat != 90 ? (rho / tm.a) * secphi * sqrt(tm.e2m + tm.e2 / (secphi*secphi)) : k0)
    x = sind(lon)
    y = cosd(lon)
    x *= rho
    y *= (northp ? -rho : rho)
    gamma = AngNormalize(northp ? lon : -lon)

    return (x, y, gamma, k)
end

"""
    (lat, lon, gamma, k) = polarst_inv(northp::Bool, k0::Float64, tm::TransverseMercator, x, y)

Invert polar-stereographic projection of `x` and `y` with respect to north
or south pole `northpole` and horizontal scaling `k0` (`= 0.994` for UPS).
`γ` and `k` are the local convergence and scaling factors, respectively.
"""
function polarst_inv(northp::Bool, k0::Float64, tm::TransverseMercator, x, y) # k0 is scale factor...
    rho = hypot(x, y)
    t = (rho > 0 ? rho / (2 * k0 * tm.a / tm.c) : eps(Float64)*eps(Float64))
    taup = (1 / t - t) / 2
    tau = tauf(taup, tm.e2) # TODO revert to C++ es here?
    secphi = hypot(1.0, tau)
    k = (rho > 0 ? (rho / tm.a) * secphi * sqrt(tm.e2m + tm.e2 / (secphi*secphi)) : k0)
    lat = (northp ? 1 : -1) * atand(tau)
    lon = atan2(x, northp ? -y : y ) * 180/pi
    gamma = AngNormalize(northp ? lon : -lon)

    return (lat, lon, gamma, k)
end
