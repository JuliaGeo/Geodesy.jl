# Code adapted from MATLAB port of geograhiclib, which was licensed under
# a permissive MIT-like license and the following copyright notice
#     Copyright (c) 2016, Charles Karney
#     All rights reserved.

"""
eatanhe(x, e2) returns e*atanh(e*x) where e = sqrt(e2)
e2 is a scalar; x can be any shape.
"""
function eatanhe(x, e2)
    e = sqrt(abs(e2))
    if (e2 >= 0)
        y = e * atanh(e * x)
    else
        y = -e * atan(e * x)
    end
end

"""
Evaluate polynomial with coefficients p[i] at x
"""
function polyval(p, x)
    y = 1
    out = p[end]
    for i = length(p)-1:-1:1
        y *= x
        out += p[i] * y
    end
    return out
end

"""
a1m1(epsi) evaluates A_1 - 1 using Eq. (17).
"""
function A1m1f(epsi)
    coeff = (1, 4, 64, 0, 256)

    eps2 = epsi.^2
    t = polyval(coeff[1 : end - 1], eps2) / coeff[end]
    return (t + epsi) ./ (1 - epsi)
end

function alpf(n)
    alpcoeff = (31564, -66675, 34440, 47250, -100800, 75600, 151200,
        -1983433, 863232, 748608, -1161216, 524160, 1935360,
        670412, 406647, -533952, 184464, 725760,
        6601661, -7732800, 2230245, 7257600,
        -13675556, 3438171, 7983360,
        212378941, 319334400)

    maxpow = 6
    alp = zeros(maxpow)
    o = 1
    d = n
    for l = 1:maxpow
        m = maxpow - l
        alp[l] = d * polyval(alpcoeff[o:o+m], n) / alpcoeff[o+m+1]
        o = o + m + 2
        d = d * n
    end
    return alp
end

"""
LatFix(x) returns x is it is in the range [-90, 90]; otherwise it
returns NaN.
"""
function LatFix(x::Number)
    y = x
    if abs(x) > 90
        y = NaN
    end
    return y
end

"""
AngDiff  Compute angle difference accurately

   (d, e) = AngDiff(x, y) computes z = y - x, reduced to (-180,180].  d =
   round(z) and e = z - round(z).
"""
function AngDiff(x, y)
    (d, t) = sumx(AngNormalize(x), AngNormalize(-y))
    d = - AngNormalize(d)
    if d == 180 && t < 0
        d = -180
    end
    return sumx(d, -t)
end

"""
Error free sum

   (s, t) = sumx(u, v) returns the rounded sum u + v in s and the error in
   t, such that s + t = u + v, exactly.
"""
function sumx(u, v)
  s = u + v;
  up = s - v;
  vpp = s - up;
  up = up - u;
  vpp = vpp -  v;
  t = -(up + vpp);

  return (s,t)
end

"""
%AngNormalize:  Reduce angle to range [-180, 180)
%
%   x = AngNormalize(x) reduces angles to the range [-180, 180).
"""
function AngNormalize(x::Number)
    x = rem(x, 360)
    if x >= 180
        x = x - 360
    elseif x < -180
        x = x + 360
    end

    return x
end

"""
Simultaneous sin and cos (in degrees)
"""
function sincosdx(x)
    sind(x), cosd(x)
end

"""
z = atan2dx(y, x) compute atan2(y, x) with result in degrees in
[-180,180) and quadrant symmetries enforced.
"""
function atan2dx(y, x)
    if y > x
        tmp = x
        x = y
        y = tmp

        if x < 0 # q = 3
            x = -x

            z = atan2(y, x) * (180 / pi);       # z in [-45, 45]
            return -90 + z
        else # q = 2
            z = atan2(y, x) * (180 / pi);       # z in [-45, 45]
            return 90 - z
        end
    else
        if x < 0 # q = 1
            x = -x

            z = atan2(y, x) * (180 / pi);       # z in [-45, 45]
            if y > 0
                return 180 - z
            else
                return -180 - z
            end
        else # q = 0
            z = atan2(y, x) * (180 / pi);       # z in [-45, 45]
            return z
        end
    end
end

"""
taupf(tau, e2) returns tangent of chi in terms of tau the tangent of
phi.  e2, the square of the eccentricity, is a scalar; tau can be any
shape.
"""
function taupf(tau, e2)
  tau1 = hypot(1, tau);
  sig = sinh( eatanhe( tau ./ tau1, e2 ) );
  return hypot(1, sig) .* tau - sig .* tau1;
end

"""
(x, y) = norm2(x, y) normalize x and y so that x^2 + y^2 = 1.
"""
function norm2(x, y)
  r = hypot(x, y)
  return (x/r, y/r)
end

"""
SinCosSeries:  Evaluate a sine or cosine series using Clenshaw summation

    y = SinCosSeries(sinp, sinx, cosx, c) evaluate
    y = sum(c[i] * sin( 2*i    * x), i, 1, n), if  sinp
    y = sum(c[i] * cos((2*i-1) * x), i, 1, n), if ~sinp

where n is the size of c.  x is given via its sine and cosine in sinx
and cosx.  sinp is a scalar.  sinx, cosx, and y are scalars.  c is
a vector or tuple of length n.
"""
function SinCosSeries(sinp::Bool, sinx::Number, cosx::Number, c)
    n = length(c)
    ar = 2 * (cosx - sinx) * (cosx + sinx)
    y1 = 0.0
    if mod(n, 2) == 1
        y0 = c[n]
        n = n - 1
    else
        y0 = y1
    end

    for k = n : -2 : 1 # n:-2:2 ??
        y1 = ar * y0 - y1 + c[k]
        y0 = ar * y1 - y0 + c[k-1]
    end
    if sinp
        y = 2 * sinx * cosx * y0
    else
        y = cosx * (y0 - y1)
    end
    return y
end

"""
C1f:  Evaluate C_{1,k}

C1 = C1f(epsi) evaluates C_{1,l} using Eq. (18).
"""
function C1f(epsi)
    nC1 = 6
    coeff = (-1, 6, -16, 32,
        -9, 64, -128, 2048,
        9, -16, 768,
        3, -5, 512,
        -7, 1280,
        -7, 2048)
    C1 = zeros(nC1)
    eps2 = epsi^2
    d = epsi
    o = 1
    for l = 1:nC1
        m = div(nC1 - l, 2)
        C1[l] = d .* polyval(coeff[o : o + m], eps2) / coeff[o + m + 1]
        o = o + m + 2
        d = d * epsi
    end
    return C1
end

function transverse_mercator(lat0, lon0, lat, lon, ellipsoid::Ellipsoid) #, scale::Union{Type{Val{true}}, Type{Val{false}}} = Val{false})
    Z = 0.0
    maxpow = 6

    a = ellipsoid.a
    e2 = ellipsoid.e²
    f = e2 / (1.0 + sqrt(1.0 - e2))
    e2m = 1.0 - e2
    cc = sqrt(e2m) * exp(eatanhe(1.0, e2))
    n = f / (2.0 - f)
    alp = alpf(n)
    b1 = (1.0 - f) * (A1m1f(n) + 1.0)
    a1 = b1 * a

    lat = LatFix(lat)
    lon = AngDiff(lon0, lon)[1]

    latsign = 1 - 2 * (lat < 0)
    lonsign = 1 - 2 * (lon < 0)
    lon = lon .* lonsign
    lat = lat .* latsign
    backside = lon > 90
    if backside
        if lat == 0
           latsign = -1
        end
        lon = 180.0 - lon
    end

    (sphi, cphi) = sincosdx(lat)
    (slam, clam) = sincosdx(lon)
    tau = sphi ./ max(sqrt(realmin(Float64)), cphi)
    taup = taupf(tau, e2)

    if lat == 90
        xip = pi/2
        etap = 0.0
        gam = lon
        k = cc
    else
        xip = atan2(taup, clam)
        etap = asinh(slam ./ hypot(taup, clam))
        gam = atan2dx(slam .* taup, clam .* hypot(1, taup))
        k = sqrt(e2m + e2 * cphi.^2) .* hypot(1, tau) ./ hypot(taup, clam)
    end

    c0 = cos(2 * xip)
    ch0 = cosh(2 * etap)
    s0 = sin(2 * xip)
    sh0 = sinh(2 * etap)
    ar = 2 * c0 .* ch0
    ai = -2 * s0 .* sh0
    j = maxpow
    xi0 = Z
    yr0 = Z
    if mod(j, 2) == 1
        xi0 = xi0 + alp[j]
        yr0 = yr0 + 2 * maxpow * alp[j]
        j = j - 1;
    end
    xi1 = Z
    eta0 = Z
    eta1 = Z
    yi0 = Z
    yr1 = Z
    yi1 = Z
    for j = j : -2 : 1
        xi1  = ar .* xi0 - ai .* eta0 - xi1 + alp[j]
        eta1 = ai .* xi0 + ar .* eta0 - eta1
        yr1 = ar .* yr0 - ai .* yi0 - yr1 + 2 * j * alp[j]
        yi1 = ai .* yr0 + ar .* yi0 - yi1
        xi0  = ar .* xi1 - ai .* eta1 - xi0 + alp[j-1]
        eta0 = ai .* xi1 + ar .* eta1 - eta0
        yr0 = ar .* yr1 - ai .* yi1 - yr0 + 2 * (j-1) * alp[j-1]
        yi0 = ai .* yr1 + ar .* yi1 - yi0
    end
    ar = ar/2
    ai = ai/2
    yr1 = 1 - yr1 + ar .* yr0 - ai .* yi0
    yi1 =   - yi1 + ai .* yr0 + ar .* yi0
    ar = s0 .* ch0; ai = c0 .* sh0
    xi  = xip  + ar .* xi0 - ai .* eta0
    eta = etap + ai .* xi0 + ar .* eta0
    gam = gam - atan2dx(yi1, yr1)
    k = k .* (b1 * hypot(yr1, yi1))
    if backside
        xi = pi - xi
    end
    y = a1 * xi .* latsign
    x = a1 * eta .* lonsign
    if backside
        gam = 180 - gam
    end
    gam = AngNormalize(gam .* latsign .* lonsign)

    if lat0 != 0
        (sbet0, cbet0) = sincosdx(LatFix(lat0))
        (sbet0, cbet0) = norm2((1-f) * sbet0, cbet0)
        y0 = a1 * (atan2(sbet0, cbet0) +
                   SinCosSeries(true, sbet0, cbet0, C1f(n)));
        y = y - y0
    end

    return (x,y,gam,k)
end

# Helper for datum->ellipsoid conversion
transverse_mercator(ref::LatLon, latlon::LatLon, datum) = transverse_mercator(ref, latlon, ellipsoid(datum))

function betf(n)
     betcoeff = (384796, -382725, -6720, 932400, -1612800, 1209600, 2419200,
        -1118711, 1695744, -1174656, 258048, 80640, 3870720,
        22276, -16929, -15984, 12852, 362880,
        -830251, -158400, 197865, 7257600,
        -435388, 453717, 15966720,
        20648693, 638668800)
    maxpow = 6
    bet = zeros(maxpow)
    o = 1
    d = n
    for l = 1 : maxpow
        m = maxpow - l
        bet[l] = d * polyval(betcoeff[o : o + m], n) / betcoeff[o + m + 1]
        o = o + m + 2
        d = d * n
    end
    return bet
end

"""
%TAUF   tan(phi)
%
%   TAUF(taup, e2) returns tangent of phi in terms of taup the tangent of
%   chi.  e2, the square of the eccentricity, is a scalar; taup can be any
%   shape.
"""
function tauf(taup, e2)
    numit = 5
    e2m = 1 - e2
    tau = taup / e2m
    stol = 0.1 * sqrt(eps(Float64)) * max(1.0, abs(taup))
    g = isfinite(tau)
    for i = 1 : numit
        if !g
            break
        end
        tau1 = hypot(1, tau)
        sig = sinh( eatanhe( tau / tau1, e2 ) )
        taupa = hypot(1, sig) * tau - sig * tau1
        dtau = (taup - taupa) * (1 + e2m * tau^2) / (e2m * tau1 * hypot(1, taupa))
        tau = tau + dtau
        g = g & (abs(dtau) >= stol)
    end
    return tau
end


"""
%TRANMERC_INV  Inverse transverse Mercator projection
%
%   [lat, lon] = TRANMERC_INV(lat0, lon0, x, y)
%   [lat, lon, gam, k] = TRANMERC_INV(lat0, lon0, x, y, ellipsoid)
%
%   performs the inverse transverse Mercator projection of points (x,y) to
%   (lat,lon) using (lat0,lon0) as the center of projection.  These input
%   arguments can be scalars or arrays of equal size.  The ellipsoid vector
%   is of the form [a, e], where a is the equatorial radius in meters, e is
%   the eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by defaultellipsoid) is used.  The common
%   case of lat0 = 0 is treated efficiently provided that lat0 is specified
%   as a scalar.  projdoc defines the projection and gives the restrictions
%   on the allowed ranges of the arguments.  The forward projection is
%   given by tranmerc_fwd.  The scale on the central meridian is 1.
%
%   gam and K give metric properties of the projection at (lat,lon); gam is
%   the meridian convergence at the point and k is the scale.
%
%   lat0, lon0, lat, lon, gam are in degrees.  The projected coordinates x,
%   y are in meters (more precisely the units used for the equatorial
%   radius).  k is dimensionless.
%
%   This implementation of the projection is based on the series method
%   described in
%
%     C. F. F. Karney, Transverse Mercator with an accuracy of a few
%     nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011);
%     Addenda: http://geographiclib.sourceforge.net/tm-addenda.html
%
%   This extends the series given by Krueger (1912) to sixth order in the
%   flattening.  This is a substantially better series than that used by
%   the MATLAB mapping toolbox.  In particular the errors in the projection
%   are less than 5 nanometers withing 3900 km of the central meridian (and
%   less than 1 mm within 7600 km of the central meridian).  The mapping
%   can be continued accurately over the poles to the opposite meridian.
%
%   See also PROJDOC, TRANMERC_FWD, UTMUPS_FWD, UTMUPS_INV,
%     DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2012-2015) <charles@karney.com>.
"""
# output = [lat, lon, gam, k]
function transverse_mercator_inv(lat0, lon0, x, y, ellipsoid)
    Z = 0.0
    maxpow = 6

    a = ellipsoid.a
    e2 = ellipsoid.e²
    f = e2 / (1 + sqrt(1 - e2))
    e2m = 1 - e2
    cc = sqrt(e2m) * exp(eatanhe(1, e2))
    n = f / (2 -f)
    bet = betf(n)
    b1 = (1 - f) * (A1m1f(n) + 1)
    a1 = b1 * a

    if lat0 == 0
        y0 = 0
    else
        (sbet0, cbet0) = sincosdx(LatFix(lat0))
        (sbet0, cbet0) = norm2((1-f) * sbet0, cbet0)
        y0 = a1 * (atan2(sbet0, cbet0) + SinCosSeries(true, sbet0, cbet0, C1f(n)))
    end

    y = y + y0

    xi = y / a1
    eta = x / a1
    xisign = 1 - 2 * (xi < 0 )
    etasign = 1 - 2 * (eta < 0 )
    xi = xi * xisign
    eta = eta * etasign
    backside = xi > pi/2
    if backside
        xi = pi - xi
    end

    c0 = cos(2 * xi)
    ch0 = cosh(2 * eta)
    s0 = sin(2 * xi)
    sh0 = sinh(2 * eta)
    ar = 2 * c0 .* ch0
    ai = -2 * s0 .* sh0

    j = maxpow
    xip0 = Z
    yr0 = Z
    if mod(j, 2) == 1
        xip0 = xip0 + bet[j]
        yr0 = yr0 - 2 * maxpow * bet[j]
        j = j - 1
    end

    xip1 = Z
    etap0 = Z
    etap1 = Z
    yi0 = Z
    yr1 = Z
    yi1 = Z
    for j = j : -2 : 1
        xip1  = ar * xip0 - ai * etap0 - xip1 - bet[j]
        etap1 = ai * xip0 + ar * etap0 - etap1
        yr1 = ar * yr0 - ai * yi0 - yr1 - 2 * j * bet[j]
        yi1 = ai * yr0 + ar * yi0 - yi1
        xip0  = ar * xip1 - ai * etap1 - xip0 - bet[j-1]
        etap0 = ai * xip1 + ar * etap1 - etap0
        yr0 = ar * yr1 - ai * yi1 - yr0 - 2 * (j-1) * bet[j-1]
        yi0 = ai * yr1 + ar * yi1 - yi0
    end

    ar = ar/2
    ai = ai/2
    yr1 = 1 - yr1 + ar * yr0 - ai * yi0
    yi1 =   - yi1 + ai * yr0 + ar * yi0
    ar = s0 * ch0
    ai = c0 * sh0
    xip  = xi  + ar * xip0 - ai * etap0
    etap = eta + ai * xip0 + ar * etap0
    gam = atan2dx(yi1, yr1)
    k = b1 / hypot(yr1, yi1)
    s = sinh(etap)
    c = max(0.0, cos(xip))
    r = hypot(s, c)
    lon = atan2dx(s, c)
    sxip = sin(xip)
    tau = tauf(sxip/r, e2)
    lat = atan2dx(tau, 1)
    gam = gam + atan2dx(sxip .* tanh(etap), c)

    if r != 0
        k = k * sqrt(e2m + e2 / (1 + tau^2)) * hypot(1, tau) * r
    else
        lat = 90.0
        lon = 0.0
        k = k * cc
    end

    lat = lat * xisign
    if backside
        lon = 180 - lon
    end
    lon = lon * etasign
    lon = AngNormalize(lon + AngNormalize(lon0))
    if backside == true
        gam = 180 - gam
    end
    gam = AngNormalize(gam * xisign * etasign)

    return (lat, lon, gam, k)
end






#=
function [x, y, gam, k] = tranmerc_fwd(lat0, lon0, lat, lon, ellipsoid)
%TRANMERC_FWD  Forward transverse Mercator projection
%
%   [x, y] = TRANMERC_FWD(lat0, lon0, lat, lon)
%   [x, y, gam, k] = TRANMERC_FWD(lat0, lon0, lat, lon, ellipsoid)
%
%   performs the forward transverse Mercator projection of points (lat,lon)
%   to (x,y) using (lat0,lon0) as the center of projection.  These input
%   arguments can be scalars or arrays of equal size.  The ellipsoid vector
%   is of the form [a, e], where a is the equatorial radius in meters, e is
%   the eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by defaultellipsoid) is used.  The common
%   case of lat0 = 0 is treated efficiently provided that lat0 is specified
%   as a scalar.  projdoc defines the projection and gives the restrictions
%   on the allowed ranges of the arguments.  The inverse projection is
%   given by tranmerc_inv.  The scale on the central meridian is 1.
%
%   gam and k give metric properties of the projection at (lat,lon); gam is
%   the meridian convergence at the point and k is the scale.
%
%   lat0, lon0, lat, lon, gam are in degrees.  The projected coordinates x,
%   y are in meters (more precisely the units used for the equatorial
%   radius).  k is dimensionless.
%
%   This implementation of the projection is based on the series method
%   described in
%
%     C. F. F. Karney, Transverse Mercator with an accuracy of a few
%     nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011);
%     Addenda: http://geographiclib.sourceforge.net/tm-addenda.html
%
%   This extends the series given by Krueger (1912) to sixth order in the
%   flattening.  This is a substantially better series than that used by
%   the MATLAB mapping toolbox.  In particular the errors in the projection
%   are less than 5 nanometers withing 3900 km of the central meridian (and
%   less than 1 mm within 7600 km of the central meridian).  The mapping
%   can be continued accurately over the poles to the opposite meridian.
%
%   See also PROJDOC, TRANMERC_INV, UTMUPS_FWD, UTMUPS_INV,
%     DEFAULTELLIPSOID.

% Copyright (c) Charles Karney (2012-2015) <charles@karney.com>.

  narginchk(4, 5)
  if nargin < 5, ellipsoid = defaultellipsoid; end
  try
    S = size(lat0 + lon0 + lat + lon);
  catch
    error('lat0, lon0, lat, lon have incompatible sizes')
  end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end

  Z = zeros(prod(S),1);
  maxpow = 6;(1)

  a = ellipsoid(1);
  e2 = real(ellipsoid(2)^2);
  f = e2 / (1 + sqrt(1 - e2));
  e2m = 1 - e2;
  cc = sqrt(e2m) * exp(eatanhe(1, e2));
  n = f / (2 -f);
  alp = alpf(n);
  b1 = (1 - f) * (A1m1f(n) + 1);
  a1 = b1 * a;

  lat = LatFix(lat(:)) + Z;
  lon = AngDiff(lon0(:), lon(:)) + Z;

  latsign = 1 - 2 * (lat < 0);
  lonsign = 1 - 2 * (lon < 0);
  lon = lon .* lonsign;
  lat = lat .* latsign;
  backside = lon > 90;
  latsign(backside & lat == 0) = -1;
  lon(backside) = 180 - lon(backside);
  [sphi, cphi] = sincosdx(lat);
  [slam, clam] = sincosdx(lon);
  tau = sphi ./ max(sqrt(realmin), cphi);
  taup = taupf(tau, e2);
  xip = atan2(taup, clam);
  etap = asinh(slam ./ hypot(taup, clam));
  gam = atan2dx(slam .* taup, clam .* hypot(1, taup));
  k = sqrt(e2m + e2 * cphi.^2) .* hypot(1, tau) ./ hypot(taup, clam);
  c = ~(lat ~= 90);
  if any(c)
    xip(c) = pi/2;
    etap(c) = 0;
    gam(c) = lon(c);
    k(c) = cc;
  end
  c0 = cos(2 * xip); ch0 = cosh(2 * etap);
  s0 = sin(2 * xip); sh0 = sinh(2 * etap);
  ar = 2 * c0 .* ch0; ai = -2 * s0 .* sh0;
  j = maxpow;
  xi0 = Z; yr0 = Z;
  if mod(j, 2)
    xi0 = xi0 + alp(j);
    yr0 = yr0 + 2 * maxpow * alp(j);
    j = j - 1;
  end
  xi1 = Z; eta0 = Z; eta1 = Z;
  yi0 = Z; yr1 = Z; yi1 = Z;
  for j = j : -2 : 1
    xi1  = ar .* xi0 - ai .* eta0 - xi1 + alp(j);
    eta1 = ai .* xi0 + ar .* eta0 - eta1;
    yr1 = ar .* yr0 - ai .* yi0 - yr1 + 2 * j * alp(j);
    yi1 = ai .* yr0 + ar .* yi0 - yi1;
    xi0  = ar .* xi1 - ai .* eta1 - xi0 + alp(j-1);
    eta0 = ai .* xi1 + ar .* eta1 - eta0;
    yr0 = ar .* yr1 - ai .* yi1 - yr0 + 2 * (j-1) * alp(j-1);
    yi0 = ai .* yr1 + ar .* yi1 - yi0;
  end
  ar = ar/2; ai = ai/2;
  yr1 = 1 - yr1 + ar .* yr0 - ai .* yi0;
  yi1 =   - yi1 + ai .* yr0 + ar .* yi0;
  ar = s0 .* ch0; ai = c0 .* sh0;
  xi  = xip  + ar .* xi0 - ai .* eta0;
  eta = etap + ai .* xi0 + ar .* eta0;
  gam = gam - atan2dx(yi1, yr1);
  k = k .* (b1 * hypot(yr1, yi1));
  xi(backside) = pi - xi(backside);
  y = a1 * xi .* latsign;
  x = a1 * eta .* lonsign;
  gam(backside) = 180 - gam(backside);
  gam = AngNormalize(gam .* latsign .* lonsign);

  if isscalar(lat0) && lat0 == 0
    y0 = 0;
  else
    [sbet0, cbet0] = sincosdx(LatFix(lat0(:)));
    [sbet0, cbet0] = norm2((1-f) * sbet0, cbet0);
    y0 = a1 * (atan2(sbet0, cbet0) + ...
               SinCosSeries(true, sbet0, cbet0, C1f(n)));
  end
  y = y - y0;

  x = reshape(x, S); y = reshape(y, S);
  gam = reshape(gam, S); k = reshape(k, S);
end

function alp = alpf(n)
  persistent alpcoeff
  if isempty(alpcoeff)
    alpcoeff = [ ...
        31564, -66675, 34440, 47250, -100800, 75600, 151200, ...
        -1983433, 863232, 748608, -1161216, 524160, 1935360, ...
        670412, 406647, -533952, 184464, 725760, ...
        6601661, -7732800, 2230245, 7257600, ...
        -13675556, 3438171, 7983360, ...
        212378941, 319334400, ...
               ];
  end
  maxpow = 6;
  alp = zeros(1, maxpow);
  o = 1;
  d = n;
  for l = 1 : maxpow
    m = maxpow - l;
    alp(l) = d * polyval(alpcoeff(o : o + m), n) / alpcoeff(o + m + 1);
    o = o + m + 2;
    d = d * n;
  end
end
=#
