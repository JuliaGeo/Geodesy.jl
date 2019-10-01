module Geodesics

using StaticArrays: @SVector

import ..Math, ..Constants, ..GeodesicCapability, ..Result

export Geodesic, ArcDirect, Direct, Inverse

const GEOGRAPHICLIB_GEODESIC_ORDER = 6
const nA1_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nC1_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nC1p_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nA2_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nC2_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nA3_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nA3x_ = nA3_
const nC3_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nC3x_ = (nC3_ * (nC3_ - 1)) ÷ 2
const nC4_ = GEOGRAPHICLIB_GEODESIC_ORDER
const nC4x_ = (nC4_ * (nC4_ + 1)) ÷ 2
const maxit1_ = 20
const maxit2_ = maxit1_ + Math.digits + 10

const tiny_ = sqrt(Math.minval)
const tol0_ = Math.epsilon
const tol1_ = 200 * tol0_
const tol2_ = sqrt(tol0_)
const tolb_ = tol0_ * tol2_
const xthresh_ = 1000 * tol2_

const CAP_NONE = GeodesicCapability.CAP_NONE
const CAP_C1   = GeodesicCapability.CAP_C1
const CAP_C1p  = GeodesicCapability.CAP_C1p
const CAP_C2   = GeodesicCapability.CAP_C2
const CAP_C3   = GeodesicCapability.CAP_C3
const CAP_C4   = GeodesicCapability.CAP_C4
const CAP_ALL  = GeodesicCapability.CAP_ALL
const CAP_MASK = GeodesicCapability.CAP_MASK
const OUT_ALL  = GeodesicCapability.OUT_ALL
const OUT_MASK = GeodesicCapability.OUT_MASK

"""No capabilities, no output."""
const EMPTY         = GeodesicCapability.EMPTY
"""Calculate latitude `lat2`."""
const LATITUDE      = GeodesicCapability.LATITUDE
"""Calculate longitude `lon2`."""
const LONGITUDE     = GeodesicCapability.LONGITUDE
"""Calculate azimuths `azi1` and `azi2`."""
const AZIMUTH       = GeodesicCapability.AZIMUTH
"""Calculate distance `s12`."""
const DISTANCE      = GeodesicCapability.DISTANCE
"""All of the above."""
const STANDARD      = GeodesicCapability.STANDARD
"""Allow distance `s12` to be used as input in the direct geodesic problem."""
const DISTANCE_IN   = GeodesicCapability.DISTANCE_IN
"""Calculate reduced length `m12`."""
const REDUCEDLENGTH = GeodesicCapability.REDUCEDLENGTH
"""Calculate geodesic scales `M12` and `M21`."""
const GEODESICSCALE = GeodesicCapability.GEODESICSCALE
"""Calculate area `S12`."""
const AREA          = GeodesicCapability.AREA
"""All of the above."""
const ALL           = GeodesicCapability.ALL
"""Unroll longitudes, rather than reducing them to the range [-180°,180°]."""
const LONG_UNROLL   = GeodesicCapability.LONG_UNROLL

struct Geodesic
    a::Float64
    f::Float64
    _f1::Float64
    _e2::Float64
    _ep2::Float64
    _n::Float64
    _b::Float64
    _c2::Float64
    _etol2::Float64
    _A3x::Vector{Float64}
    _C3x::Vector{Float64}
    _C4x::Vector{Float64}
end

# Treat the Geodesic struct as a scalar
Base.Broadcast.broadcastable(geod::Geodesic) = Ref(geod)

Base.:(==)(g1::Geodesic, g2::Geodesic) =
    all(getfield(g1, f) == getfield(g2, f) for f in fieldnames(Geodesic))

"""
    Geodesic(a, f) -> geodesic

Set up an ellipsoid for geodesic calculations.  `a` is the semimajor radius of the
ellipsoid, whilst flattening is given by `f`.
"""
function Geodesic(a, f)
    f1 = 1 - f
    e2 = f*(2 - f)
    ep2 = e2/Math.sq(f1) # e2 / (1 - e2)
    n = f/(2 - f)
    b = a*f1
    # authalic radius squared
    term = if e2 == 0
        1.0
    elseif e2 > 0
        Math.atanh(sqrt(e2))/sqrt(abs(e2))
    else
        atan(sqrt(-e2))/sqrt(abs(e2))
    end
    c2 = (Math.sq(a) + Math.sq(b) * term) / 2
    # The sig12 threshold for "really short".  Using the auxiliary sphere
    # solution with dnm computed at (bet1 + bet2) / 2, the relative error in
    # the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
    # (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given
    # f and sig12, the max error occurs for lines near the pole.  If the old
    # rule for computing dnm = (dn1 + dn2)/2 is used, then the error increases
    # by a factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.
    # Here 0.1 is a safety factor (error decreased by 100) and max(0.001,
    # abs(f)) stops etol2 getting too large in the nearly spherical case.
    etol2 = 0.1 * tol2_ / sqrt(max(0.001, abs(f)) * min(1.0, 1-f/2) / 2)
    !(Math.isfinite(a) && a > 0) &&
      throw(ArgumentError("Major radius is not positive"))
    !(Math.isfinite(b) && b > 0) &&
      throw(ArgumentError("Minor radius is not positive"))
    A3x = Vector{Float64}(undef, nA3x_)
    C3x = Vector{Float64}(undef, nC3x_)
    C4x = Vector{Float64}(undef, nC4x_)
    self = Geodesic(a, f, f1, e2, ep2, n, b, c2, etol2, A3x, C3x, C4x)
    _A3coeff(self)
    _C3coeff(self)
    _C4coeff(self)
    self
end

"""
    Inverse(geodesic, lat1, lon1, lat2, lon2, outmask=STANDARD) -> result::Result

Solve the inverse geodesic problem and return a `Result` containing the
parameters of interest.

Input arguments:
- `lat1`: latitude of the first point in degrees
- `lon1`: longitude of the first point in degrees
- `lat2`: latitude of the second point in degrees
- `lon2`: longitude of the second point in degrees
- `outmask`: a mask setting which output values are computed (see note below)

Compute geodesic between (`lat1`, `lon1`) and (`lat2`, `lon2`).
The default value of `outmask` is `Geodesics.STANDARD`, i.e., the `lat1`,
`lon1`, `azi1`, `lat2`, `lon2`, `azi2`, `s12`, `a12` entries are returned.

### Output mask

May be any combination of:
`Geodesics.EMPTY`, `Geodesics.LATITUDE`, `Geodesics.LONGITUDE`,
`Geodesics.AZIMUTH`, `Geodesics.DISTANCE`, `Geodesics.STANDARD`,
`Geodesics.DISTANCE_IN`, `Geodesics.REDUCEDLENGTH`, `Geodesics.GEODESICSCALE`,
`Geodesics.AREA`, `Geodesics.ALL` or `Geodesics.LONG_UNROLL`.
See the docstring for each for more information.

Flags are combined by bitwise or-ing values together, e.g.
`Geodesics.AZIMUTH | Geodesics.DISTANCE`.
"""
function Inverse(self::Geodesic, lat1::T1, lon1::T2, lat2::T3, lon2::T4,
            outmask = GeodesicCapability.STANDARD) where {T1,T2,T3,T4}

  T = float(promote_type(Float64, T1, T2, T3, T4))
  lat1, lon1, lat2, lon2 = T.((lat1, lon1, lat2, lon2))

  a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 = _GenInverse(self,
    lat1, lon1, lat2, lon2, outmask)
  outmask &= OUT_MASK
  if (outmask & LONG_UNROLL) > 0
    lon12, e = Math.AngDiff(lon1, lon2)
    lon2 = (lon1 + lon12) + e
  else
    lon2 = Math.AngNormalize(lon2)
  end
  result = Result()
  result.lat1 = Math.LatFix(lat1)
  result.lon1 = (outmask & LONG_UNROLL) > 0 ? lon1 : Math.AngNormalize(lon1)
  result.lat2 = Math.LatFix(lat2)
  result.lon2 = lon2
  result.a12 = a12
  if (outmask & DISTANCE) > 0
      result.s12 = s12
  end
  if (outmask & AZIMUTH) > 0
    result.azi1 = atand(salp1, calp1)
    result.azi2 = atand(salp2, calp2)
  end
  if (outmask & REDUCEDLENGTH) > 0
    result.m12 = m12
  end
  if (outmask & GEODESICSCALE) > 0
    result.M12 = M12
    result.M21 = M21
  end
  if (outmask & AREA) > 0
    result.S12 = S12
  end
  result
end

"""Private: Evaluate a trig series using Clenshaw summation."""
function _SinCosSeries(sinp, sinx::T1, cosx::T2, c) where {T1,T2}
    # Evaluate
    # y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
    #            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
    # using Clenshaw summation.  N.B. c[0] is unused for sin series
    # Approx operation count = (n + 5) mult and (2 * n + 2) add
    T = float(promote_type(T1, T2, eltype(c)))
    # T = identity
    k = length(c)               # Point to one beyond last element
    n = k - Int(sinp)
    ar = 2 * (cosx - sinx) * (cosx + sinx) # 2 * cos(2 * x)
    y1 = T(0)                                 # accumulators for sum
    if n & 1 == 1
        k -= 1
        y0 = T(c[k+1])
    else
        y0 = T(0)
    end
    # Now n is even
    n = n ÷ 2
    while n != 0                # while n--:
        n -= 1
        # Unroll loop x 2, so accumulators return to their original role
        k -= 1
        y1 = ar * y0 - y1 + c[k+1]
        k -= 1
        y0 = ar * y1 - y0 + c[k+1]
    end
    sinp ? 2 * sinx * cosx * y0 : # sin(2 * x) * y0
           cosx * (y0 - y1)       # cos(x) * (y0 - y1)
end

"""Private: solve astroid equation."""
function _Astroid(x, y)
    # Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
    # This solution is adapted from Geocentric::Reverse.
    p = Math.sq(x)
    q = Math.sq(y)
    r = (p + q - 1) / 6
    if !(q == 0 && r <= 0)
        # Avoid possible division by zero when r = 0 by multiplying equations
        # for s and t by r^3 and r, resp.
        S = p * q / 4            # S = r^3 * s
        r2 = Math.sq(r)
        r3 = r * r2
        # The discrimant of the quadratic equation for T3.  This is zero on
        # the evolute curve p^(1/3)+q^(1/3) = 1
        disc = S * (S + 2 * r3)
        u = r
        if disc >= 0
            T3 = S + r3
            # Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
            # of precision due to cancellation.  The result is unchanged because
            # of the way the T is used in definition of u.
            T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc) # T3 = (r * t)^3
            # N.B. cbrt always returns the real root.  cbrt(-8) = -2.
            T = Math.cbrt(T3)       # T = r * t
            # T can be zero; but then r2 / T -> 0.
            u += T + (T != 0 ? (r2 / T) : 0.0)
        else
            # T is complex, but the way u is defined the result is real.
            ang = atan(sqrt(-disc), -(S + r3))
            # There are three possible cube roots.  We choose the root which
            # avoids cancellation.  Note that disc < 0 implies that r < 0.
            u += 2 * r * cos(ang / 3)
        end
        v = sqrt(Math.sq(u) + q) # guaranteed positive
        # Avoid loss of accuracy when u < 0.
        uv = u < 0 ? q/(v - u) : u + v # u+v, guaranteed positive
        w = (uv - q) / (2 * v)               # positive?
        # Rearrange expression for k to avoid loss of accuracy due to
        # subtraction.  Division by 0 not possible because uv > 0, w >= 0.
        k = uv / (sqrt(uv + Math.sq(w)) + w) # guaranteed positive
    else                                        # q == 0 && r <= 0
        # y = 0 with |x| <= 1.  Handle this case directly.
        # for y small, positive root is k = abs(y)/sqrt(1-x^2)
        k = 0.0
    end
    k
end

"""Private: return A1-1."""
function _A1m1f(eps)
    coeff = @SVector[1, 4, 64, 0, 256]
    m = nA1_ ÷ 2
    t = Math.polyval(m, coeff, 0, Math.sq(eps)) / coeff[m + 2]
    (t + eps)/(1 - eps)
end

"""Private: return C1."""
function _C1f(eps, c)
    coeff = @SVector [-1, 6, -16, 32, -9, 64, -128, 2048, 9, -16, 768, 3,  -5, 512,
                      -7, 1280, -7, 2048]
    eps2 = Math.sq(eps)
    d = eps
    o = 0
    for l in 1:nC1_  # l is index of C1p[l]
        m = (nC1_ - l) ÷ 2        # order of polynomial in eps^2
        c[l+1] = d * Math.polyval(m, coeff, o, eps2) / coeff[o + m + 2]
        o += m + 2
        d *= eps
    end
    c
end

"""Private: return C1'"""
function _C1pf(eps, c)
    coeff = @SVector [205, -432, 768, 1536, 4005, -4736, 3840, 12288, -225, 116, 384,
                      -7173, 2695, 7680, 3467, 7680, 38081, 61440]
    eps2 = Math.sq(eps)
    d = eps
    o = 0
    for l in 1:nC1p_ # l is index of C1p[l]
        m = (nC1p_ - l) ÷ 2 # order of polynomial in eps^2
        c[l+1] = d * Math.polyval(m, coeff, o, eps2) / coeff[o + m + 2]
        o += m + 2
        d *= eps
    end
    c
end

# TODO: Double check that this is correct.  Differs from Python library by
#       1.0e-17 for eps = 0.01
#       8.8e-10 for eps = 0.1
"""Private: return A2-1"""
function _A2m1f(eps)
    coeff = @SVector [-11, -28, -192, 0, 256]
    m = nA2_ ÷ 2
    t = Math.polyval(m, coeff, 0, Math.sq(eps)) / coeff[m + 2]
    (t - eps) / (1 + eps)
end

"""Private: return C2"""
function _C2f(eps, c)
    coeff = @SVector [1, 2, 16, 32, 35, 64, 384, 2048, 15, 80, 768, 7, 35, 512,
             63, 1280, 77, 2048]
    eps2 = Math.sq(eps)
    d = eps
    o = 0
    for l in 1:nC2_ # l is index of C2[l]
        m = (nC2_ - l) ÷ 2        # order of polynomial in eps^2
        c[l+1] = d * Math.polyval(m, coeff, o, eps2) / coeff[o + m + 2]
        o += m + 2
        d *= eps
    end
    c
end

"""Private: return coefficients for A3"""
function _A3coeff(self::Geodesic)
    coeff = (-3, 128, -2, -3, 64, -1, -3, -1, 16, 3, -1, -2, 8, 1, -1, 2, 1, 1)
    o = k = 0
    for j in (nA3_-1):-1:0 # coeff of eps^j
        m = min(nA3_ - j - 1, j) # order of polynomial in n
        self._A3x[k+1] = Math.polyval(m, coeff, o, self._n) / coeff[o + m + 2]
        k += 1
        o += m + 2
    end
    self._A3x
end

"""Private: return coefficients for C3"""
function _C3coeff(self::Geodesic)
    coeff = [3, 128, 2, 5, 128, -1, 3, 3, 64, -1, 0, 1, 8, -1, 1, 4, 5, 256,
             1, 3, 128, -3, -2, 3, 64, 1, -3, 2, 32, 7, 512, -10, 9, 384, 5,
             -9, 5, 192, 7, 512, -14, 7, 512, 21, 2560]
    o = k = 0
    for l in 1:(nC3_ - 1) # l is index of C3[l]
        for j in (nC3_ - 1):-1:l # coeff of eps^j
            m = min(nC3_ - j - 1, j) # order of polynomial in n
            self._C3x[k+1] = Math.polyval(m, coeff, o, self._n) / coeff[o + m + 2]
            k += 1
            o += m + 2
        end
    end
    self._C3x
end

"""Private: return coefficients for C4"""
function _C4coeff(self::Geodesic)
    coeff = [97, 15015, 1088, 156, 45045, -224, -4784, 1573, 45045, -10656, 14144,
             -4576, -858, 45045, 64, 624, -4576, 6864, -3003, 15015, 100, 208, 572,
             3432, -12012, 30030, 45045, 1, 9009, -2944, 468, 135135, 5792, 1040,
             -1287, 135135, 5952, -11648, 9152, -2574, 135135, -64, -624, 4576,
             -6864, 3003, 135135, 8, 10725, 1856, -936, 225225, -8448, 4992, -1144,
             225225, -1440, 4160, -4576, 1716, 225225, -136, 63063, 1024, -208,
             105105, 3584, -3328, 1144, 315315, -128, 135135, -2560, 832, 405405,
             128, 99099]
    o = k = 0
    for l in 0:(nC4_ - 1) # l is index of C4[l]
        for j in (nC4_ - 1):-1:l # coeff of eps^j
            m = nC4_ - j - 1 # order of polynomial in n
            self._C4x[k+1] = Math.polyval(m, coeff, o, self._n) / coeff[o + m + 2]
            k += 1
            o += m + 2
        end
    end
    self._C4x
end

"""Private: return A3"""
_A3f(self::Geodesic, eps) = Math.polyval(nA3_ - 1, self._A3x, 0, eps)

"""Private: return C3"""
function _C3f(self::Geodesic, eps, c)
    # Evaluate C3
    # Elements c[1] thru c[nC3_ - 1] are set
    mult = oneunit(eps)
    o = 0
    for l in 1:(nC3_ - 1) # l is index of C3[l]
        m = nC3_ - l - 1       # order of polynomial in eps
        mult *= eps
        c[l+1] = mult * Math.polyval(m, self._C3x, o, eps)
        o += m + 1
    end
    c
end

"""Private: return C4"""
function _C4f(self::Geodesic, eps, c)
    # Evaluate C4 coeffs by Horner's method
    # Elements c[0] thru c[nC4_ - 1] are set
    mult = oneunit(eps)
    o = 0
    for l in 0:(nC4_ - 1) # l is index of C4[l]
        m = nC4_ - l - 1    # order of polynomial in eps
        c[l+1] = mult * Math.polyval(m, self._C4x, o, eps)
        o += m + 1
        mult *= eps
    end
    c
end

"""Private: return a bunch of lengths"""
function _Lengths(self::Geodesic, eps, sig12,
             ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, outmask,
             # Scratch areas of the right size
             C1a, C2a)
    # Return s12b, m12b, m0, M12, M21, where
    # m12b = (reduced length)/_b; s12b = distance/_b,
    # m0 = coefficient of secular term in expression for reduced length.
    outmask &= OUT_MASK
    # outmask & DISTANCE: set s12b
    # outmask & REDUCEDLENGTH: set m12b & m0
    # outmask & GEODESICSCALE: set M12 & M21

    s12b = m12b = m0 = M12 = M21 = Math.nan
    if (outmask & (DISTANCE | REDUCEDLENGTH | GEODESICSCALE)) > 0
        A1 = _A1m1f(eps)
        _C1f(eps, C1a)
        if (outmask & (REDUCEDLENGTH | GEODESICSCALE)) > 0
            A2 = _A2m1f(eps)
            _C2f(eps, C2a)
            m0x = A1 - A2
            A2 = 1 + A2
        end
        A1 = 1 + A1
    end
    if (outmask & DISTANCE) > 0
        B1 = _SinCosSeries(true, ssig2, csig2, C1a) - _SinCosSeries(true, ssig1, csig1, C1a)
        # Missing a factor of _b
        s12b = A1 * (sig12 + B1)
        if (outmask & (REDUCEDLENGTH | GEODESICSCALE)) > 0
            B2 = _SinCosSeries(true, ssig2, csig2, C2a) - _SinCosSeries(true, ssig1, csig1, C2a)
            J12 = m0x * sig12 + (A1 * B1 - A2 * B2)
        end
    elseif (outmask & (REDUCEDLENGTH | GEODESICSCALE)) > 0
        # Assume here that nC1_ >= nC2_
        for l in 1:(nC2_ - 1)
            C2a[l+1] = A1 * C1a[l+1] - A2 * C2a[l+1]
            J12 = m0x * sig12 + _SinCosSeries(true, ssig2, csig2, C2a) -
                            _SinCosSeries(true, ssig1, csig1, C2a)
        end
    end
    if (outmask & REDUCEDLENGTH) > 0
        m0 = m0x
        # Missing a factor of _b.
        # Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
        # accurate cancellation in the case of coincident points.
        m12b = (dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) -
           csig1 * csig2 * J12)
    end
    if (outmask & GEODESICSCALE) > 0
        csig12 = csig1 * csig2 + ssig1 * ssig2
        t = self._ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
        M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1
        M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2
    end
    return s12b, m12b, m0, M12, M21
end

"""Private: Find a starting value for Newton's method."""
function _InverseStart(self, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                       lam12, slam12, clam12,
                       # Scratch areas of the right size
                       C1a, C2a)
  # Return a starting point for Newton's method in salp1 and calp1 (function
  # value is -1).  If Newton's method doesn't need to be used, return also
  # salp2 and calp2 and function value is sig12.
  sig12 = -1.0
  salp2 = calp2 = dnm = Math.nan # Return values
  # bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
  sbet12 = sbet2 * cbet1 - cbet2 * sbet1
  cbet12 = cbet2 * cbet1 + sbet2 * sbet1
  # Volatile declaration needed to fix inverse cases
  # 88.202499451857 0 -88.202499451857 179.981022032992859592
  # 89.262080389218 0 -89.262080389218 179.992207982775375662
  # 89.333123580033 0 -89.333123580032997687 179.99295812360148422
  # which otherwise fail with g++ 4.4.4 x86 -O3
  sbet12a = sbet2 * cbet1
  sbet12a += cbet2 * sbet1

  shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5
  if shortline
    sbetm2 = Math.sq(sbet1 + sbet2)
    # sin((bet1+bet2)/2)^2
    # =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
    sbetm2 /= sbetm2 + Math.sq(cbet1 + cbet2)
    dnm = sqrt(1 + self._ep2 * sbetm2)
    omg12 = lam12 / (self._f1 * dnm)
    somg12 = sin(omg12)
    comg12 = cos(omg12)
  else
    somg12 = slam12
    comg12 = clam12
  end

  salp1 = cbet2 * somg12
  calp1 = if comg12 >= 0
        sbet12 + cbet2 * sbet1 * Math.sq(somg12) / (1 + comg12)
    else
        sbet12a - cbet2 * sbet1 * Math.sq(somg12) / (1 - comg12)
    end

  ssig12 = hypot(salp1, calp1)
  csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12

  if shortline && ssig12 < self._etol2
    # really short lines
    salp2 = cbet1 * somg12
    calp2 = sbet12 - cbet1 * sbet2 * (comg12 >= 0 ?
                                          Math.sq(somg12) / (1 + comg12) :
                                          1 - comg12)
    salp2, calp2 = Math.norm(salp2, calp2)
    # Set return value
    sig12 = atan(ssig12, csig12)
  elseif (abs(self._n) >= 0.1 || # Skip astroid calc if too eccentric
        csig12 >= 0 ||
        ssig12 >= 6 * abs(self._n) * pi * Math.sq(cbet1))
    # Nothing to do, zeroth order spherical approximation is OK
  else
    # Scale lam12 and bet2 to x, y coordinate system where antipodal point
    # is at origin and singular point is at y = 0, x = -1.
    # real y, lamscale, betscale
    # Volatile declaration needed to fix inverse case
    # 56.320923501171 0 -56.320923501171 179.664747671772880215
    # which otherwise fails with g++ 4.4.4 x86 -O3
    # volatile real x
    lam12x = atan(-slam12, -clam12)
    if self.f >= 0             # In fact f == 0 does not get here
      # x = dlong, y = dlat
      k2 = Math.sq(sbet1) * self._ep2
      eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)
      lamscale = self.f * cbet1 * _A3f(self, eps) * pi
      betscale = lamscale * cbet1
      x = lam12x / lamscale
      y = sbet12a / betscale
    else                      # _f < 0
      # x = dlat, y = dlong
      cbet12a = cbet2 * cbet1 - sbet2 * sbet1
      bet12a = atan(sbet12a, cbet12a)
      # real m12b, m0, dummy
      # In the case of lon12 = 180, this repeats a calculation made in
      # Inverse.
      dummy, m12b, m0, dummy, dummy = _Lengths(self,
        self._n, pi + bet12a, sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
        cbet1, cbet2, REDUCEDLENGTH, C1a, C2a)
      x = -1 + m12b / (cbet1 * cbet2 * m0 * pi)
      betscale = (x < -0.01 ? sbet12a / x :
                  -self.f * Math.sq(cbet1) * pi)
      lamscale = betscale / cbet1
      y = lam12x / lamscale
    end

    if y > -tol1_ && x > -1 - xthresh_
      # strip near cut
      if self.f >= 0
        salp1 = min(1.0, -x)
        calp1 = - sqrt(1 - Math.sq(salp1))
      else
        calp1 = max((x > -tol1_ ? 0.0 : -1.0), x)
        salp1 = sqrt(1 - Math.sq(calp1))
      end
    else
      # Estimate alp1, by solving the astroid problem.
      #
      # Could estimate alpha1 = theta + pi/2, directly, i.e.,
      #   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
      #   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
      #
      # However, it's better to estimate omg12 from astroid and use
      # spherical formula to compute alp1.  This reduces the mean number of
      # Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
      # (min 0 max 5).  The changes in the number of iterations are as
      # follows:
      #
      # change percent
      #    1       5
      #    0      78
      #   -1      16
      #   -2       0.6
      #   -3       0.04
      #   -4       0.002
      #
      # The histogram of iterations is (m = number of iterations estimating
      # alp1 directly, n = number of iterations estimating via omg12, total
      # number of trials = 148605):
      #
      #  iter    m      n
      #    0   148    186
      #    1 13046  13845
      #    2 93315 102225
      #    3 36189  32341
      #    4  5396      7
      #    5   455      1
      #    6    56      0
      #
      # Because omg12 is near pi, estimate work with omg12a = pi - omg12
      k = _Astroid(x, y)
      omg12a = lamscale * (self.f >= 0 ? -x * k/(1 + k) :
                            -y * (1 + k)/k )
      somg12 = sin(omg12a)
      comg12 = -cos(omg12a)
      # Update spherical estimate of alp1 using omg12 instead of lam12
      salp1 = cbet2 * somg12
      calp1 = sbet12a - cbet2 * sbet1 * Math.sq(somg12) / (1 - comg12)
    end
  end
  # Sanity check on starting guess.  Backwards check allows NaN through.
  if !(salp1 <= 0)
    salp1, calp1 = Math.norm(salp1, calp1)
  else
    salp1 = 1.0
    calp1 = 0.0
  end
  sig12, salp1, calp1, salp2, calp2, dnm
end

"""Private: Solve hybrid problem"""
function _Lambda12(self, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
              slam120, clam120, diffp,
              # Scratch areas of the right size
              C1a, C2a, C3a)
  if sbet1 == 0 && calp1 == 0
    # Break degeneracy of equatorial line.  This case has already been
    # handled.
    calp1 = -tiny_
  end

  # sin(alp1) * cos(bet1) = sin(alp0)
  salp0 = salp1 * cbet1
  calp0 = hypot(calp1, salp1 * sbet1) # calp0 > 0

  # real somg1, comg1, somg2, comg2, lam12
  # tan(bet1) = tan(sig1) * cos(alp1)
  # tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
  ssig1 = sbet1; somg1 = salp0 * sbet1
  csig1 = comg1 = calp1 * cbet1
  ssig1, csig1 = Math.norm(ssig1, csig1)
  # Math.norm(somg1, comg1); -- don't need to normalize!

  # Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
  # about this case, since this can yield singularities in the Newton
  # iteration.
  # sin(alp2) * cos(bet2) = sin(alp0)
  salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1
  # calp2 = sqrt(1 - sq(salp2))
  #       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
  # and subst for calp0 and rearrange to give (choose positive sqrt
  # to give alp2 in [0, pi/2]).
  calp2 = if cbet2 != cbet1 || abs(sbet2) != -sbet1
      sqrt(Math.sq(calp1 * cbet1) +
                     (cbet1 < -sbet1 ?
                      (cbet2 - cbet1) * (cbet1 + cbet2) :
                      (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2
      else
          abs(calp1)
      end
  # tan(bet2) = tan(sig2) * cos(alp2)
  # tan(omg2) = sin(alp0) * tan(sig2).
  ssig2 = sbet2; somg2 = salp0 * sbet2
  csig2 = comg2 = calp2 * cbet2
  ssig2, csig2 = Math.norm(ssig2, csig2)
  # Math.norm(somg2, comg2); -- don't need to normalize!

  # sig12 = sig2 - sig1, limit to [0, pi]
  sig12 = atan(max(0.0, csig1 * ssig2 - ssig1 * csig2),
               csig1 * csig2 + ssig1 * ssig2)

  # omg12 = omg2 - omg1, limit to [0, pi]
  somg12 = max(0.0, comg1 * somg2 - somg1 * comg2)
  comg12 =          comg1 * comg2 + somg1 * somg2
  # eta = omg12 - lam120
  eta = atan(somg12 * clam120 - comg12 * slam120,
             comg12 * clam120 + somg12 * slam120)

  # real B312
  k2 = Math.sq(calp0) * self._ep2
  eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)
  _C3f(self, eps, C3a)
  B312 = (_SinCosSeries(true, ssig2, csig2, C3a) -
          _SinCosSeries(true, ssig1, csig1, C3a))
  domg12 =  -self.f * _A3f(self, eps) * salp0 * (sig12 + B312)
  lam12 = eta + domg12

  if diffp
    if calp2 == 0
      dlam12 = - 2 * self._f1 * dn1 / sbet1
    else
      dummy, dlam12, dummy, dummy, dummy = _Lengths(self,
        eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
        REDUCEDLENGTH, C1a, C2a)
      dlam12 *= self._f1 / (calp2 * cbet2)
    end
  else
    dlam12 = Math.nan
  end

  lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12
end

"""Private: General version of the inverse problem"""
function _GenInverse(self, lat1::T1, lon1::T2, lat2::T3, lon2::T4, outmask) where {T1,T2,T3,T4}
  T = float(promote_type(Float64, T1, T2, T3, T4))
  lat1, lon1, lat2, lon2 = T.((lat1, lon1, lat2, lon2))
  a12 = s12 = m12 = M12 = M21 = S12 = Math.nan # return vals

  outmask &= OUT_MASK
  # Compute longitude difference (AngDiff does this carefully).  Result is
  # in [-180, 180] but -180 is only for west-going geodesics.  180 is for
  # east-going and meridional geodesics.
  lon12, lon12s = Math.AngDiff(lon1, lon2)
  # Make longitude difference positive.
  lonsign = lon12 >= 0  ? 1 : -1
  # If very close to being on the same half-meridian, then make it so.
  lon12 = lonsign * Math.AngRound(lon12)
  lon12s = Math.AngRound((180 - lon12) - lonsign * lon12s)
  lam12 = deg2rad(lon12)
  if lon12 > 90
    slam12, clam12 = Math.sincosd(lon12s)
    clam12 = -clam12
  else
    slam12, clam12 = Math.sincosd(lon12)
  end

  # If really close to the equator, treat as on equator.
  lat1 = Math.AngRound(Math.LatFix(lat1))
  lat2 = Math.AngRound(Math.LatFix(lat2))
  # Swap points so that point with higher (abs) latitude is point 1
  # If one latitude is a nan, then it becomes lat1.
  swapp = abs(lat1) < abs(lat2) ? -1 : 1
  if swapp < 0
    lonsign *= -1
    lat2, lat1 = lat1, lat2
  end
  # Make lat1 <= 0
  latsign = lat1 < 0 ? 1 : -1
  lat1 *= latsign
  lat2 *= latsign
  # Now we have
  #
  #     0 <= lon12 <= 180
  #     -90 <= lat1 <= 0
  #     lat1 <= lat2 <= -lat1
  #
  # longsign, swapp, latsign register the transformation to bring the
  # coordinates to this canonical form.  In all cases, 1 means no change was
  # made.  We make these transformations so that there are few cases to
  # check, e.g., on verifying quadrants in atan2.  In addition, this
  # enforces some symmetries in the results returned.

  # real phi, sbet1, cbet1, sbet2, cbet2, s12x, m12x

  sbet1, cbet1 = Math.sincosd(lat1)
  sbet1 *= self._f1
  # Ensure cbet1 = +epsilon at poles
  sbet1, cbet1 = Math.norm(sbet1, cbet1)
  cbet1 = max(tiny_, cbet1)

  sbet2, cbet2 = Math.sincosd(lat2); sbet2 *= self._f1
  # Ensure cbet2 = +epsilon at poles
  sbet2, cbet2 = Math.norm(sbet2, cbet2)
  cbet2 = max(tiny_, cbet2)

  # If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
  # |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
  # a better measure.  This logic is used in assigning calp2 in Lambda12.
  # Sometimes these quantities vanish and in that case we force bet2 = +/-
  # bet1 exactly.  An example where is is necessary is the inverse problem
  # 48.522876735459 0 -48.52287673545898293 179.599720456223079643
  # which failed with Visual Studio 10 (Release and Debug)

  if cbet1 < -sbet1
    if cbet2 == cbet1
      sbet2 = sbet2 < 0 ? sbet1 : -sbet1
    end
  else
    if abs(sbet2) == -sbet1
      cbet2 = cbet1
    end
  end

  dn1 = sqrt(1 + self._ep2 * Math.sq(sbet1))
  dn2 = sqrt(1 + self._ep2 * Math.sq(sbet2))

  # real a12, sig12, calp1, salp1, calp2, salp2
  # index zero elements of these arrays are unused
  C1a = Vector{Float64}(undef, nC1_ + 1)
  C2a = Vector{Float64}(undef, nC2_ + 1)
  C3a = Vector{Float64}(undef, nC3_)

  meridian = lat1 == -90 || slam12 == 0

  if meridian

    # Endpoints are on a single full meridian, so the geodesic might lie on
    # a meridian.

    calp1 = clam12
    salp1 = slam12 # Head to the target longitude
    calp2 = 1.0
    salp2 = 0.0       # At the target we're heading north

    # tan(bet) = tan(sig) * cos(alp)
    ssig1 = sbet1; csig1 = calp1 * cbet1
    ssig2 = sbet2; csig2 = calp2 * cbet2

    # sig12 = sig2 - sig1
    sig12 = atan(max(0.0, csig1 * ssig2 - ssig1 * csig2),
                          csig1 * csig2 + ssig1 * ssig2)

    s12x, m12x, dummy, M12, M21 = _Lengths(self,
      self._n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
      outmask | DISTANCE | REDUCEDLENGTH, C1a, C2a)

    # Add the check for sig12 since zero length geodesics might yield m12 <
    # 0.  Test case was
    #
    #    echo 20.001 0 20.001 0 | GeodSolve -i
    #
    # In fact, we will have sig12 > pi/2 for meridional geodesic which is
    # not a shortest path.
    if sig12 < 1 || m12x >= 0
      if sig12 < 3 * tiny_
        sig12 = m12x = s12x = 0.0
      end
      m12x *= self._b
      s12x *= self._b
      a12 = rad2deg(sig12)
    else
      # m12 < 0, i.e., prolate and too close to anti-podal
      meridian = false
    end
  end
  # end if meridian:

  # somg12 > 1 marks that it needs to be calculated
  somg12 = 2.0
  comg12 = 0.0
  omg12 = 0.0
  if (! meridian &&
      sbet1 == 0 &&   # and sbet2 == 0
      # Mimic the way Lambda12 works with calp1 = 0
      (self.f <= 0 || lon12s >= self.f * 180))

    # Geodesic runs along equator
    calp1 = calp2 = 0.0
    salp1 = salp2 = 1.0
    s12x = self.a * lam12
    sig12 = omg12 = lam12 / self._f1
    m12x = self._b * sin(sig12)
    if (outmask & GEODESICSCALE) > 0
      M12 = M21 = cos(sig12)
    end
    a12 = lon12 / self._f1

  elseif ! meridian

    # Now point1 and point2 belong within a hemisphere bounded by a
    # meridian and geodesic is neither meridional or equatorial.

    # Figure a starting point for Newton's method
    sig12, salp1, calp1, salp2, calp2, dnm = _InverseStart(self,
      sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, C1a, C2a)

    if sig12 >= 0
      # Short lines (InverseStart sets salp2, calp2, dnm)
      s12x = sig12 * self._b * dnm
      m12x = (Math.sq(dnm) * self._b * sin(sig12 / dnm))
      if (outmask & GEODESICSCALE) > 0
        M12 = M21 = cos(sig12 / dnm)
      end
      a12 = rad2deg(sig12)
      omg12 = lam12 / (self._f1 * dnm)
    else

      # Newton's method.  This is a straightforward solution of f(alp1) =
      # lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
      # root in the interval (0, pi) and its derivative is positive at the
      # root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
      # alp1.  During the course of the iteration, a range (alp1a, alp1b) is
      # maintained which brackets the root and with each evaluation of f(alp)
      # the range is shrunk if possible.  Newton's method is restarted
      # whenever the derivative of f is negative (because the new value of
      # alp1 is then further from the solution) or if the new estimate of
      # alp1 lies outside (0,pi); in this case, the new starting guess is
      # taken to be (alp1a + alp1b) / 2.
      # real ssig1, csig1, ssig2, csig2, eps
      numit = 0
      tripn = tripb = false
      # Bracketing range
      salp1a = tiny_
      calp1a = 1.0
      salp1b = tiny_
      calp1b = -1.0

      while numit < maxit2_
        # the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
        # WGS84 and random input: mean = 2.85, sd = 0.60
        (v, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
         eps, domg12, dv) = _Lambda12(self,
           sbet1, cbet1, dn1, sbet2, cbet2, dn2,
           salp1, calp1, slam12, clam12, numit < maxit1_,
           C1a, C2a, C3a)
        # 2 * tol0 is approximately 1 ulp for a number in [0, pi].
        # Reversed test to allow escape with NaNs
        if tripb || !(abs(v) >= (tripn ? 8 : 1) * tol0_)
          break
        end
        # Update bracketing values
        if v > 0 && (numit > maxit1_ ||
                      calp1/salp1 > calp1b/salp1b)
          salp1b = salp1; calp1b = calp1
        elseif v < 0 && (numit > maxit1_ ||
                        calp1/salp1 < calp1a/salp1a)
          salp1a = salp1; calp1a = calp1
        end

        numit += 1
        if numit < maxit1_ && dv > 0
          dalp1 = -v/dv
          sdalp1 = sin(dalp1)
          cdalp1 = cos(dalp1)
          nsalp1 = salp1 * cdalp1 + calp1 * sdalp1
          if nsalp1 > 0 && abs(dalp1) < pi
            calp1 = calp1 * cdalp1 - salp1 * sdalp1
            salp1 = nsalp1
            salp1, calp1 = Math.norm(salp1, calp1)
            # In some regimes we don't get quadratic convergence because
            # slope -> 0.  So use convergence conditions based on epsilon
            # instead of sqrt(epsilon).
            tripn = abs(v) <= 16 * tol0_
            continue
          end
        end
        # Either dv was not positive or updated value was outside
        # legal range.  Use the midpoint of the bracket as the next
        # estimate.  This mechanism is not needed for the WGS84
        # ellipsoid, but it does catch problems with more eccentric
        # ellipsoids.  Its efficacy is such for
        # the WGS84 test set with the starting guess set to alp1 = 90deg:
        # the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
        # WGS84 and random input: mean = 4.74, sd = 0.99
        salp1 = (salp1a + salp1b)/2
        calp1 = (calp1a + calp1b)/2
        salp1, calp1 = Math.norm(salp1, calp1)
        tripn = false
        tripb = (abs(salp1a - salp1) + (calp1a - calp1) < tolb_ ||
                 abs(salp1 - salp1b) + (calp1 - calp1b) < tolb_)
      end

      lengthmask = outmask
      if (outmask & (REDUCEDLENGTH | GEODESICSCALE)) > 0
          lengthmask |= DISTANCE
      else
          lengthmask |= EMPTY
      end
      s12x, m12x, dummy, M12, M21 = _Lengths(self,
        eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
        lengthmask, C1a, C2a)

      m12x *= self._b
      s12x *= self._b
      a12 = rad2deg(sig12)
      if (outmask & AREA) > 0
        # omg12 = lam12 - domg12
        sdomg12 = sin(domg12)
        cdomg12 = cos(domg12)
        somg12 = slam12 * cdomg12 - clam12 * sdomg12
        comg12 = clam12 * cdomg12 + slam12 * sdomg12
      end
    end
  end
  # end elif not meridian

  if (outmask & DISTANCE) > 0
    s12 = 0.0 + s12x          # Convert -0 to 0
  end

  if (outmask & REDUCEDLENGTH) > 0
    m12 = 0.0 + m12x          # Convert -0 to 0
  end

  if (outmask & AREA) > 0
    # From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
    salp0 = salp1 * cbet1
    calp0 = hypot(calp1, salp1 * sbet1) # calp0 > 0
    # real alp12
    if calp0 != 0 && salp0 != 0
      # From Lambda12: tan(bet) = tan(sig) * cos(alp)
      ssig1 = sbet1; csig1 = calp1 * cbet1
      ssig2 = sbet2; csig2 = calp2 * cbet2
      k2 = Math.sq(calp0) * self._ep2
      eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)
      # Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
      A4 = Math.sq(self.a) * calp0 * salp0 * self._e2
      ssig1, csig1 = Math.norm(ssig1, csig1)
      ssig2, csig2 = Math.norm(ssig2, csig2)
      C4a = Vector{Float64}(undef, nC4_)
      _C4f(self, eps, C4a)
      B41 = _SinCosSeries(false, ssig1, csig1, C4a)
      B42 = _SinCosSeries(false, ssig2, csig2, C4a)
      S12 = A4 * (B42 - B41)
    else
      # Avoid problems with indeterminate sig1, sig2 on equator
      S12 = 0.0
    end

    if ! meridian && somg12 > 1
      somg12 = sin(omg12)
      comg12 = cos(omg12)
    end

    if (!meridian &&
        # omg12 < 3/4 * pi
        comg12 > -0.7071 &&   # Long difference not too big
        sbet2 - sbet1 < 1.75) # Lat difference not too big
      # Use tan(Gamma/2) = tan(omg12/2)
      # * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
      # with tan(x/2) = sin(x)/(1+cos(x))
      domg12 = 1 + comg12; dbet1 = 1 + cbet1; dbet2 = 1 + cbet2
      alp12 = 2 * atan( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                        domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) )
    else
      # alp12 = alp2 - alp1, used in atan2 so no need to normalize
      salp12 = salp2 * calp1 - calp2 * salp1
      calp12 = calp2 * calp1 + salp2 * salp1
      # The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
      # salp12 = -0 and alp12 = -180.  However this depends on the sign
      # being attached to 0 correctly.  The following ensures the correct
      # behavior.
      if salp12 == 0 && calp12 < 0
        salp12 = tiny_ * calp1
        calp12 = -1.0
      end
      alp12 = atan(salp12, calp12)
    end
    S12 += self._c2 * alp12
    S12 *= swapp * lonsign * latsign
    # Convert -0 to 0
    S12 += 0.0
  end

  # Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
  if swapp < 0
    salp2, salp1 = salp1, salp2
    calp2, calp1 = calp1, calp2
    if (outmask & GEODESICSCALE) > 0
      M21, M12 = M12, M21
    end
  end

  salp1 *= swapp * lonsign
  calp1 *= swapp * latsign
  salp2 *= swapp * lonsign
  calp2 *= swapp * latsign

  a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12
end

end # module
