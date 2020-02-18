module Math

const digits = 53
const epsilon = eps()
const minval = 2.0^-1022
const maxval = 2.0^1023 * (2 - epsilon)
const inf = Inf
const nan = NaN

"""Square a number"""
sq(x) = x*x

"""Real cube root of a number"""
function cbrt(x)
    y = abs(x)^(1/3.0)
    x >= 0 ? y : -y
end

"""log(1 + x) accurate for small x (missing from python 2.5.2)"""
log1p(x) = Base.log1p(x)

"""atanh(x) (missing from python 2.5.2)"""
function atanh(x)
    y = abs(x)
    y = log1p(2 * y/(1 - y))/2
    x < 0 ? -y : y
end

"""return x with the sign of y (missing from python 2.5.2)"""
copysign(x, y) = Base.copysign(x, y)

"""Private: Normalize a two-vector"""
norm(x, y) = (r = hypot(x, y); (x/r, y/r))

"""Error free transformation of a sum.
Note that t can be the same as one of the first two arguments."""
function sum(u, v)
    s = u + v
    up = s - v
    vpp = s - up
    up -= u
    vpp -= v
    t = -(up + vpp)
    s, t 
end

"""Evaluate a polynomial"""
function polyval(N, p, s, x)
    y = float(N < 0 ? 0 : p[s+1])
    while N > 0
        N -= 1
        s += 1
        y = y*x + p[s+1]
    end
    y
end

"""
Private: Round an angle so that small values underflow to zero.
"""
function AngRound(x)
    # The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
    # for reals = 0.7 pm on the earth if x is an angle in degrees.  (This
    # is about 1000 times more resolution than we get with angles around 90
    # degrees.)  We use this to avoid having to deal with near singular
    # cases when x is non-zero but tiny (e.g., 1.0e-200).
    z = 1/16.0
    y = abs(x)
    # The compiler mustn't "simplify" z - (z - y) to y
    y < z && (y = z - (z - y))
    x < 0 ? zero(x) - y : y
end

"""reduce angle in (-180,180]"""
function AngNormalize(x)
    y = x % 360
    y = x == 0 ? x : y
    if y <= -180
        y + 360
    elseif y <= 180
        y
    else
        y - 360
    end
end

"""replace angles outside [-90,90] by NaN"""
LatFix(x) = abs(x) > 90 ? nan : x

"Compute y - x and reduce to [-180, 180]° accurately"
function AngDiff(x, y)
    d, t = sum(AngNormalize(-x), AngNormalize(y))
    d = AngNormalize(d)
    sum(d == 180 && t > 0 ? -180 : d, t)
end

"""Compute sine and cosine of x in degrees."""
function sincosd(x)
    r = x % 360
    q = isnan(r) ? nan : floor(Int, r / 90 + 0.5)
    r -= 90 * q
    r = deg2rad(r)
    s = sin(r)
    c = cos(r)
    q = mod(q, 4)
    if q == 1
        s, c = c, -s
    elseif q == 2
        s, c = -s, -c
    elseif q == 3
        s, c = -c, s
    end
    s, c = x == 0 ? (x, c) : (0.0+s, 0.0+c)
end

"""compute atan2(y, x) with the result in degrees"""
function atan2d(y::T1, x::T2) where {T1,T2}
    T = float(promote_type(T1, T2))
    if abs(y) > abs(x)
        q = 2
        x, y = y, x
    else
        q = 0
    end
    if x < 0
        q += 1
        x = -x
    end
    ang = rad2deg(atan(y, x))
    if q == 1
        ang = (y >= 0 ? 180 : -180) - ang
    elseif q == 2
        ang =  90 - ang
    elseif q == 3
        ang = -90 + ang
    end
    T(ang)
end

isfinite(x) = Base.isfinite(x)

isnan(x) = Base.isnan(x)

end # module
