# Ported to Julia by Andy Ferris, 2016, and re-released under an MIT license.
#/**
# * Copyright (c) Charles Karney (2008-2015) <charles@karney.com> and licensed
# * under the MIT/X11 License.  For more information, see
# * http://geographiclib.sourceforge.net/
# **********************************************************************/

"""
    immutable TransverseMercator{MaxPow}
    TransverseMercator(datum)
    TransverseMercator(ellipsoid)
    TransverseMercator(a, f, [::Type{Val{MaxPow}} = Val{6}])

Cache of ellipsoidal calculations necessary for transverse-Mercator and
polar-stereographic transformations. Series expansion coefficients up to order
`MaxPow` (between 4 and 8, default 6) are calculated and stored for fast
transverse-Mercator and UTM calculations.
"""
immutable TransverseMercator{MaxPow}
    a::Float64
    f::Float64
    e2::Float64
    es::Float64
    e2m::Float64
    c::Float64
    n::Float64
    a1::Float64
    b1::Float64
    alp::NTuple{MaxPow, Float64}
    bet::NTuple{MaxPow, Float64}
end

TransverseMercator(datum) = TransverseMercator(ellipsoid(datum))
TransverseMercator(el::Ellipsoid) = TransverseMercator(el.a, 1-el.b/el.a)
TransverseMercator(a,f) = TransverseMercator(a,f,Val{6}) # Default to sixth-order expansion
function TransverseMercator{MaxPow}(a::Float64, f::Float64, ::Type{Val{MaxPow}})
    e2 = f * (2 - f)
    es = (f < 0 ? -1 : 1) * sqrt(abs(e2))
    e2m = (1 - e2)
    # c = sqrt( pow(1 + e, 1 + e) * pow(1 - e, 1 - e) ) )
    # See, for example, Lee (1976), p 100.
    c = sqrt(e2m) * exp(eatanhe(1.0, e2)) # Maybe revert to C++ version with es
    n = f / (2 - f)

    (1 - f) * exp(eatanhe(1.0, es))

    if !(isa(MaxPow,Int)) || MaxPow < 4 || MaxPow > 8
        error("MaxPow must be 4, 5, 6, 7 or 8")
    end

    # Stack-allocated constants
    if MaxPow == 4
        b1coeff = (1, 16, 64, 64)

        alpcoeff = (# alp[1]/n^1, polynomial in n of order 3
                    164, 225, -480, 360, 720,
                    # alp[2]/n^2, polynomial in n of order 2
                    557, -864, 390, 1440,
                    # alp[3]/n^3, polynomial in n of order 1
                    -1236, 427, 1680,
                    # alp[4]/n^4, polynomial in n of order 0
                    49561, 161280)

        betcoeff = (# bet[1]/n^1, polynomial in n of order 3
                    -4, 555, -960, 720, 1440,
                    # bet[2]/n^2, polynomial in n of order 2
                    -437, 96, 30, 1440,
                    # bet[3]/n^3, polynomial in n of order 1
                    -148, 119, 3360,
                    # bet[4]/n^4, polynomial in n of order 0
                    4397, 161280)
    elseif MaxPow == 5
        b1coeff = (1, 16, 64, 64)

        alpcoeff = (# alp[1]/n^1, polynomial in n of order 4
                    -635, 328, 450, -960, 720, 1440,
                    # alp[2]/n^2, polynomial in n of order 3
                    4496, 3899, -6048, 2730, 10080,
                    # alp[3]/n^3, polynomial in n of order 2
                    15061, -19776, 6832, 26880,
                    # alp[4]/n^4, polynomial in n of order 1
                    -171840, 49561, 161280,
                    # alp[5]/n^5, polynomial in n of order 0
                    34729, 80640)

        betcoeff = (# bet[1]/n^1, polynomial in n of order 4
                    -3645, -64, 8880, -15360, 11520, 23040,
                    # bet[2]/n^2, polynomial in n of order 3
                    4416, -3059, 672, 210, 10080,
                    # bet[3]/n^3, polynomial in n of order 2
                    -627, -592, 476, 13440,
                    # bet[4]/n^4, polynomial in n of order 1
                    -3520, 4397, 161280,
                    # bet[5]/n^5, polynomial in n of order 0
                    4583, 161280)
    elseif MaxPow == 6
        b1coeff = (1, 4, 64, 256, 256)

        alpcoeff = (# alp[1]/n^1, polynomial in n of order 5
                    31564, -66675, 34440, 47250, -100800, 75600, 151200,
                    # alp[2]/n^2, polynomial in n of order 4
                    -1983433, 863232, 748608, -1161216, 524160, 1935360,
                    # alp[3]/n^3, polynomial in n of order 3
                    670412, 406647, -533952, 184464, 725760,
                    # alp[4]/n^4, polynomial in n of order 2
                    6601661, -7732800, 2230245, 7257600,
                    # alp[5]/n^5, polynomial in n of order 1
                    -13675556, 3438171, 7983360,
                    # alp[6]/n^6, polynomial in n of order 0
                    212378941, 319334400)

        betcoeff = (# bet[1]/n^1, polynomial in n of order 5
                    384796, -382725, -6720, 932400, -1612800, 1209600, 2419200,
                    # bet[2]/n^2, polynomial in n of order 4
                    -1118711, 1695744, -1174656, 258048, 80640, 3870720,
                    # bet[3]/n^3, polynomial in n of order 3
                    22276, -16929, -15984, 12852, 362880,
                    # bet[4]/n^4, polynomial in n of order 2
                    -830251, -158400, 197865, 7257600,
                    # bet[5]/n^5, polynomial in n of order 1
                    -435388, 453717, 15966720,
                    # bet[6]/n^6, polynomial in n of order 0
                    20648693, 638668800)
    elseif MaxPow == 7
        b1coeff = (1, 4, 64, 256, 256)

        alpcoeff = (# alp[1]/n^1, polynomial in n of order 6
                    1804025, 2020096, -4267200, 2204160, 3024000, -6451200, 4838400, 9676800,
                    # alp[2]/n^2, polynomial in n of order 5
                    4626384, -9917165, 4316160, 3743040, -5806080, 2620800, 9676800,
                    # alp[3]/n^3, polynomial in n of order 4
                    -67102379, 26816480, 16265880, -21358080, 7378560, 29030400,
                    # alp[4]/n^4, polynomial in n of order 3
                    155912000, 72618271, -85060800, 24532695, 79833600,
                    # alp[5]/n^5, polynomial in n of order 2
                    102508609, -109404448, 27505368, 63866880,
                    # alp[6]/n^6, polynomial in n of order 1
                    -12282192400, 2760926233, 4151347200,
                    # alp[7]/n^7, polynomial in n of order 0
                    1522256789, 1383782400)

        betcoeff = (# bet[1]/n^1, polynomial in n of order 6
                    -5406467, 6156736, -6123600, -107520, 14918400, -25804800, 19353600,
                    38707200,
                    # bet[2]/n^2, polynomial in n of order 5
                    829456, -5593555, 8478720, -5873280, 1290240, 403200, 19353600,
                    # bet[3]/n^3, polynomial in n of order 4
                    9261899, 3564160, -2708640, -2557440, 2056320, 58060800,
                    # bet[4]/n^4, polynomial in n of order 3
                    14928352, -9132761, -1742400, 2176515, 79833600,
                    # bet[5]/n^5, polynomial in n of order 2
                    -8005831, -1741552, 1814868, 63866880,
                    # bet[6]/n^6, polynomial in n of order 1
                    -261810608, 268433009, 8302694400,
                    # bet[7]/n^7, polynomial in n of order 0
                    219941297, 5535129600)
    elseif MaxPow == 8
        b1coeff = (25, 64, 256, 4096, 16384, 16384)

        alpcoeff = (# alp[1]/n^1, polynomial in n of order 7
                    -75900428, 37884525, 42422016, -89611200, 46287360, 63504000, -135475200,
                    101606400, 203212800,
                    # alp[2]/n^2, polynomial in n of order 6
                    148003883, 83274912, -178508970, 77690880, 67374720, -104509440,
                    47174400, 174182400,
                    # alp[3]/n^3, polynomial in n of order 5
                    318729724, -738126169, 294981280, 178924680, -234938880, 81164160,
                    319334400,
                    # alp[4]/n^4, polynomial in n of order 4
                    -40176129013, 14967552000, 6971354016, -8165836800, 2355138720,
                    7664025600,
                    # alp[5]/n^5, polynomial in n of order 3
                    10421654396, 3997835751, -4266773472, 1072709352, 2490808320,
                    # alp[6]/n^6, polynomial in n of order 2
                    175214326799, -171950693600, 38652967262, 58118860800,
                    # alp[7]/n^7, polynomial in n of order 1
                    -67039739596, 13700311101, 12454041600,
                    # alp[8]/n^8, polynomial in n of order 0
                    1424729850961, 743921418240)

        betcoeff = (# bet[1]/n^1, polynomial in n of order 7
                    31777436, -37845269, 43097152, -42865200, -752640, 104428800, -180633600,
                    135475200, 270950400,
                    # bet[2]/n^2, polynomial in n of order 6
                    24749483, 14930208, -100683990, 152616960, -105719040, 23224320, 7257600,
                    348364800,
                    # bet[3]/n^3, polynomial in n of order 5
                    -232468668, 101880889, 39205760, -29795040, -28131840, 22619520,
                    638668800,
                    # bet[4]/n^4, polynomial in n of order 4
                    324154477, 1433121792, -876745056, -167270400, 208945440, 7664025600,
                    # bet[5]/n^5, polynomial in n of order 3
                    457888660, -312227409, -67920528, 70779852, 2490808320,
                    # bet[6]/n^6, polynomial in n of order 2
                    -19841813847, -3665348512, 3758062126, 116237721600,
                    # bet[7]/n^7, polynomial in n of order 1
                    -1989295244, 1979471673, 49816166400,
                    # bet[8]/n^8, polynomial in n of order 0
                    191773887257, 3719607091200)
    end

    m = div(MaxPow, 2)
    b1 = polyval(b1coeff[1:m+1], n*n) / (b1coeff[m+2] * (1 + n))

    # a1 is the equivalent radius for computing the circumference of ellipse.
    a1 = b1 * a
    o = 1
    d = n
    _alp = Vector{Float64}(MaxPow)
    _bet = Vector{Float64}(MaxPow)
    for l = 1:MaxPow
        m = MaxPow - l
        _alp[l] = d * polyval(alpcoeff[o:o+m], n) / alpcoeff[o + m + 1]
        _bet[l] = d * polyval(betcoeff[o:o+m], n) / betcoeff[o + m + 1]
        o += m + 2
        d *= n
    end
    # Post condition: o == sizeof(alpcoeff) / sizeof(real) &&
    # o == sizeof(betcoeff) / sizeof(real)
    alp = (_alp...)
    bet = (_bet...)

    return TransverseMercator{MaxPow}(a,f,e2,es,e2m,c,n,a1,b1,alp,bet)
end


# // Engsager and Poder (2007) use trigonometric series to convert between phi
# // and phip.  Here are the series...
# //
# // Conversion from phi to phip:
# //
# //     phip = phi + sum(c[j] * sin(2*j*phi), j, 1, 6)
# //
# //       c[1] = - 2 * n
# //              + 2/3 * n^2
# //              + 4/3 * n^3
# //              - 82/45 * n^4
# //              + 32/45 * n^5
# //              + 4642/4725 * n^6;
# //       c[2] =   5/3 * n^2
# //              - 16/15 * n^3
# //              - 13/9 * n^4
# //              + 904/315 * n^5
# //              - 1522/945 * n^6;
# //       c[3] = - 26/15 * n^3
# //              + 34/21 * n^4
# //              + 8/5 * n^5
# //              - 12686/2835 * n^6;
# //       c[4] =   1237/630 * n^4
# //              - 12/5 * n^5
# //              - 24832/14175 * n^6;
# //       c[5] = - 734/315 * n^5
# //              + 109598/31185 * n^6;
# //       c[6] =   444337/155925 * n^6;
# //
# // Conversion from phip to phi:
# //
# //     phi = phip + sum(d[j] * sin(2*j*phip), j, 1, 6)
# //
# //       d[1] =   2 * n
# //              - 2/3 * n^2
# //              - 2 * n^3
# //              + 116/45 * n^4
# //              + 26/45 * n^5
# //              - 2854/675 * n^6;
# //       d[2] =   7/3 * n^2
# //              - 8/5 * n^3
# //              - 227/45 * n^4
# //              + 2704/315 * n^5
# //              + 2323/945 * n^6;
# //       d[3] =   56/15 * n^3
# //              - 136/35 * n^4
# //              - 1262/105 * n^5
# //              + 73814/2835 * n^6;
# //       d[4] =   4279/630 * n^4
# //              - 332/35 * n^5
# //              - 399572/14175 * n^6;
# //       d[5] =   4174/315 * n^5
# //              - 144838/6237 * n^6;
# //       d[6] =   601676/22275 * n^6;
# //
# // In order to maintain sufficient relative accuracy close to the pole use
# //
# //     S = sum(c[i]*sin(2*i*phi),i,1,6)
# //     taup = (tau + tan(S)) / (1 - tau * tan(S))
#
# // In Math::taupf and Math::tauf we evaluate the forward transform explicitly
# // and solve the reverse one by Newton's method.
# //
# // There are adapted from TransverseMercatorExact (taup and taupinv).  tau =
# // tan(phi), taup = sinh(psi)

"""
    (x, y, γ, k) = transverse_mercator_forward(lon0, lat, lon, k0, tm::TransverseMercator)

Perform transverse-Mercator projection of `lat` and `lon` with respect to reference
meridian `lat0` and horizontal scaling `k0` (`= 0.9996` for UTM) using a
series expansion approach (see `TransverseMercator`). `γ` and `k` are the local
convergence and scaling factors, respectively.
"""
function transverse_mercator_forward{MaxPow}(lon0, lat, lon, k0, tm::TransverseMercator{MaxPow})
    lat = LatFix(lat)
    lon = AngDiff(lon0, lon)
    # Explicitly enforce the parity
    latsign = lat < 0 ? -1 : 1
    lonsign = lon < 0 ? -1 : 1
    lon *= lonsign
    lat *= latsign
    backside = lon > 90
    if backside
        if (lat == 0)
            latsign = -1
        end
        lon = 180 - lon
    end

    sphi = sind(lat)
    cphi = cosd(lat)
    slam = sind(lon)
    clam = cosd(lon)

    # phi = latitude
    # phi' = conformal latitude
    # psi = isometric latitude
    # tau = tan(phi)
    # tau' = tan(phi')
    # [xi', eta'] = Gauss-Schreiber TM coordinates
    # [xi, eta] = Gauss-Krueger TM coordinates
    #
    # We use
    #   tan(phi') = sinh(psi)
    #   sin(phi') = tanh(psi)
    #   cos(phi') = sech(psi)
    #   denom^2    = 1-cos(phi')^2*sin(lam)^2 = 1-sech(psi)^2*sin(lam)^2
    #   sin(xip)   = sin(phi')/denom          = tanh(psi)/denom
    #   cos(xip)   = cos(phi')*cos(lam)/denom = sech(psi)*cos(lam)/denom
    #   cosh(etap) = 1/denom                  = 1/denom
    #   sinh(etap) = cos(phi')*sin(lam)/denom = sech(psi)*sin(lam)/denom
    if (lat != 90)
        tau = sphi / cphi
        taup = taupf(tau, tm.e2) # TODO maybe switch to the C++ version of taupf using tm.es?
        xip = atan2(taup, clam)
        # Used to be
        #   etap = Math::atanh(sin(lam) / cosh(psi));
        etap = asinh(slam / hypot(taup, clam))
        # convergence and scale for Gauss-Schreiber TM (xip, etap) -- gamma0 =
        # atan(tan(xip) * tanh(etap)) = atan(tan(lam) * sin(phi'));
        # sin(phi') = tau'/sqrt(1 + tau'^2)
        # Krueger p 22 (44)
        gamma = atan2d(slam * taup, clam * hypot(1.0, taup))
        # k0 = sqrt(1 - _e2 * sin(phi)^2) * (cos(phi') / cos(phi)) * cosh(etap)
        # Note 1/cos(phi) = cosh(psip);
        # and cos(phi') * cosh(etap) = 1/hypot(sinh(psi), cos(lam))
        #
        # This form has cancelling errors.  This property is lost if cosh(psip)
        # is replaced by 1/cos(phi), even though it's using "primary" data (phi
        # instead of psip).
        k = sqrt(tm.e2m + tm.e2 * cphi*cphi) * hypot(1.0, tau) / hypot(taup, clam)
    else
        xip = pi/2
        etap = 0.0
        gamma = lon
        k = tm.c
    end

    # {xi',eta'} is {northing,easting} for Gauss-Schreiber transverse Mercator
    # (for eta' = 0, xi' = bet). {xi,eta} is {northing,easting} for transverse
    # Mercator with constant scale on the central meridian (for eta = 0, xip =
    # rectifying latitude).  Define
    #
    #   zeta = xi + i*eta
    #   zeta' = xi' + i*eta'
    #
    # The conversion from conformal to rectifying latitude can be expressed as
    # a series in _n:
    #
    #   zeta = zeta' + sum(h[j-1]' * sin(2 * j * zeta'), j = 1..maxpow_)
    #
    # where h[j]' = O(_n^j).  The reversion of this series gives
    #
    #   zeta' = zeta - sum(h[j-1] * sin(2 * j * zeta), j = 1..maxpow_)
    #
    # which is used in Reverse.
    #
    # Evaluate sums via Clenshaw method.  See
    #    http:#mathworld.wolfram.com/ClenshawRecurrenceFormula.html
    #
    # Let
    #
    #    S = sum(c[k] * F[k](x), k = 0..N)
    #    F[n+1](x) = alpha(n,x) * F[n](x) + beta(n,x) * F[n-1](x)
    #
    # Evaluate S with
    #
    #    y[N+2] = y[N+1] = 0
    #    y[k] = alpha(k,x) * y[k+1] + beta(k+1,x) * y[k+2] + c[k]
    #    S = c[0] * F[0](x) + y[1] * F[1](x) + beta(1,x) * F[0](x) * y[2]
    #
    # Here we have
    #
    #    x = 2 * zeta'
    #    F[n](x) = sin(n * x)
    #    a(n, x) = 2 * cos(x)
    #    b(n, x) = -1
    #    [ sin(A+B) - 2*cos(B)*sin(A) + sin(A-B) = 0, A = n*x, B = x ]
    #    N = maxpow_
    #    c[k] = _alp[k]
    #    S = y[1] * sin(x)
    #
    # For the derivative we have
    #
    #    x = 2 * zeta'
    #    F[n](x) = cos(n * x)
    #    a(n, x) = 2 * cos(x)
    #    b(n, x) = -1
    #    [ cos(A+B) - 2*cos(B)*cos(A) + cos(A-B) = 0, A = n*x, B = x ]
    #    c[0] = 1; c[k] = 2*k*_alp[k]
    #    S = (c[0] - y[2]) + y[1] * cos(x)
    c0 = cos(2 * xip)
    ch0 = cosh(2 * etap)
    s0 = sin(2 * xip)
    sh0 = sinh(2 * etap)
    ar = 2 * c0 * ch0
    ai = -2 * s0 * sh0 # 2 * cos(2*zeta')

    n = MaxPow

    xi0 = (isodd(n) ? tm.alp[n] : 0.0)
    eta0 = 0.0
    xi1 = 0.0
    eta1 = 0.0

    # Accumulators for dzeta/dzeta'
    if isodd(n)
        yr0 = 2 * MaxPow * tm.alp[n]
        n = n - 1
    else
        yr0 = 0.0
    end
    yi0 = 0.0
    yr1 = 0.0
    yi1 = 0.0
    while n > 0
        xi1  = ar * xi0 - ai * eta0 - xi1 + tm.alp[n]
        eta1 = ai * xi0 + ar * eta0 - eta1
        yr1 = ar * yr0 - ai * yi0 - yr1 + 2 * n * tm.alp[n]
        yi1 = ai * yr0 + ar * yi0 - yi1

        n = n - 1

        xi0  = ar * xi1 - ai * eta1 - xi0 + tm.alp[n]
        eta0 = ai * xi1 + ar * eta1 - eta0
        yr0 = ar * yr1 - ai * yi1 - yr0 + 2 * n * tm.alp[n]
        yi0 = ai * yr1 + ar * yi1 - yi0

        n = n - 1
    end

    ar *= 0.5
    ai *= 0.5      # cos(2*zeta')
    yr1 = 1 - yr1 + ar * yr0 - ai * yi0
    yi1 =   - yi1 + ai * yr0 + ar * yi0
    ar = s0 * ch0
    ai = c0 * sh0  # sin(2*zeta')

    xi  = xip  + ar * xi0 - ai * eta0
    eta = etap + ai * xi0 + ar * eta0

    # Fold in change in convergence and scale for Gauss-Schreiber TM to
    # Gauss-Krueger TM.

    gamma -= atan2d(yi1, yr1)
    k *= tm.b1 * hypot(yr1, yi1)
    y = tm.a1 * k0 * (backside ? pi - xi : xi) * latsign
    x = tm.a1 * k0 * eta * lonsign
    if backside
        gamma = 180 - gamma
    end
    gamma *= latsign * lonsign
    gamma = AngNormalize(gamma)
    k *= k0

    return (x, y, gamma, k)
end

"""
    (lat, lon, γ, k) = transverse_mercator_reverse(lon0, x, y, k0, tm::TransverseMercator)

Invert transverse-Mercator projection of `x` and `y` with respect to reference
meridian `lat0` and horizontal scaling `k0` (`= 0.9996` for UTM) using a
series expansion approach (see `TransverseMercator`). `γ` and `k` are the local
convergence and scaling factors, respectively.
"""
function transverse_mercator_reverse{MaxPow}(lon0, x, y, k0, tm::TransverseMercator{MaxPow})
    # This undoes the steps in transverse_mercator_forward.  The wrinkles are: (1) Use of the
    # reverted series to express zeta' in terms of zeta. (2) Newton's method
    # to solve for phi in terms of tan(phi).
    xi = y / (tm.a1 * k0)
    eta = x / (tm.a1 * k0)

    # Explicitly enforce the parity
    xisign = xi < 0 ? -1 : 1
    etasign = eta < 0 ? -1 : 1
    xi *= xisign
    eta *= etasign
    backside = xi > pi/2
    if backside
        xi = pi - xi
    end

    c0 = cos(2 * xi)
    ch0 = cosh(2 * eta)
    s0 = sin(2 * xi)
    sh0 = sinh(2 * eta)

    ar = 2 * c0 * ch0
    ai = -2 * s0 * sh0 # 2 * cos(2*zeta)
    n = MaxPow

    # Accumulators for zeta'
    xip0 = (isodd(n) ? -tm.bet[n] : 0.0)
    etap0 = 0.0
    xip1 = 0.0
    etap1 = 0.0

    # Accumulators for dzeta'/dzeta
    if isodd(n)
        yr0 =  2 * MaxPow * tm.bet[n]
        n -= 1
    else
        yr0 = 0.0
    end
    yi0 = 0.0
    yr1 = 0.0
    yi1 = 0.0

    while (n > 0)
        xip1  = ar * xip0 - ai * etap0 - xip1 - tm.bet[n]
        etap1 = ai * xip0 + ar * etap0 - etap1
        yr1 = ar * yr0 - ai * yi0 - yr1 - 2 * n * tm.bet[n]
        yi1 = ai * yr0 + ar * yi0 - yi1

        n = n - 1

        xip0  = ar * xip1 - ai * etap1 - xip0 - tm.bet[n]
        etap0 = ai * xip1 + ar * etap1 - etap0
        yr0 = ar * yr1 - ai * yi1 - yr0 - 2 * n * tm.bet[n]
        yi0 = ai * yr1 + ar * yi1 - yi0

        n = n - 1
    end

    ar *= 0.5
    ai *= 0.5  # cos(2*zeta')
    yr1 = 1 - yr1 + ar * yr0 - ai * yi0
    yi1 =   - yi1 + ai * yr0 + ar * yi0
    ar = s0 * ch0
    ai = c0 * sh0 # sin(2*zeta)

    xip  = xi  + ar * xip0 - ai * etap0
    etap = eta + ai * xip0 + ar * etap0

    # Convergence and scale for Gauss-Schreiber TM to Gauss-Krueger TM.
    gamma = atan2(yi1, yr1) * 180/pi
    k = tm.b1 / hypot(yr1, yi1)
    # JHS 154 has
    #
    #   phi' = asin(sin(xi') / cosh(eta')) (Krueger p 17 (25))
    #   lam = asin(tanh(eta') / cos(phi')
    #   psi = asinh(tan(phi'))

    s = sinh(etap)
    c = max(0.0, cos(xip))  # cos(pi/2) might be negative
    r = hypot(s, c)

    if (r != 0)
        lon = atan2(s, c) * 180/pi # Krueger p 17 (25)

        # Use Newton's method to solve for tau
        sxip = sin(xip)
        tau = tauf(sxip/r, tm.e2) # TODO maybe change to C++ es version
        gamma += atan2(sxip * tanh(etap), c) * 180/pi # Krueger p 19 (31)
        lat = atand(tau)
        # Note cos(phi') * cosh(eta') = r
        k *= sqrt(tm.e2m + tm.e2 / (1 + tau*tau)) * hypot(1.0, tau) * r
    else
        lat = 90.0
        lon = 0.0
        k *= tm.c
    end

    lat *= xisign
    if backside
        lon = 180 - lon
    end
    lon *= etasign
    lon = AngNormalize(lon + lon0)
    if backside
        gamma = 180 - gamma
    end
    gamma *= xisign * etasign
    gamma = AngNormalize(gamma)
    k *= k0

    return (lat, lon, gamma, k)
end


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
%TAUF   tan(phi)
%
%   TAUF(taup, e2) returns tangent of phi in terms of taup the tangent of
%   chi.  e2, the square of the eccentricity, is a scalar; taup can be any
%   shape.
"""
function tauf(taup, e2)
    numit = 5
    e2m = 1 - e2     # This possibly disagrees with the C++ code in the prolate case
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
        dtau = (taup - taupa) * (1 + e2m * tau*tau) / (e2m * tau1 * hypot(1, taupa))
        tau = tau + dtau
        g = g & (abs(dtau) >= stol)
    end
    return tau
end


"""
taupf(tau, e2) returns tangent of chi in terms of tau the tangent of
phi.  e2, the square of the eccentricity, is a scalar; tau can be any
shape.
"""
function taupf(tau, e2)
  tau1 = hypot(1, tau)
  sig = sinh( eatanhe( tau / tau1, e2 ) )
  return hypot(1, sig) * tau - sig * tau1
end









# Code adapted from MATLAB port of geograhiclib, which was licensed under
# a permissive MIT-like license and the following copyright notice
#     Copyright (c) 2016, Charles Karney
#     All rights reserved.

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


function AngDiff(x,y)
    AngNormalize(AngNormalize(y) - AngNormalize(x))
    #(d, t) = sumx(AngNormalize(x), AngNormalize(-y))
    #d = - AngNormalize(d)
    #if d == 180 && t < 0
    #    d = -180
    #end
    #return sumx(d, -t)[1]
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

function atan2d(y, x)
    # In order to minimize round-off errors, this function rearranges the
    # arguments so that result of atan2 is in the range [-pi/4, pi/4] before
    # converting it to degrees and mapping the result to the correct
    # quadrant.
    q = 0
    if (abs(y) > abs(x))
        tmp = x
        x = y
        y = tmp
        q = 2
    end
    if (x < 0)
        x = -x
        q = q + 1
    end
    # here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4]
    ang = atan2(y, x) * 180 / pi
    # Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
    # atan2d will not be called with y = -0.  If need be, include
    #
    #   case 0: ang = 0 + ang; break;
    #
    # and handle mpfr as in AngRound.
    if q == 1
        ang = (y > 0 ? 180 : -180) - ang
    elseif q == 2
        ang =  90 - ang
    elseif q == 3
        ang = -90 + ang
    end
    return ang
end




#=
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
=#
#=
"""
AngDiffE  Compute angle difference accurately

   (d, e) = AngDiff(x, y) computes z = y - x, reduced to (-180,180].  d =
   round(z) and e = z - round(z).
"""
function AngDiffE(x, y)
    (d, t) = sumx(AngNormalize(x), AngNormalize(-y))
    d = - AngNormalize(d)
    if d == 180 && t < 0
        d = -180
    end
    return sumx(d, -t)
end
=#
#=
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
=#

#=
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
=#
#=
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
=#
#=
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
=#
#=
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
=#
#=

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





=#
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
