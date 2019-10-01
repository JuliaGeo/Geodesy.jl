module GeodesicLines

import ..Math, ..GeodesicCapability, ..Geodesics, ..Result

export GeodesicLine, SetArc, SetDistance, Position, ArcPosition

struct GeodesicLine
    a::Float64
    f::Float64
    _b::Float64
    _c2::Float64
    _f1::Float64
    """the capabilities (readonly)"""
    caps::Int
    """the latitude of the first point in degrees (readonly)"""
    lat1::Float64
    """the longitude of the first point in degrees (readonly)"""
    lon1::Float64
    """the azimuth at the first point in degrees (readonly)"""
    azi1::Float64
    """the sine of the azimuth at the first point (readonly)"""
    salp1::Float64
    """the cosine of the azimuth at the first point (readonly)"""
    calp1::Float64
    _dn1::Float64
    _salp0::Float64
    _calp0::Float64
    _ssig1::Float64
    _csig1::Float64
    _somg1::Float64
    _comg1::Float64
    _k2::Float64
    _A1m1::Float64
    _C1a::Vector{Float64}
    _B11::Float64
    _stau1::Float64
    _ctau1::Float64
    _C1pa::Vector{Float64}
    _A2m1::Float64
    _C2a::Vector{Float64}
    _B21::Float64
    _C3a::Vector{Float64}
    _A3c::Float64
    _B31::Float64
    _C4a::Vector{Float64}
    _A4::Float64
    _B41::Float64
    """the distance between point 1 and point 3 in meters (readonly)"""
    s13::Float64
    """the arc length between point 1 and point 3 in degrees (readonly)"""
    a13::Float64
end

# Treat the GeodesicLine struct as a scalar
Base.Broadcast.broadcastable(line::GeodesicLine) = Ref(line)

"""
    GeodesicLine(geod::Geodesic, lat1, lon1, azi1; caps, salp1, calp1) -> line

Create a `GeodesicLine` with starting latitude `lat1`° and longitude `lon1`°,
and azimuth `azi1`°.

Control the capabilities of the `line` with `caps`.  Optionally specify
the sine and cosine of the azimuth at point 1, respectively `salp1` and `calp1`.
"""
function GeodesicLine(geod::Geodesics.Geodesic, lat1, lon1, azi1;
                      caps = GeodesicCapability.STANDARD | GeodesicCapability.DISTANCE_IN,
                      salp1 = Math.nan, calp1 = Math.nan)
    self_a = geod.a
    self_f = geod.f
    self__b = geod._b
    self__c2 = geod._c2
    self__f1 = geod._f1
    self_caps = (caps | Geodesics.LATITUDE | Geodesics.AZIMUTH |
                  Geodesics.LONG_UNROLL)

    # Guard against underflow in salp0
    self_lat1 = Math.LatFix(lat1)
    self_lon1 = lon1
    if Math.isnan(salp1) || Math.isnan(calp1)
      self_azi1 = Math.AngNormalize(azi1)
      self_salp1, self_calp1 = Math.sincosd(Math.AngRound(azi1))
    else
      self_azi1 = azi1
      self_salp1 = salp1
      self_calp1 = calp1
    end

    # real cbet1, sbet1
    sbet1, cbet1 = Math.sincosd(Math.AngRound(lat1))
    sbet1 *= self__f1
    # Ensure cbet1 = +epsilon at poles
    sbet1, cbet1 = Math.norm(sbet1, cbet1)
    cbet1 = max(Geodesics.tiny_, cbet1)
    self__dn1 = sqrt(1 + geod._ep2 * Math.sq(sbet1))

    # Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    self__salp0 = self_salp1 * cbet1 # alp0 in [0, pi/2 - |bet1|]
    # Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    # is slightly better (consider the case salp1 = 0).
    self__calp0 = hypot(self_calp1, self_salp1 * sbet1)
    # Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    # sig = 0 is nearest northward crossing of equator.
    # With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    # With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    # With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    # Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    # With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    # No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    # With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    self__ssig1 = sbet1
    self__somg1 = self__salp0 * sbet1
    self__csig1 = self__comg1 = ((sbet1 != 0 || self_calp1 != 0) ? cbet1 * self_calp1 : 1)
    # sig1 in (-pi, pi]
    self__ssig1, self__csig1 = Math.norm(self__ssig1, self__csig1)
    # No need to normalize
    # self__somg1, self__comg1 = Math.norm(self__somg1, self__comg1)

    self__k2 = Math.sq(self__calp0) * geod._ep2
    eps = self__k2 / (2 * (1 + sqrt(1 + self__k2)) + self__k2)

    if (self_caps & Geodesics.CAP_C1) > 0
      self__A1m1 = Geodesics._A1m1f(eps)
      self__C1a = Vector{Float64}(undef, Geodesics.nC1_ + 1) #list(range(Geodesics.nC1_ + 1))
      Geodesics._C1f(eps, self__C1a)
      self__B11 = Geodesics._SinCosSeries(true, self__ssig1, self__csig1, self__C1a)
      s = sin(self__B11)
      c = cos(self__B11)
      # tau1 = sig1 + B11
      self__stau1 = self__ssig1 * c + self__csig1 * s
      self__ctau1 = self__csig1 * c - self__ssig1 * s
      # Not necessary because C1pa reverts C1a
      #    _B11 = -_SinCosSeries(true, _stau1, _ctau1, _C1pa)
    else
      self__A1m1 = self__B11 = self__stau1 = self__ctau1 = Math.nan
      self__C1a = Float64[]
    end

    if (self_caps & Geodesics.CAP_C1p) > 0
      self__C1pa = Vector{Float64}(undef, Geodesics.nC1p_ + 1) #list(range(Geodesics.nC1p_ + 1))
      Geodesics._C1pf(eps, self__C1pa)
    else
      self__C1pa = Float64[]
    end

    if (self_caps & Geodesics.CAP_C2) > 0
      self__A2m1 = Geodesics._A2m1f(eps)
      self__C2a = Vector{Float64}(undef, Geodesics.nC2_ + 1) #list(range(Geodesics.nC2_ + 1))
      Geodesics._C2f(eps, self__C2a)
      self__B21 = Geodesics._SinCosSeries(true, self__ssig1, self__csig1, self__C2a)
    else
      self__A2m1 = self__B21 = Math.nan
      self__C2a = Float64[]
    end

    if (self_caps & Geodesics.CAP_C3) > 0
      self__C3a = Vector{Float64}(undef, Geodesics.nC3_) #list(range(Geodesics.nC3_))
      Geodesics._C3f(geod, eps, self__C3a)
      self__A3c = -self_f * self__salp0 * Geodesics._A3f(geod, eps)
      self__B31 = Geodesics._SinCosSeries(true, self__ssig1, self__csig1, self__C3a)
    else
        self__C3a = Float64[]
        self__A3c = self__B31 = Math.nan
    end

    if (self_caps & Geodesics.CAP_C4) > 0
      self__C4a = Vector{Float64}(undef, Geodesics.nC4_) #list(range(Geodesics.nC4_))
      Geodesics._C4f(geod, eps, self__C4a)
      # Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
      self__A4 = Math.sq(self_a) * self__calp0 * self__salp0 * geod._e2
      self__B41 = Geodesics._SinCosSeries(false, self__ssig1, self__csig1, self__C4a)
    else
      self__C4a = Float64[]
      self__A4 = self__B41 = Math.nan
    end

    self_s13 = Math.nan
    self_a13 = Math.nan

    GeodesicLine(self_a, self_f, self__b, self__c2, self__f1, self_caps,
                 self_lat1, self_lon1, self_azi1, self_salp1, self_calp1,
                 self__dn1, self__salp0, self__calp0, self__ssig1, self__csig1,
                 self__somg1, self__comg1, self__k2, self__A1m1, self__C1a,
                 self__B11, self__stau1, self__ctau1, self__C1pa, self__A2m1,
                 self__C2a, self__B21, self__C3a, self__A3c, self__B31, self__C4a,
                 self__A4, self__B41, self_s13, self_a13)
end

"""Private: General solution of position along geodesic"""
function _GenPosition(self::GeodesicLine, arcmode, s12_a12, outmask)
  s12_a12 = Float64(s12_a12)

  a12 = lat2 = lon2 = azi2 = s12 = m12 = M12 = M21 = S12 = Math.nan
  outmask &= self.caps & Geodesics.OUT_MASK
  if ! (arcmode ||
          (self.caps & (Geodesics.OUT_MASK & Geodesics.DISTANCE_IN)) > 0)
    # Uninitialized or impossible distance calculation requested
    return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12
  end

  # Avoid warning about uninitialized B12.
  B12 = 0.0
  AB1 = 0.0
  if arcmode
    # Interpret s12_a12 as spherical arc length
    sig12 = deg2rad(s12_a12)
    ssig12, csig12 = Math.sincosd(s12_a12)
  else
    # Interpret s12_a12 as distance
    tau12 = s12_a12 / (self._b * (1 + self._A1m1))
    s = sin(tau12)
    c = cos(tau12)
    # tau2 = tau1 + tau12
    B12 = - Geodesics._SinCosSeries(true,
                                    self._stau1 * c + self._ctau1 * s,
                                    self._ctau1 * c - self._stau1 * s,
                                    self._C1pa)
    sig12 = tau12 - (B12 - self._B11)
    ssig12 = sin(sig12)
    csig12 = cos(sig12)
    if abs(self.f) > 0.01
      # Reverted distance series is inaccurate for |f| > 1/100, so correct
      # sig12 with 1 Newton iteration.  The following table shows the
      # approximate maximum error for a = WGS_a() and various f relative to
      # GeodesicExact.
      #     erri = the error in the inverse solution (nm)
      #     errd = the error in the direct solution (series only) (nm)
      #     errda = the error in the direct solution (series + 1 Newton) (nm)
      #
      #       f     erri  errd errda
      #     -1/5    12e6 1.2e9  69e6
      #     -1/10  123e3  12e6 765e3
      #     -1/20   1110 108e3  7155
      #     -1/50  18.63 200.9 27.12
      #     -1/100 18.63 23.78 23.37
      #     -1/150 18.63 21.05 20.26
      #      1/150 22.35 24.73 25.83
      #      1/100 22.35 25.03 25.31
      #      1/50  29.80 231.9 30.44
      #      1/20   5376 146e3  10e3
      #      1/10  829e3  22e6 1.5e6
      #      1/5   157e6 3.8e9 280e6
      ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12
      csig2 = self._csig1 * csig12 - self._ssig1 * ssig12
      B12 = Geodesics._SinCosSeries(true, ssig2, csig2, self._C1a)
      serr = ((1 + self._A1m1) * (sig12 + (B12 - self._B11)) - s12_a12 / self._b)
      sig12 = sig12 - serr / sqrt(1 + self._k2 * Math.sq(ssig2))
      ssig12 = sin(sig12)
      csig12 = cos(sig12)
      # Update B12 below
    end
  end

  # real omg12, lam12, lon12
  # real ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2
  # sig2 = sig1 + sig12
  ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12
  csig2 = self._csig1 * csig12 - self._ssig1 * ssig12
  dn2 = sqrt(1 + self._k2 * Math.sq(ssig2))
  if outmask & (Geodesics.DISTANCE | Geodesics.REDUCEDLENGTH | Geodesics.GEODESICSCALE) > 0
    if arcmode || abs(self.f) > 0.01
      B12 = Geodesics._SinCosSeries(true, ssig2, csig2, self._C1a)
    end
    AB1 = (1 + self._A1m1) * (B12 - self._B11)
  end
  # sin(bet2) = cos(alp0) * sin(sig2)
  sbet2 = self._calp0 * ssig2
  # Alt: cbet2 = hypot(csig2, salp0 * ssig2)
  cbet2 = hypot(self._salp0, self._calp0 * csig2)
  if cbet2 == 0
    # I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
    cbet2 = csig2 = Geodesics.tiny_
  end
  # tan(alp0) = cos(sig2)*tan(alp2)
  salp2 = self._salp0
  calp2 = self._calp0 * csig2 # No need to normalize

  if (outmask & Geodesics.DISTANCE) > 0
    s12 = arcmode ? self._b * ((1 + self._A1m1) * sig12 + AB1) : s12_a12
  end

  if (outmask & Geodesics.LONGITUDE) > 0
    # tan(omg2) = sin(alp0) * tan(sig2)
    somg2 = self._salp0 * ssig2
    comg2 = csig2 # No need to normalize
    E = Math.copysign(1, self._salp0)          # East or west going?
    # omg12 = omg2 - omg1
    omg12 = ((outmask & Geodesics.LONG_UNROLL) > 0 ?
             E * (sig12
                  - (atan(          ssig2,       csig2) -
                     atan(    self._ssig1, self._csig1))
                  + (atan(E *       somg2,       comg2) -
                     atan(E * self._somg1, self._comg1))) :
             atan(somg2 * self._comg1 - comg2 * self._somg1,
                  comg2 * self._comg1 + somg2 * self._somg1))
    lam12 = omg12 + self._A3c * (
      sig12 + (Geodesics._SinCosSeries(true, ssig2, csig2, self._C3a) - self._B31))
    lon12 = rad2deg(lam12)
    lon2 = ((outmask & Geodesics.LONG_UNROLL) > 0 ?
            self.lon1 + lon12 :
            Math.AngNormalize(Math.AngNormalize(self.lon1) +
                              Math.AngNormalize(lon12)))
  end

  if (outmask & Geodesics.LATITUDE) > 0
    lat2 = Math.atan2d(sbet2, self._f1 * cbet2)
  end

  if (outmask & Geodesics.AZIMUTH) > 0
    azi2 = Math.atan2d(salp2, calp2)
  end

  if (outmask & (Geodesics.REDUCEDLENGTH | Geodesics.GEODESICSCALE)) > 0
    B22 = Geodesics._SinCosSeries(true, ssig2, csig2, self._C2a)
    AB2 = (1 + self._A2m1) * (B22 - self._B21)
    J12 = (self._A1m1 - self._A2m1) * sig12 + (AB1 - AB2)
    if (outmask & Geodesics.REDUCEDLENGTH) > 0
      # Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
      # accurate cancellation in the case of coincident points.
      m12 = self._b * ((      dn2 * (self._csig1 * ssig2) -
                        self._dn1 * (self._ssig1 * csig2))
                       - self._csig1 * csig2 * J12)
    end
    if (outmask & Geodesics.GEODESICSCALE) > 0
      t = (self._k2 * (ssig2 - self._ssig1) *
           (ssig2 + self._ssig1) / (self._dn1 + dn2))
      M12 = csig12 + (t * ssig2 - csig2 * J12) * self._ssig1 / self._dn1
      M21 = csig12 - (t * self._ssig1 - self._csig1 * J12) * ssig2 / dn2
    end
  end

  if (outmask & Geodesics.AREA) > 0
    B42 = Geodesics._SinCosSeries(false, ssig2, csig2, self._C4a)
    # real salp12, calp12
    if self._calp0 == 0 || self._salp0 == 0
      # alp12 = alp2 - alp1, used in atan2 so no need to normalize
      salp12 = salp2 * self.calp1 - calp2 * self.salp1
      calp12 = calp2 * self.calp1 + salp2 * self.salp1
    else
      # tan(alp) = tan(alp0) * sec(sig)
      # tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
      # = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
      # If csig12 > 0, write
      #   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
      # else
      #   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
      # No need to normalize
      salp12 = self._calp0 * self._salp0 * (
        csig12 <= 0 ?
        self._csig1 * (1 - csig12) + ssig12 * self._ssig1 :
        ssig12 * (self._csig1 * ssig12 / (1 + csig12) + self._ssig1))
      calp12 = (Math.sq(self._salp0) +
                Math.sq(self._calp0) * self._csig1 * csig2)
    end
    S12 = (self._c2 * atan(salp12, calp12) +
           self._A4 * (B42 - self._B41))
  end

  a12 = arcmode ? s12_a12 : rad2deg(sig12)
  a12, lat2, lon2, azi2, s12, m12, M12, M21, S12
end

"""
    Position(line::GeodesicLine, s12, outmask=STANDARD) -> result::Result

Find the position on the line given `s12`, the distance from the
first point to the second in metres.

The default value of `outmask` is `STANDARD`, i.e., the `lat1`,
`lon1`, `azi1`, `lat2`, `lon2`, `azi2`, `s12`, `a12` entries are
returned.  The `GeodesicLine` object must have been constructed with
the `DISTANCE_IN` capability.
"""
function Position(self::GeodesicLine, s12; outmask = GeodesicCapability.STANDARD)
  result = Result()
  result.lat1 = self.lat1
  result.lon1 = (outmask & Geodesics.LONG_UNROLL) > 0 ?
                 self.lon1 :
                 Math.AngNormalize(self.lon1)
  result.azi1 = self.azi1
  result.s12 = s12
  a12, lat2, lon2, azi2, s12, m12, M12, M21, S12 = _GenPosition(self, false, s12, outmask)
  outmask &= Geodesics.OUT_MASK
  result.a12 = a12
  (outmask & Geodesics.LATITUDE) > 0 && (result.lat2 = lat2)
  (outmask & Geodesics.LONGITUDE) > 0 && (result.lon2 = lon2)
  (outmask & Geodesics.AZIMUTH) > 0 && (result.azi2 = azi2)
  (outmask & Geodesics.REDUCEDLENGTH) > 0 && (result.m12 = m12)
  if (outmask & Geodesics.GEODESICSCALE) > 0
    result.M12 = M12
    result.M21 = M21
  end
  (outmask & Geodesics.AREA) > 0 && (result.S12 = S12)
  result
end

"""
    ArcPosition(line::GeodesicLine, a12, outmask=STANDARD) -> result::Result

Find the position on the line given `a12`, the spherical arc length from
the first point to the second in degrees.

The default value of `outmask` is `STANDARD`, i.e., the `lat1`,
`lon1`, `azi1`, `lat2`, `lon2`, `azi2`, `s12`, `a12` entries are
returned.
"""
function ArcPosition(self::GeodesicLine, a12; outmask = GeodesicCapability.STANDARD)
  result = Result()
  result.lat1 = self.lat1
  result.lon1 = (outmask & Geodesics.LONG_UNROLL) > 0 ?
                 self.lon1 :
                 Math.AngNormalize(self.lon1)
  result.azi1 = self.azi1
  result.a12 = a12
  a12, lat2, lon2, azi2, s12, m12, M12, M21, S12 = _GenPosition(self, true, a12, outmask)
  outmask &= Geodesics.OUT_MASK
  (outmask & Geodesics.DISTANCE) > 0  && (result.s12 = s12)
  (outmask & Geodesics.LATITUDE) > 0  && (result.lat2 = lat2)
  (outmask & Geodesics.LONGITUDE) > 0 && (result.lon2 = lon2)
  (outmask & Geodesics.AZIMUTH) > 0  && (result.azi2 = azi2)
  (outmask & Geodesics.REDUCEDLENGTH) > 0  && (result.m12 = m12)
  if (outmask & Geodesics.GEODESICSCALE) > 0
    result.M12 = M12
    result.M21 = M21
  end
  (outmask & Geodesics.AREA) > 0 && (result.S12 = S12)
  result
end

"""
    SetDistance(line::GeodesicLine, s13) -> line′::GeodesicLine

Specify the position of point 3 in terms of distance
from point 1 to point 3 in meters

Return a new `GeodesicLine` with `s13` and `a13` set.
"""
function SetDistance(self::GeodesicLine, s13)
  self = GeodesicLine(self.a, self.f, self._b, self._c2, self._f1, self.caps,
                      self.lat1, self.lon1, self.azi1, self.salp1, self.calp1,
                      self._dn1, self._salp0, self._calp0, self._ssig1, self._csig1,
                      self._somg1, self._comg1, self._k2, self._A1m1, self._C1a,
                      self._B11, self._stau1, self._ctau1, self._C1pa, self._A2m1,
                      self._C2a, self._B21, self._C3a, self._A3c, self._B31, self._C4a,
                      self._A4, self._B41,
                      s13, self.a13)
  a13, _, _, _, _, _, _, _, _ = _GenPosition(self, false, s13, 0)
  GeodesicLine(self.a, self.f, self._b, self._c2, self._f1, self.caps,
               self.lat1, self.lon1, self.azi1, self.salp1, self.calp1,
               self._dn1, self._salp0, self._calp0, self._ssig1, self._csig1,
               self._somg1, self._comg1, self._k2, self._A1m1, self._C1a,
               self._B11, self._stau1, self._ctau1, self._C1pa, self._A2m1,
               self._C2a, self._B21, self._C3a, self._A3c, self._B31, self._C4a,
               self._A4, self._B41,
               s13, a13)
end

"""
    SetArc(line::GeodesicLine, a13) -> line′::GeodesicLine

Specify the position of point 3 in terms of spherical arc length
from point 1 to point 3 in degrees.

Return a new `GeodesicLine` with `a13` and `s13` set.
"""
function SetArc(self::GeodesicLine, a13)
  self = GeodesicLine(self.a, self.f, self._b, self._c2, self._f1, self.caps,
                      self.lat1, self.lon1, self.azi1, self.salp1, self.calp1,
                      self._dn1, self._salp0, self._calp0, self._ssig1, self._csig1,
                      self._somg1, self._comg1, self._k2, self._A1m1, self._C1a,
                      self._B11, self._stau1, self._ctau1, self._C1pa, self._A2m1,
                      self._C2a, self._B21, self._C3a, self._A3c, self._B31, self._C4a,
                      self._A4, self._B41,
                      self.s13, a13)
  _, _, _, _, s13, _, _, _, _ = _GenPosition(self, true, a13, Geodesics.DISTANCE)
  GeodesicLine(self.a, self.f, self._b, self._c2, self._f1, self.caps,
               self.lat1, self.lon1, self.azi1, self.salp1, self.calp1,
               self._dn1, self._salp0, self._calp0, self._ssig1, self._csig1,
               self._somg1, self._comg1, self._k2, self._A1m1, self._C1a,
               self._B11, self._stau1, self._ctau1, self._C1pa, self._A2m1,
               self._C2a, self._B21, self._C3a, self._A3c, self._B31, self._C4a,
               self._A4, self._B41,
               s13, a13)
end

end # module
