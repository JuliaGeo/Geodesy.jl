module Polygons

import ..Math, ..Geodesics, ..Accumulators, .._GenDirect, ..WGS84

export Polygon, AddPoint!, AddEdge!, Compute

#"""Area of a geodesic polygon"""
mutable struct Polygon
	  """The geodesic object (readonly)"""
	  earth::Geodesics.Geodesic
	  """Is this a polyline? (readonly)"""
    polyline::Bool
	  """The total area of the ellipsoid in meter^2 (readonly)"""
    area0::Float64
	  _mask::Int
	  _areasum::Accumulators.Accumulator{Float64}
	  _perimetersum::Accumulators.Accumulator{Float64}
	  """The current number of points in the polygon (readonly)"""
	  num::Int
	  _lon0::Float64
	  _lat0::Float64
	  """The current latitude in degrees (readonly)"""
    lat1::Float64
    """The current longitude in degrees (readonly)"""
	  lon1::Float64
	  _crossings::Int
end

Base.:(==)(p1::Polygon, p2::Polygon) =
  all(x -> (isbits(x[1]) && isnan(x[1]) && isnan(x[2])) || x[1] == x[2],
      (getfield(p1, f), getfield(p2, f)) for f in fieldnames(Polygon))

# Treat the Geodesic struct as a scalar
Base.Broadcast.broadcastable(p::Polygon) = Ref(p)

"""
    Polygon([ellipsoid::Geodesic=WGS84,] polyline=false)

Construct a `Polygon`, which contains a set of points on a certain
`ellipsoid`.

With this construction, the `Polygon` contains no points.

If `polyline` is true, then the `Polygon` will not accumulate area
and instead only its perimeter can be calculated.
"""
function Polygon(ellipsoid::Geodesics.Geodesic, polyline::Bool=false)
	area0 = 4π*ellipsoid._c2
	_mask = (Geodesics.LATITUDE | Geodesics.LONGITUDE |
           Geodesics.DISTANCE | (polyline ? Geodesics.EMPTY :
                                 Geodesics.AREA | Geodesics.LONG_UNROLL))
	_areasum = Accumulators.Accumulator()
	_perimetersum = Accumulators.Accumulator()
	num = 0
	_lat0 = _lon0 = lat1 = lon1 = Math.nan
	_crossings = 0
	Polygon(ellipsoid, polyline, area0, _mask, _areasum, _perimetersum,
		num, _lon0, _lat0, lat1, lon1, _crossings)
end

Polygon(polyline::Bool=false) = Polygon(WGS84, polyline)

"""Count crossings of prime meridian for AddPoint."""
function _transit(lon1, lon2)
  # Return 1 or -1 if crossing prime meridian in east or west direction.
  # Otherwise return zero.
  # Compute lon12 the same way as Geodesic::Inverse.
  lon1 = Math.AngNormalize(lon1)
  lon2 = Math.AngNormalize(lon2)
  lon12, _ = Math.AngDiff(lon1, lon2)
  cross = (lon1 <= 0 && lon2 > 0 && lon12 > 0 ? 1 :
           (lon2 <= 0 && lon1 > 0 && lon12 < 0 ? -1 : 0))
  cross
end

"""Count crossings of prime meridian for AddEdge."""
function _transitdirect(lon1, lon2)
  # We want to compute exactly
  #   int(floor(lon2 / 360)) - int(floor(lon1 / 360))
  # Since we only need the parity of the result we can use std::remquo but
  # this is buggy with g++ 4.8.3 and requires C++11.  So instead we do
  lon1 = lon1 % 720.0
  lon2 = lon2 % 720.0
  (((lon2 >= 0 && lon2 < 360) || lon2 < -360) ? 0 : 1) -
         (((lon1 >= 0 && lon1 < 360) || lon1 < -360) ? 0 : 1)
end

"""Reset to empty polygon."""
function Clear!(self)
    self.num = 0
    self._crossings = 0
    if !self.polyline
    	Accumulators.Set!(self._areasum, 0)
    end
    Accumulators.Set!(self._perimetersum, 0)
    self._lat0 = self._lon0 = self.lat1 = self.lon1 = Math.nan
    self
end

"""
    AddPoint!(polygon, lat, lon) -> polygon

Add the next vertex to the polygon

`lat`: the latitude of the point in degrees
`lon`: the longitude of the point in degrees

This adds an edge from the current vertex to the new vertex.
"""
function AddPoint!(self::Polygon, lat, lon)
  if self.num == 0
    self._lat0 = self.lat1 = lat
    self._lon0 = self.lon1 = lon
  else
    _, s12, _, _, _, _, _, _, _, S12 = Geodesics._GenInverse(self.earth,
      self.lat1, self.lon1, lat, lon, self._mask)
    Accumulators.Add!(self._perimetersum, s12)
    if !self.polyline
      Accumulators.Add!(self._areasum, S12)
      self._crossings += _transit(self.lon1, lon)
	  end
    self.lat1 = lat
    self.lon1 = lon
  end
  self.num += 1
  self
end

"""
    AddEdge!(polygon, azimuth, distance) -> polygon

Add the next edge to the polygon

`azi`: the azimuth at the current the point in degrees
`s`: the length of the edge in meters

This specifies the new vertex in terms of the edge from the current
vertex.
"""
function AddEdge!(self::Polygon, azi, s)
  if self.num != 0
    _, lat, lon, _, _, _, _, _, S12 = _GenDirect(self.earth,
      self.lat1, self.lon1, azi, false, s, self._mask)
    Accumulators.Add!(self._perimetersum, s)
    if !self.polyline
      Accumulators.Add!(self._areasum, S12)
      self._crossings += _transitdirect(self.lon1, lon)
    end
    self.lat1 = lat
    self.lon1 = lon
    self.num += 1
  else
    error("cannot add an edge to a polygon with no points")
  end
  self
end

# return number, perimeter, area
"""
    Compute(polygon, reverse=false, sign=true) -> n, perimeter, area

Compute the properties of the polygon

- `reverse`: if true then clockwise (instead of
  counter-clockwise) traversal counts as a positive area
- `sign`: if true then return a signed result for the area if the
  polygon is traversed in the "wrong" direction instead of returning
  the area for the rest of the ellispoid

Return a tuple of number, perimeter (meters), area (meters^2).

If the object is a polygon (and not a polygon), the perimeter
includes the length of a final edge connecting the current point to
the initial point.  If the object is a polyline, then area is nan.

More points can be added to the polygon after this call.

"""
function Compute(self::Polygon, reverse = false, sign = true)
  if self.polyline
  	area = Math.nan
  end
  if self.num < 2
    perimeter = 0.0
    if !self.polyline
      area = 0.0
    end
    return self.num, perimeter, area
  end

  if self.polyline
    perimeter = Accumulators.Sum(self._perimetersum)
    return self.num, perimeter, area
  end

  _, s12, _, _, _, _, _, _, _, S12 = Geodesics._GenInverse(self.earth,
    self.lat1, self.lon1, self._lat0, self._lon0, self._mask)
  perimeter = Accumulators.Sum(self._perimetersum, s12)
  tempsum = Accumulators.Accumulator(self._areasum)
  Accumulators.Add!(tempsum, S12)
  crossings = self._crossings + _transit(self.lon1, self._lon0)
  if (crossings & 1) > 0
    Accumulators.Add!(tempsum, (Accumulators.Sum(tempsum) < 0 ? 1 : -1) * self.area0/2)
  end
  # area is with the clockwise sense.  If !reverse convert to
  # counter-clockwise convention.
  if !reverse
  	Accumulators.Negate!(tempsum)
  end
  # If sign put area in (-area0/2, area0/2], else put area in [0, area0)
  if sign
    if Accumulators.Sum(tempsum) > self.area0/2
      Accumulators.Add!(tempsum, -self.area0 )
    elseif Accumulators.Sum(tempsum) <= -self.area0/2
      Accumulators.Add!(tempsum,  self.area0 )
    end
  else
    if Accumulators.Sum(tempsum) >= self.area0
      Accumulators.Add!(tempsum, -self.area0 )
    elseif Accumulators.Sum(tempsum) < 0
      Accumulators.Add!(tempsum,  self.area0 )
    end
  end

  area = 0.0 + Accumulators.Sum(tempsum)
  self.num, perimeter, area
end


"""Compute the properties for a tentative additional vertex

:param lat: the latitude of the point in degrees
:param lon: the longitude of the point in degrees
:param reverse: if true then clockwise (instead of
  counter-clockwise) traversal counts as a positive area
:param sign: if true then return a signed result for the area if the
  polygon is traversed in the "wrong" direction instead of returning
  the area for the rest of the earth
:return: a tuple of number, perimeter (meters), area (meters^2)

"""
function TestPoint(self::Polygon, lat, lon, reverse = false, sign = true)
  if self.polyline
  	area = Math.nan
  end
  if self.num == 0
    perimeter = 0.0
    if !self.polyline
      area = 0.0
    end
    return 1, perimeter, area
  end

  perimeter = Accumulators.Sum(self._perimetersum)
  tempsum = self.polyline ? 0.0 : Accumulators.Sum(self._areasum)
  crossings = self._crossings
  num = self.num + 1
  for i in (self.polyline ? (0,) : (0, 1))
    _, s12, _, _, _, _, _, _, _, S12 = Geodesics._GenInverse(self.earth,
      i == 0 ? self.lat1 : lat,
      i == 0 ? self.lon1 : lon,
      i != 0 ? self._lat0 : lat,
      i != 0 ? self._lon0 : lon,
      self._mask)
    perimeter += s12
    if !self.polyline
      tempsum += S12
      crossings += _transit(i == 0 ? self.lon1 : lon,
                            i != 0 ? self._lon0 : lon)
    end
  end

  if self.polyline
    return num, perimeter, area
  end

  if (crossings & 1) > 0
    tempsum += (tempsum < 0 ? 1 : -1) * self.area0/2
  end
  # area is with the clockwise sense.  If !reverse convert to
  # counter-clockwise convention.
  if !reverse
  	tempsum *= -1
  end
  # If sign put area in (-area0/2, area0/2], else put area in [0, area0)
  if sign
    if tempsum > self.area0/2
      tempsum -= self.area0
    elseif tempsum <= -self.area0/2
      tempsum += self.area0
    end
  else
    if tempsum >= self.area0
      tempsum -= self.area0
    elseif tempsum < 0
      tempsum += self.area0
    end
  end

  area = 0.0 + tempsum
  num, perimeter, area
end

# return num, perimeter, area
"""Compute the properties for a tentative additional edge

:param azi: the azimuth at the current the point in degrees
:param s: the length of the edge in meters
:param reverse: if true then clockwise (instead of
  counter-clockwise) traversal counts as a positive area
:param sign: if true then return a signed result for the area if the
  polygon is traversed in the "wrong" direction instead of returning
  the area for the rest of the earth
:return: a tuple of number, perimeter (meters), area (meters^2)

"""
function TestEdge(self::Polygon, azi, s, reverse = False, sign = True)
  if self.num == 0           # we don't have a starting point!
    return 0, Math.nan, Math.nan
  end
  num = self.num + 1
  perimeter = self._perimetersum.Sum() + s
  if self.polyline
    return num, perimeter, Math.nan
  end

  tempsum =  Accumulators.Sum(self._areasum)
  crossings = self._crossings
  _, lat, lon, _, _, _, _, _, S12 = _GenDirect(self.earth,
    self.lat1, self.lon1, azi, false, s, self._mask)
  tempsum += S12
  crossings += _transitdirect(self.lon1, lon)
  _, s12, _, _, _, _, _, _, _, S12 = Geodesics._GenInverse(self.earth,
    lat, lon, self._lat0, self._lon0, self._mask)
  perimeter += s12
  tempsum += S12
  crossings += _transit(lon, self._lon0)

  if (crossings & 1) > 0
    tempsum += (tempsum < 0 ? 1 : -1) * self.area0/2
  end
  # area is with the clockwise sense.  If !reverse convert to
  # counter-clockwise convention.
  if !reverse
  	tempsum *= -1
  end
  # If sign put area in (-area0/2, area0/2], else put area in [0, area0)
  if sign
    if tempsum > self.area0/2
      tempsum -= self.area0
    elseif tempsum <= -self.area0/2
      tempsum += self.area0
    end
  else
    if tempsum >= self.area0
      tempsum -= self.area0
    elseif tempsum < 0
      tempsum += self.area0
    end
  end

  area = 0.0 + tempsum
  num, perimeter, area
end

end # module
