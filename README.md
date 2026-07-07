# Geodesy

[![Build Status](https://github.com/JuliaGeo/Geodesy.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeo/Geodesy.jl/actions)

**Geodesy** is a Julia package for working with points in various world and
local coordinate systems. The primary feature of *Geodesy* is to define and
perform coordinate transformations in a convenient and safe framework,
leveraging the *CoordinateTransformations* [package](https://github.com/FugroRoames/CoordinateTransformations.jl).
Transformations are accurate and efficient and implemented in native Julia code
(with many functions being ported from Charles Karney's *GeographicLib*
[C++ library](http://geographiclib.sourceforge.net/)), and some common geodetic
datums are provided for convenience.

*This package is currently in maintenance mode and new features are not being developed. We endeavour to fix any bugs that are reported.*

## Quick start

Lets define a 3D point by its latitude, longitude and altitude (LLA):
```julia
x_lla = LLA(-27.468937, 153.023628, 0.0) # City Hall, Brisbane, Australia
```

This can be converted to a Cartesian Earth-Centered-Earth-Fixed (ECEF)
coordinate simply by calling the constructor
```julia
x_ecef = ECEF(x_lla, wgs84)
```

Here we have used the WGS-84 ellipsoid to calculate the transformation, but other
datums such as `osgb36`, `nad27` and `grs80` are provided. All transformations
use the *CoordinateTransformations*' interface, and the above is short for
```julia
x_ecef = ECEFfromLLA(wgs84)(x_lla)
```
where `ECEFfromLLA` is a type inheriting from *CoordinateTransformations*'
`Transformation`. (Similar names `XfromY` exist for each of the
coordinate types.)

Often, points are measured or required in a *local* frame, such as the north-east-up
coordinates with respect to a given origin. The `ENU` type represents points in this
coordinate system and we may transform between ENU and globally referenced
coordinates using `ENUfromLLA`, etc.
```julia
origin_lla = LLA(-27.468937, 153.023628, 0.0) # City Hall, Brisbane, Australia
point_lla = LLA(-27.465933, 153.025900, 0.0)  # Central Station, Brisbane, Australia

# Define the transformation and execute it
trans = ENUfromLLA(origin_lla, wgs84)
point_enu = trans(point_lla)

# Equivalently
point_enu = ENU(point_lla, origin_lla, wgs84)
```
The `NED` coordinate system is also available as an alternative to `ENU`, with 
the same interface.

Similarly, we could convert to UTM/UPS coordinates, and two types are provided
for this - `UTM` stores 3D coordinates `x`, `y`, and `z` in an unspecified zone,
while `UTMZ` includes the `zone` number and `hemisphere` bool (where `true` =
northern, `false` = southern). To get the canonical zone for your coordinates,
simply use:
```julia
x_utmz = UTMZ(x_lla, wgs84)
```

If you are transforming a large number of points to or from a given zone, it may
be more effective to define the transformation explicitly and use the lighter
`UTM` storage type.
```julia
points_lla::Vector{LLA{Float64}}
utm_from_lla = UTMfromLLA(56, false, wgs84) # Zone 56-South
points_utm = map(utm_from_lla, points_lla) # A new vector of UTM coordinates
```

*Geodesy* becomes particularly powerful when you chain together transformations.
For example, you can define a single transformation from your data on disk in UTM
coordinates to a local frame in ENU coordinates. Internally, this will perform
UTM (+ zone) → LLA → ECEF → ENU via composing transformations with `∘` into a
`ComposedTransformation`:
```julia
julia> origin = LLA(-27.468937, 153.023628, 0.0) # City Hall, Brisbane, Australia
LLA(lat=-27.468937°, lon=153.023628°, alt=0.0)

julia> trans = ENUfromUTMZ(origin, wgs84)
(ENUfromECEF(ECEF(-5.046925124630393e6, 2.5689157252069353e6, -2.924416653602336e6), lat=-27.468937°, lon=153.023628°) ∘ (ECEFfromLLA(wgs84) ∘ LLAfromUTMZ(wgs84)))
```

This transformation can then be composed with rotations and translations in
*CoordinateTransformations* (or your own custom-defined `AbstractTransformation`
to define further reference frames. For example, in this way, a point measured
by a scanner on a moving vehicle at a particular time may be globally
georeferenced with a single call to the `Transformation`!

Finally, the Cartesian distance between world points can be calculated via
automatic transformation to a Cartesian frame:
```julia
x_lla = LLA(-27.468937, 153.023628, 0.0) # City Hall, Brisbane, Australia
y_lla = LLA(-27.465933, 153.025900, 0.0) # Central Station, Brisbane, Australia
euclidean_distance(x_lla, y_lla)                   # 401.54 meters
```
(assuming the `wgs84` datum, which can be configured in `distance(x, y, datum)`).


## Alternatives
- [CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl) - Native Julia implementation with built-in support for units, widely tested against the PROJ C library.
- [MapMaths.jl](https://github.com/subnero1/MapMaths.jl) - Native Julia implementation with support for the WGS84 datum and local coordinates.
- [Proj.jl](https://github.com/JuliaGeo/Proj.jl) - Julia wrapper for the PROJ C library, which provides a wide range of coordinate system transformations and projections.


## Basic Terminology

This section describes some terminology and concepts that are relevant to
*Geodesy.jl*, attempting to define Geodesy-specific jargon where possible.  For
a longer, less technical discussion with more historical context, ICSM's
[Fundamentals of Mapping page](https://www.icsm.gov.au/education/fundamentals-mapping)
is highly recommended.

### Coordinate Reference Systems and Spatial Reference Identifiers

A position on the Earth can be given by some numerical coordinate values, but
those don't mean much without more information.  The extra information is called
the **Coordinate Reference System** or **CRS** (also known as a *Spatial
Reference System* or SRS).  A CRS tells you two main things:

* The measurement procedure: which real world objects were used to
  define the frame of reference or *datum* of the measurement?
* The *coordinate system*: how do coordinate numerical values relate to the
  reference frame defined by the datum?

The full specification of a CRS can be complex, so a short label called a
**Spatial Reference IDentifier** or **SRID** is usually used instead.  For
example, [EPSG:4326](http://epsg.io/4326) is one way to refer to the 2D WGS84
latitude and longitude you'd get from a mobile phone GPS device.  An SRID
is of the form `AUTHORITY:CODE`, where the code is a number and the authority is
the name of an organization maintaining a list of codes with associated CRS
information.  There are services where you can look up a CRS, for example,
<http://epsg.io> is a convenient interface to the SRIDs maintained by the
*European Petroleum Survey Group* (EPSG) authority.  Likewise,
<http://spatialreference.org> is an open registry to which anyone can
contribute.

When maintaining a spatial database, it's typical to define an internal list of
SRIDs (effectively making your organization the authority), and a mapping from
these to CRS information.  A link back to a definitive SRID from an external
authority should also be included where possible.

### Datums

In spatial measurement and positioning, a **datum** is a set of reference
objects with given coordinates, *relative to which* other objects may be
positioned.  For example, in traditional surveying a datum might comprise a
pair of pegs in the ground, separated by a carefully measured distance.  When
surveying the position of an unknown but nearby point, the angle back to the
original datum objects can be measured using a theodolite.  After this, the
relative position of the new point can be computed using simple triangulation.
Repeating this trick with any of the now three known points, an entire
triangulation network of surveyed objects can be extended outward.  Any point
surveyed relative to the network is said to be measured *in the datum* of the
original objects.  Datums are often named with an acronym, for example OSGB36 is
the Ordnance Survey of Great Britain, 1936.

In the era of satellite geodesy, coordinates are determined for an object
by timing signals from a satellite constellation (eg, the GPS satellites) and
computing position relative to those satellites.  Where is the datum here? At
first glance the situation seems quite different from the traditional setup
described above.  However, the satellite positions as a function of time
(*ephemerides*, in the jargon) must themselves be defined relative to some
frame. This is done by continuously observing the satellites from a set of
highly stable ground stations equipped with GPS receivers. It is the full set of
these ground stations and their assigned coordinates which form the datum.

Let's inspect the flow of positional information in both cases:
* For traditional surveying,
  ```
  datum object positions -> triangulation network -> newly surveyed point
  ```
* For satellite geodesy,
  ```
  datum object positions -> satellite ephemerides -> newly surveyed point
  ```

We see that the basic nature of a datum is precisely the same regardless of
whether we're doing a traditional survey or using a GPS receiver.


### Terrestrial reference systems and frames

Coordinates for new points are measured by transferring coordinates from the
datum objects, as described above.  However, how do we decide
on coordinates for the datum objects themselves?  This is purely a matter of
convention, consistency and measurement.

For example, the **International Terrestrial Reference System** (**ITRS**) is a
reference system that rotates with the Earth so that the average velocity of
the crust is zero. That is, in this reference system the only crust movement is
geophysical.  Roughly speaking, the *defining conventions* for the ITRS are:

* Space is modeled as a three-dimensional Euclidean affine space.
* The origin is at the center of mass of the Earth (it is *geocentric*).
* The z-axis is the axis of rotation of the Earth.
* The scale is set to 1 SI meter.
* The x-axis is orthogonal to the z-axis and aligns with the international
  reference meridian through Greenwich.
* The y-axis is set to the cross product of the z and x axes, forming a right
  handed coordinate frame.
* Various rates of change of the above must also be specified, for example, the
  scale should stay constant in time.

The precise conventions are defined in chapter 4 of the
[IERS conventions](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)
published by the International Earth Rotation and Reference Service (IERS).
These conventions define an ideal reference *system*, but they're useless
without physical measurements that give coordinates for a set of real world
datum objects.  The process of measuring and computing coordinates for datum
objects is called *realizing* the reference system and the result is called a
*reference frame*.  For example, the **International Terrestrial Reference Frame
of 2014** (**ITRF2014**) realizes the ITRS conventions using raw measurement
data gathered in the 25 years prior to 2014.

To measure and compute coordinates, several space geodesy techniques are used to
gather raw measurement data; currently the IERS includes
[VLBI (very long baseline interferometry)](https://en.wikipedia.org/wiki/Very-long-baseline_interferometry) of distant astronomical radio sources,
[SLR (satellite laser ranging)](https://en.wikipedia.org/wiki/Satellite_laser_ranging),
[GPS (global positioning system)](https://en.wikipedia.org/wiki/Global_Positioning_System) and
[DORIS (gosh these acronyms are tiring)](https://en.wikipedia.org/wiki/DORIS_(geodesy)).
The raw data is not in the form of positions, but must be condensed down in a
large scale fitting problem, ideally by requiring physical and statistical
consistency of all measurements, tying measurements at different sites together
with physical models.


### Coordinate systems

In geometry, a **coordinate system**
[is a system](https://en.wikipedia.org/wiki/Coordinate_system)
which uses one or more numbers, or **coordinates** to uniquely
determine the position of a point in a mathematical space such as Euclidean
space.  For example, in geodesy a point is commonly referred to using geodetic
latitude, longitude and height relative to a given reference ellipsoid; this is
called a **geodetic coordinate system**.

An [**ellipsoid**](https://en.wikipedia.org/wiki/Ellipsoid) is chosen because
it's a reasonable model for the shape of the Earth and its gravitational field
without being overly complex; it has only a few parameters, and a simple
mathematical form.  The term [**spheroid**](https://en.wikipedia.org/wiki/Spheroid)
is also used because the ellipsoids in use today are rotationally symmetric
around the pole. Note that there's several ways to define
[latitude](https://en.wikipedia.org/wiki/Latitude) on an ellipsoid. The most
natural for geodesy is **geodetic latitude**, used by default because it's
physically accessible in any location as a good approximation to the angle
between the gravity vector and the equatorial plane.  (This type of latitude
is not an angle measured at the centre of the ellipsoid, which may be surprising
if you're used to spherical coordinates!)

There are usually several useful coordinate systems for the same space.  As well
as the geodetic coordinates mentioned above, it's common to see

* The x,y,z components in an Earth-Centred Cartesian coordinate system rotating
  with the Earth.  This is conventionally called an
  **Earth-Centred Earth-Fixed** (**ECEF**) coordinate system. This is a natural
  coordinate system in which to define coordinates for the datum objects
  defining a terrestrial reference frame.
* The east,north and up **ENU** components of a Cartesian coordinate frame at a
  particular point on the ellipsoid.  This coordinate system is useful as a
  local frame for navigation. It is also common for navigation systems to use
  north,east and down (**NED**) coordinates.
* Easting,northing and vertical components of a **projected coordinate system** or
  [**map projection**](http://www.icsm.gov.au/mapping/about_projections.html).
  There's an entire zoo of these, designed to represent the curved surface of an
  ellipsoid with a flat map.

Different coordinates systems provide different coordinates for the same point,
so it's obviously important to specify exactly which coordinate system you're
using.  In particular, you should specify which ellipsoid parameters are in
use if you deal with latitude and longitude, as in principle you could have more
than one ellipsoid.  This is a point of confusion, because a datum in geodesy
also comes with a reference ellipsoid as a very strong matter of convention
(thus being called a **geodetic datum**).

With its conventional ellipsoid, a geodetic datum also defines a conventional
geodetic coordinate system, thus bringing together concepts which are
interconnected but conceptually distinct.  To emphasize:

* A coordinate system is a mathematical abstraction allowing us to manipulate
  *geometric* quantities using numeric and algebraic techniques.  By itself,
  mathematical geometry is pure abstraction without a connection to the physical
  world.
* A datum is a set of physical objects with associated coordinates, thereby
  *defining* a reference frame in a way which is physically accessible.  A datum
  is the bridge which connects physical reality to the abstract ideal of
  mathematical geometry, via the algebraic mechanism of a coordinate system.


## The API

### Coordinate types

Geodesy provides several in-built coordinate storage types for convenience and
safety. The philosophy is to avoid carrying around raw data in generic containers
like `Vector`s with no concept of what coordinate system it is in.

##### `LLA{T}` - latitude, longitude and altitude

The global `LLA` type stores data in a lat-lon-alt order, where latitude and longitude
are expected in degrees (not radians). A keyword constructor, `LLA(lat=x, lon=y, alt=z)`,
is also provided to help with having to remember the storage order.

##### `LatLon{T}` - latitude and longitude

The 2D `LatLon` type stores data in a lat-lon order, where latitude and longitude
are expected in degrees (not radians). A keyword constructor, `LatLon(lat=x, lon=y)`,
is also provided. `LatLon` is currently the only supported 2D coordinate.

##### `ECEF{T}` - Earth-centered, Earth-fixed

The global `ECEF` type stores Cartesian coordinates `x`, `y`, `z`, according to the
[usual convention](https://en.wikipedia.org/wiki/ECEF). Being a Cartesian frame,
`ECEF` is a subtype of [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl)'
`StaticVector` and they can be added and subtracted with themselves and other
vectors.

##### `UTM{T}` - universal transverse-Mercator

The `UTM` type encodes the easting `x`, northing `y` and height `z` of a UTM
coordinate in an unspecified zone. This data type is also used to encode
universal polar-stereographic (UPS) coordinates (where the zone is `0`).

##### `UTMZ{T}` - universal transverse-Mercator + zone

In addition to the easting `x`, northing `y` and height `z`, the global `UTMZ` type
also encodes the UTM `zone` and `hemisphere`, where `zone` is a `UInt8` and
`hemisphere` is a `Bool` for compact storage. The northern hemisphere is
denoted as `true`, and the southern as `false`. Zone `0` corresponds to the UPS
projection about the corresponding pole, otherwise `zone` is an integer between
`1` and `60`.

##### `ENU{T}` - east-north-up

The `ENU` type stores local Cartesian coordinates that encode a point's distance
towards east `e`, towards north `n` and upwards `u` with respect to an
unspecified origin. Like `ECEF`, `ENU` is also a subtype of `StaticVector`.

##### `NED{T}` - north-east-down

The `NED` type is an alternative convention to `ENU`. It stores local Cartesian 
coordinates that encode a point's distance towards north `n`, towards east `e` 
and downwards `d` with respect to an unspecified origin. Like `ECEF`, `NED` is 
also a subtype of `StaticVector`.

### Geodetic Datums

Geodetic datums are modelled as subtypes of the abstract type `Datum`.  The
associated ellipsoid may be obtained by calling the `ellipsoid()` function, for
example, `ellipsoid(NAD83())`.

There are several pre-defined datums.  Worldwide datums include

* `WGS84` - standard GPS datum for moderate precision work (representing both
  the latest frame realization, or if time is supplied a discontinuous dynamic
  datum where time looks up the frame implementation date in the broadcast
  ephemerides.)
* `WGS84{GpsWeek}` - specific realizations of the WGS84 frame.
* `ITRF{Year}` - Realizations of the International Terrestrial Reference System
  for high precision surveying.

National datums include

* `OSGB36` - Ordnance Survey of Great Britain of 1936.
* `NAD27`, `NAD83` - North American Datums of 1927 and 1983, respectively
* `GDA94` - Geocentric Datum of Australia, 1994.

Datums may also be passed to coordinate transformation constructors such as
transverse-Mercator and polar-stereographic projections in which case the
associated ellipsoid will be extracted to form the transformation.  For datums
without extra parameters (everything except `ITRF` and `WGS84{Week}`) there is a
standard instance defined to reduce the amount of brackets you have to type.
For example, `LLAfromECEF(NAD83())` and `LLAfromECEF(nad83)` are equivalent.

### Transformations and conversions

*Geodesy* provides two interfaces changing coordinate systems.

"Transformations" are based on *CoordinateTransformations* interface for defining
`AbstractTransformation`s and allow the user to apply them by calling them,
invert them with `inv()` and compose them with `compose()` or `∘`. The transformations
cache any possible pre-calculations for efficiency when the same transformation
is applied to many points.

"Conversions" are based on type-constructors, obeying simple syntax like `LLA(ecef, datum)`.
The `datum` or other information is *always* necessary, as no assumptions are
made by *Geodesy* for safety and consistency reasons. Similarly, `Base.convert`
is not defined because, without assumptions, it would require additional
information. The main drawback of this approach is that some calculations may not
be pre-cached (for instance, the origin of an ENU transformation).

#### Between `LLA` and `ECEF`

The `LLAfromECEF` and `ECEFfromLLA` transformations require an ellipsoidal datum
to perform the conversion. The exact transformation is performed in both directions,
using a port of the ECEF → LLA transformation from *GeographicLib*.

Note that in some cases where points are very close to the centre of the ellipsoid,
multiple equivalent `LLA` points are valid solutions to the transformation problem.
Here, as in *GeographicLib*, the point with the greatest altitude is chosen.

#### Between `LLA` and `UTM`/`UTMZ`

The `LLAfromUTM(Z)` and `UTM(Z)fromLLA` transformations also require an
ellipsoidal datum to perform the conversion. The transformation retains a cache
of the parameters used in the transformation, which in the case of the
transverse-Mercator projection leads to a significant saving.

In all cases zone `0` corresponds to the UPS coordinate system, and the
polar-stereographic projection of *GeographicLib* has been ported to Julia to
perform the transformation.

An approximate, 6th-order expansion is used by default for the transverse-Mercator
projection and its inverse (though orders 4-8 are defined). The algorithm is a
native Julia port of that used in *GeographicLib*, and is accurate to nanometers
for up to several UTM zones away from the reference meridian. However, the
series expansion diverges at ±90° from the reference meridian. While the `UTMZ`-methods
will automatically choose the canonical zone and hemisphere for the input,
extreme care must be taken to choose an appropriate zone for the `UTM` methods
(In the future, we will implement the exact UTM transformation as a fallback —
contributions welcome!).

There is also `UTMfromUTMZ` and `UTMZfromUTM` transformations that are helpful
for converting between these two formats and putting data into the same `UTM`
zone.

#### To and from local `ENU` and `NED` frames

The `ECEFfromENU` and `ENUfromECEF` transformations define the transformation
around a specific origin. Both the origin coordinates as an `ECEF` as well as
its corresponding latitude and longitude are stored in the transformation for
maximal efficiency when performing multiple `transform`s. The transformation can
be inverted with `inv` to perform the reverse transformation with respect to the
same origin. Moreover, the `ENUfromNED` transformation and its inverse `NEDfromENU` 
can be used to change a local frame's orientation convention between `ENU` and `NED`, 
while keeping the same origin.

#### Web Mercator support

We support the Web Mercator / Pseudo Mercator projection with the
`WebMercatorfromLLA` and `LLAfromWebMercator` transformations for
interoperability with many web mapping systems. The scaling of the
northing and easting is defined to be meters at the Equator, the same as how
proj handles this (see https://proj.org/operations/projections/webmerc.html ).

If you need to deal with web mapping tile coordinate systems (zoom levels and
pixel coordinates, etc) these could be added by composing another
transformation on top of the web mercator projection defined in this package.

#### Composed transformations

Many other methods are defined as convenience constructors for composed
transformations, to go between any two of the coordinate types defined here.
These include:

* `ECEFfromUTMZ(datum) = ECEFfromLLA(datum) ∘ LLAfromUTMZ(datum)`
* `UTMZfromECEF(datum) = UTMZfromLLA(datum) ∘ LLAfromECEF(datum)`
* `UTMfromECEF(zone, hemisphere, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromECEF(datum)`
* `ECEFfromUTM(zone, hemisphere, datum) = ECEFfromLLA(datum) ∘ LLAfromUTM(zone, hemisphere, datum)`
* `ENUfromLLA(origin, datum) = ENUfromECEF(origin, datum) ∘ ECEFfromLLA(datum)`
* `LLAfromENU(origin, datum) = LLAfromECEF(datum) ∘ ECEFfromENU(origin, datum)`
* `NEDfromLLA(origin, datum) = NEDfromENU() ∘ ENUfromECEF(origin, datum) ∘ ECEFfromLLA(datum)`
* `LLAfromNED(origin, datum) = LLAfromECEF(datum) ∘ ECEFfromENU(origin,datum) ∘ ENUfromNED()`
* `ECEFfromUTMZ(datum) = ECEFfromLLA(datum) ∘ LLAfromUTMZ(datum)`
* `ENUfromUTMZ(origin, datum)  = ENUfromLLA(origin, datum) ∘ LLAfromUTMZ(datum`
* `UTMZfromENU(origin, datum)  = UTMZfromLLA(datum) ∘ LLAfromENU(origin, datum)`
* `NEDfromUTMZ(origin, datum)  = NEDfromLLA(origin, datum) ∘ LLAfromUTMZ(datum`
* `UTMZfromNED(origin, datum)  = UTMZfromLLA(datum) ∘ LLAfromNED(origin, datum)`
* `UTMfromENU(origin, zone, hemisphere, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromENU(origin, datum)`
* `ENUfromUTM(origin, zone, hemisphere, datum) = ENUfromLLA(origin, datum) ∘ LLAfromUTM(zone, hemisphere, datum)`
* `UTMfromNED(origin, zone, hemisphere, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromNED(origin, datum)`
* `NEDfromUTM(origin, zone, hemisphere, datum) = NEDfromLLA(origin, datum) ∘ LLAfromUTM(zone, hemisphere, datum)`

Constructor-based transforms for these are also provided, such as `UTMZ(ecef, datum)`
which converts to `LLA` as an intermediary, as above. When converting multiple
points to or from the *same* ENU reference frame, it is recommended to use the
transformation-based approach for efficiency. However, the other
constructor-based conversions should be similar in speed to their transformation
counterparts.

### Distance

Currently, the only defined distance measure is the straight-line or [Euclidean
distance](https://en.wikipedia.org/wiki/Euclidean_distance),
`euclidean_distance(x, y, [datum = wgs84])`, which works for all combinations
of types for `x` and `y` - except that the UTM zone and hemisphere must also be
provided for `UTM` types, as in `euclidean_distance(utm1, utm2, zone,
hemisphere, [datum = wgs84])` (the Cartesian distance for `UTM` types is not
approximated, but achieved via conversion to `ECEF`).

This is the only function currently in *Geodesy* which takes a default datum,
and *should* be relatively accurate for close points where Euclidean distances
are most important. Future work may focus on geodesics and related calculations
(contributions welcome!).
