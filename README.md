# Geodesy

[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

**Geodesy** is a Julia package for working with points in various world and
local coordinate systems. The primary feature of *Geodesy* is to define and
perform coordinate transformations in a convenient and safe framework,
leveraging the *CoordinateTransformations* [package](https://github.com/FugroRoames/CoordinateTransformations.jl).
Transformations are accurate and efficient and implemented in native Julia code
(with many functions being ported from Charles Karney's *GeographicLib*
[C++ library](http://geographiclib.sourceforge.net/)), and some common geodetic
datums are provided for convenience.

### Quick start

Lets define a 3D point by it's latitude, longitude and altitude (LLA):
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
point_enu = ENU(point_enu, point_origin, wgs84)
```

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
distance(x_lla, y_lla)                   # 401.54 meters
```
(assuming the `wgs84` datum, which can be configured in `distance(x, y, datum)`).

### Coordinate types

Geodesy provides several in-built coordinate storage types for convenience and
safety. The philosophy is to avoid carrying around raw data in `Vector`s,
*FixedSizeArray* `Vec`s or `Point`s with no concept of what coordinate system it
is in.

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
[usual convention](https://en.wikipedia.org/wiki/ECEF).

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

The `ENU` type is a local Cartesian coordinate that encodes a point's distance
towards east `e`, towards north `n` and upwards `u` with respect to an
unspecified origin.

### Datums

Geodesy comes with several in-built geodetic (i.e. ellipsoidal) datums.

* `wgs84`
* `osgb36`
* `nad27`
* `grs80`

These are instances of Julia singleton types (`WGS84`, etc) and are used to
overload internal methods for getting the elliptical parameters (as a `Geodesy.Ellipsoid`)
and pre-cached calculations for transverse-Mercator and polar-stereographic
projections (as a `Geodesy.TransverseMercator`).

To define your own ellipsoidal datum, it would be most efficient to define your own
singleton types and add methods to `Geodesy.ellipsoid` and `Geodesy.TransverseMercator`
to refer to pre-computed data `Ellipsoid`s and `TransverseMercator`s.

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
using a port the ECEF → LLA transformation from *GeographicLib*.

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
extreme care must be taken to choose an appropriate zone for the `UTM` methods.
(In the future, we implement the exact UTM transformation as a fallback —
contributions welcome!)

There is also `UTMfromUTMZ` and `UTMZfromUTM` transformations that are helpful
for converting between these two formats and putting data into the same `UTM`
zone.

#### To and from local `ENU` frames

The `ECEFfromENU` and `ENUfromECEF` transformations define the transformation
around a specific origin. Both the origin coordinates as an `ECEF` as well as
its corresponding latitude and longitude are stored in the transformation for
maximal efficiency when performing multiple `transform`s. The transformation can
be inverted with `inv` to perform the reverse transformation with respect to the
same origin.

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
* `ECEFfromUTMZ(datum) = ECEFfromLLA(datum) ∘ LLAfromUTMZ(datum)`
* `ENUfromUTMZ(origin, datum)  = ENUfromLLA(origin, datum) ∘ LLAfromUTMZ(datum`
* `UTMZfromENU(origin, datum)  = UTMZfromLLA(datum) ∘ LLAfromENU(origin, datum)`
* `UTMfromENU(origin, zone, hemisphere, datum) = UTMfromLLA(zone, hemisphere, datum) ∘ LLAfromENU(origin, datum)`
* `ENUfromUTM(origin, zone, hemisphere, datum) = ENUfromLLA(origin, datum) ∘ LLAfromUTM(zone, hemisphere, datum)`

Constructor-based transforms for these are also provided, such as `UTMZ(ecef, datum)`
which converts to `LLA` as an intermediary, as above. When converting multiple
points to or from the *same* ENU reference frame, it is recommended to use the
transformation-based approach for efficiency. However, the other
constructor-based conversions should be similar in speed to their transformation
counterparts.

### Distance

Currently, the only defined distance measure is the Cartesian distance,
`distance(x, y, [datum = wgs84])`, which works for all combinations of types for
`x` and `y` - except that the UTM zone and hemisphere must also be provided
for `UTM` types, as in `distance(utm1, utm2, zone, hemisphere, [datum = wgs84])`
(the Cartesian distance for `UTM` types is not approximated, but achieved via
conversion to `ECEF`).

This is the only function currently in
*Geodesy* which takes a default datum, and *should* be relatively accurate for
close points where Cartesian distances may be most important. Future work
may focus on geodesics and related calculations (contributions welcome!).
