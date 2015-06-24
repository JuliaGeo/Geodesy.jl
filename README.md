# Geodesy

[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Work with points defined in various coordinate systems.

The code has been split out from [OpenStreetMap.jl](https://github.com/tedsteiner/OpenStreetMap.jl), and functionality expanded.

Coordinate systems `LL`, `LLA`, `ECEF`, and `ENU` are supported. Transforms between between those types are supported.

`LL` and `LLA` are parameterized by a `Datum` -- constructors use `WGS84` as the default datum if none is specified, but after construction, the correct `Datum` is carried around for free for use in downstream transforms and distance calculations.

Currently, `Datum` is a bit of a placeholder, only including an `Ellipsoid`, so transforms to/from `LL` and `LLA` types are only correct for geocentric datums using the IERS reference meridian. This should be easy to extend, and a priority.

Defining and using a new `Datum` manually is trivial (a two-liner). A few common `Datum`s and `Ellipsoid`s are provided, but specifying any `Datum` or `Ellipsoid` by EPSG ID ought to happen (perhaps in an auxillary package), too.

The full list of types, constants, and methods provided is at the top of [src/Geodesy.jl](src/Geodesy.jl).
