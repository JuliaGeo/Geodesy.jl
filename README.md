# Geodesy

[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Work with points defined in various coordinate systems.

The code has been split out from [OpenStreetMap.jl](https://github.com/tedsteiner/OpenStreetMap.jl), and expanded slightly, but current functionality is basically limited to what is needed by OpenStreetMap.jl.

Coordinate systems `LLA`, `ECEF`, and `ENU` are supported. Most translations between between those types are supported, and use `WGS84` as the default datum. A couple other datums are provided, as well as an `Ellipsoid` type for defining datums. `distance` takes either two `ENU` or two `ECEF` points, and `Bounds{LLA}` and `Bounds{ENU}` are useful for operations such as `inBounds(point, bounds)`.
