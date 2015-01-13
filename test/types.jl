using Geodesy
using Base.Test

x, y = (rand(2) - .5) * 1000

@test ENU(x, y) == ENU(x, y, 0.0)

@test LLA(x, y) == LLA(x, y, 0.0)

ECEF(x, y, 0.0)

@test_throws TypeError Bounds{ECEF}(y, y + 1, x, x + 1)

Bounds{ENU}(y, y + 1, x, x + 1)

x, y = (rand() - .5) * 360, (rand() - .5) * 180

Bounds{LLA}(y, y + 1, x, x + 1)

Bounds(y, y, x, x)
@test_throws ArgumentError Bounds(y, y - 1, x, x)
@test_throws ArgumentError Bounds(y, y, x, x - 1)

@test_throws ArgumentError Bounds(-90.1, 90.1, 0, 1)
@test_throws ArgumentError Bounds(0, 1, -180.1, 180.1)
