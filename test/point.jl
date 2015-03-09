using Geodesy
using Base.Test

# Construction

x, y = (rand(2) - .5) * 10_000

@test ENU(x, y) == ENU(x, y, 0.0)

@test LLA(x, y) == LLA(x, y, 0.0)

ECEF(x, y, 0.0)

# get* methods

ll = LL(y, x)
lla = LLA(y, x, rand())
@test LLA(getY(lla), getX(lla), getZ(lla)) == lla
@test getY(ll) == y
@test getX(ll) == x
@test_throws MethodError getZ(ll)

enu = ENU(x, y, rand())
@test ENU(getX(enu), getY(enu), getZ(enu)) == enu
