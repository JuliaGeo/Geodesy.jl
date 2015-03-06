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

# Distance

@test distance(ENU(1, 1, 1), ENU(2, 2, 2)) == sqrt(3)
@test distance(ECEF(1, 1, 1), ECEF(3, 3, 3)) == sqrt(12)
@test_throws MethodError distance(lla, lla)
@test_throws MethodError distance(ll, ll)

randLLA() = LLA((rand() - .5) * 180,
                (rand() - .5) * 360,
                (rand() - .5) * 18_000)

for _ = 1:1_000
    lla1 = randLLA()
    lla2 = randLLA()

    enu1a = ENU(lla1, lla1)
    enu2a = ENU(lla2, lla1)

    ecef1a = ECEF(lla1)
    ecef2a = ECEF(lla2)

    @test_approx_eq distance(enu1a, enu2a) distance(ecef1a, ecef2a)
end
