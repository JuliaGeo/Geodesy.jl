using Geodesy
using Base.Test

# Construction

x, y = (rand(2) - .5) * 10_000

@test ENU(x, y) == ENU(x, y, 0.0)

@test LLA(x, y) == LLA(x, y, 0.0)

ECEF(x, y, 0.0)

# Distance

@test distance(ENU(1, 1, 1), ENU(2, 2, 2)) == sqrt(3)

@test distance(ECEF(1, 1, 1), ECEF(3, 3, 3)) == sqrt(12)

randLLA() = LLA((rand() - .5) * 180,
                (rand() - .5) * 360,
                (rand() - .5) * 18_000)

for _ = 1:1_000
    lla = randLLA()
    lla2 = randLLA()

    enu = ENU(lla, lla, wgs84)
    enu2 = ENU(lla2, lla, wgs84)

    ecef = ECEF(lla, wgs84)
    ecef2 = ECEF(lla2, wgs84)

    @test_throws MethodError distance(lla, lla2)
    @test_approx_eq distance(enu, enu2) distance(ecef, ecef2)
end
