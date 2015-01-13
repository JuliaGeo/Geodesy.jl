using Geodesy
using Base.Test

@test distance(ENU(1, 1, 1), ENU(2, 2, 2)) == sqrt(3)

@test distance(ECEF(1, 1, 1), ECEF(3, 3, 3)) == sqrt(12)

randLLA() = LLA((rand() - .5) * 178,
                (rand() - .5) * 358,
                (rand() - .5) * 18000)

for _ = 1:1_000
    lla = randLLA()
    lla2 = randLLA()

    enu = ENU(lla, lla)
    enu2 = ENU(lla2, lla)

    ecef = ECEF(lla)
    ecef2 = ECEF(lla2)

    @test_throws MethodError distance(lla, lla2)
    @test_approx_eq distance(enu, enu2) distance(ecef, ecef2)
end
