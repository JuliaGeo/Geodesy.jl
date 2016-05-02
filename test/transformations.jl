@testset "Coordinate system transformations" begin

    lla = LLA(lat=-24.007712762068806, lon=96.44002819812094, alt=1374.7804632852078)
    lla2 = LLA(lat=-24.1, lon=96.5, alt=1659.18948659927)

    ecef_lla = ECEFfromLLA(wgs84)
    lla_ecef = LLAfromECEF(wgs84)

    ecef = transform(ecef_lla, lla)
    ecef2 = transform(ecef_lla, lla2)

    @test ecef ≈ ECEF(-654007.2768949205,5.794061748224952e6,-2.5796231395720555e6)
    @test lla ≈ transform(lla_ecef, ecef)

    



end
