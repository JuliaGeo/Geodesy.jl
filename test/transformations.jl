@testset "Coordinate system transformations" begin

    lla = LLA(lat=-24.007712762068806, lon=96.44002819812094, alt=1374.7804632852078)
    lla2 = LLA(lat=-24.1, lon=96.5, alt=1659.18948659927)

    ecef_lla = ECEFfromLLA(wgs84)
    lla_ecef = LLAfromECEF(wgs84)

    ecef = transform(ecef_lla, lla)
    ecef2 = transform(ecef_lla, lla2)

    @test ecef ≈ ECEF(-654007.2768949205,5.794061748224952e6,-2.5796231395720555e6)
    @test lla ≈ transform(lla_ecef, ecef)

    lla_utm = LLAfromUTM(utm_zone(lla)..., wgs84)
    utm_lla = UTMfromLLA(utm_zone(lla)..., wgs84)

    utm = transform(utm_lla, lla)

    @test utm ≈ UTM(239579.67583179142,7.342551466042985e6,1374.7804632852078)
    @test lla ≈ transform(lla_utm, utm)

    lla_utmz = LLAfromUTMZ(wgs84)
    utmz_lla = UTMZfromLLA(wgs84)

    utmz = transform(utmz_lla, lla)

    @test utmz ≈ UTMZ(239579.67583179142,7.342551466042985e6,1374.7804632852078, UInt8(47), false)
    @test lla ≈ transform(lla_utmz, utmz)


end
