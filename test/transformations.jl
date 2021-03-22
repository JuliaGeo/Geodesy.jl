@testset "Coordinate system transformations" begin
    # Global coordinates
    lla = LLA(lat=-24.007712762068806, lon=96.44002819812094, alt=1374.7804632852078)
    lla2 = LLA(lat=-24.1, lon=96.5, alt=1659.18948659927)

    ecef_lla = ECEFfromLLA(wgs84)
    lla_ecef = LLAfromECEF(wgs84)

    ecef = ecef_lla(lla)
    ecef2 = ecef_lla(lla2)

    @test ecef ≈ ECEF(-654007.2768949205,5.794061748224952e6,-2.5796231395720555e6)
    @test lla ≈ lla_ecef(ecef)

    lla_utm = LLAfromUTM(utm_zone(lla)..., wgs84)
    utm_lla = UTMfromLLA(utm_zone(lla)..., wgs84)

    utm = utm_lla(lla)
    utm2 = utm_lla(lla2)

    @test utm ≈ UTM(239579.67583179142,7.342551466042985e6,1374.7804632852078)
    @test lla ≈ lla_utm(utm)
    # Test also the poles
    @test UTMfromLLA(0,true,wgs84)(LLA(89.0, 89.0, 89.0)) ≈ UTM(2.111009610242531e6, 1.9980623200455606e6, 89.0)
    @test LLAfromUTM(0,true,wgs84)(UTM(2.111009610242531e6, 1.9980623200455606e6, 89.0)) ≈ LLA(89.0, 89.0, 89.0)

    lla_utmz = LLAfromUTMZ(wgs84)
    utmz_lla = UTMZfromLLA(wgs84)

    utmz = utmz_lla(lla)
    utmz2 = utmz_lla(lla2)
    # Test also the poles go to UPS
    @test utmz_lla(LLA(89.0, 89.0, 89.0)) ≈ UTMZ(2.111009610242531e6, 1.9980623200455606e6, 89.0, 0, true)
    @test lla_utmz(UTMZ(2.111009610242531e6, 1.9980623200455606e6, 89.0, 0, true)) ≈ LLA(89.0, 89.0, 89.0)

    @test utmz ≈ UTMZ(239579.67583179142,7.342551466042985e6,1374.7804632852078, UInt8(47), false)
    @test lla ≈ lla_utmz(utmz)

    utmz_utm = UTMZfromUTM(47, false, wgs84)
    utm_utmz = UTMfromUTMZ(47, false, wgs84)

    @test utmz == utmz_utm(utm)
    @test utm == utm_utmz(utmz)
    # Sometimes this has to do something non-trivial - make sure it does
    utm_utmz2 = UTMfromUTMZ(46, false, wgs84)
    @test utm_utmz2(utmz) ≈ UTM(850009.7418418773, 7.340641097689279e6, 1374.7804632852078)

    # ENU coordinates
    enu_ecef = ENUfromECEF(lla2, wgs84)
    ecef_enu = ECEFfromENU(lla2, wgs84)

    enu = enu_ecef(ecef)

    @test enu_ecef ≈ ENUfromECEF(ecef2, wgs84) # (Can supply origin an any format)
    @test enu_ecef ≈ ENUfromECEF(utm2, 47, false, wgs84)
    @test enu_ecef ≈ ENUfromECEF(utmz2, wgs84)

    @test enu ≈ ENU(-6103.186938282723, 10222.54767562731, -295.55857190545794)
    @test ecef_enu(ENU(-6103.186938282723, 10222.54767562731, -295.55857190545794)) ≈ ecef

    # https://github.com/JuliaGeo/Geodesy.jl/issues/49
    @test LLA(ECEF(1, 0, 0), wgs84) ≈ LLA(89.99866260445, 0.00000000000, -6356752.314234)

    # Inverses
    @test inv(ecef_lla) == lla_ecef
    @test inv(lla_ecef) == ecef_lla
    @test inv(utm_lla) == lla_utm
    @test inv(lla_utm) == utm_lla
    @test inv(utmz_lla) == lla_utmz
    @test inv(lla_utmz) == utmz_lla
    @test inv(utmz_utm) == utm_utmz
    @test inv(utm_utmz) == utmz_utm
    @test inv(enu_ecef) == ecef_enu
    @test inv(ecef_enu) == enu_ecef

    # Web Mercator
    points = [LLA(0, 0, 0),
              LLA(49, 2, 0),
			  LLA(-19, -128, 4),
			  LLA(36, -77, 6)]
    # Reference generated using Proj4
    # [Proj4.transform(Proj4.Projection("+proj=longlat +datum=WGS84"),
    #                  Proj4.Projection("+proj=webmerc +datum=WGS84"),
    #                  [point.lon, point.lat, point.alt]) for point in points]
    points_wm = [[0.0, 0.0, 0.0],
                 [222638.98158654713, 6.274861394006577e6, 0.0],
                 [-1.4248894821539015e7, -2.1549359150858945e6, 4.0],
                 [-8.571600791082064e6, 4.300621372044271e6, 6.0]]
    wm_lla = WebMercatorfromLLA()
    lla_wm = LLAfromWebMercator()
    @testset "WebMercator" for (lla,wm) in zip(points, points_wm)
        @test wm_lla(lla) ≈ wm
        @test lla_wm(wm) ≈ lla
        ll = LatLon(lla.lat, lla.lon)
        wm_2d = wm[1:2]
        @test wm_lla(ll) ≈ wm_2d
        @test lla_wm(wm_2d) ≈ ll
	end

    # Composed transformation functions
    @test LLAfromENU(ecef, wgs84) == LLAfromECEF(wgs84) ∘ ECEFfromENU(ecef, wgs84)
    @test ENUfromLLA(ecef, wgs84) == ENUfromECEF(ecef, wgs84) ∘ ECEFfromLLA(wgs84)

    @test UTMfromECEF(47, false, wgs84) == UTMfromLLA(47, false, wgs84) ∘ LLAfromECEF(wgs84)
    @test ECEFfromUTM(47, false, wgs84) == ECEFfromLLA(wgs84) ∘ LLAfromUTM(47, false, wgs84)

    @test UTMZfromECEF(wgs84) == UTMZfromLLA(wgs84) ∘ LLAfromECEF(wgs84)
    @test ECEFfromUTMZ(wgs84) == ECEFfromLLA(wgs84) ∘ LLAfromUTMZ(wgs84)

    @test ENUfromUTMZ(ecef, wgs84) == ENUfromLLA(ecef, wgs84) ∘ LLAfromUTMZ(wgs84)
    @test UTMZfromENU(ecef, wgs84) == UTMZfromLLA(wgs84) ∘ LLAfromENU(ecef, wgs84)

    @test ENUfromUTM(ecef, 47, false, wgs84) == ENUfromLLA(ecef, wgs84) ∘ LLAfromUTM(47, false, wgs84)
    @test UTMfromENU(ecef, 47, false, wgs84) == UTMfromLLA(47, false, wgs84) ∘ LLAfromENU(ecef, wgs84)
end
