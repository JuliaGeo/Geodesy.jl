@testset "Co-ordinate system conversion" begin
    @testset "Fixed conversions" begin
        ###################################
        ### Testing fixed relationships ###
        ###################################

        lla = LLA(42.3673, -71.0960, 0.)
        lla_ref = LLA(42.36299, -71.09183, 0.)

        # LLA <-> ECEF
        ecef = ECEF(1529073.1560519305, -4465040.019013103, 4275835.339260307)
        @test ECEF(lla, wgs84) ≈ ecef
        @test LLA(ecef, wgs84) ≈ lla

        # LLA <-> ENU
        enu = ENU(-343.493749083977, 478.764855466788, -0.027242885224325164)
        @test ENU(lla, lla_ref, wgs84) ≈ enu
        @test LLA(enu, lla_ref, wgs84) ≈ lla

        # ECEF <-> ENU
        @test ENU(ecef, lla_ref, wgs84) ≈ enu
        @test ECEF(enu, lla_ref, wgs84) ≈ ecef

        # LLA <-> UTM
        (z, h) = (19, true)
        utm = UTM(327412.48528248386, 4.692686244318043e6, 0.0)
        @test UTM(lla, z, h, wgs84) ≈ utm
        @test LLA(utm, z, h, wgs84) ≈ lla

        # ECEF <-> UTM
        @test UTM(ecef, z, h, wgs84) ≈ utm
        @test ECEF(utm, z, h, wgs84) ≈ ecef

        # ENU <-> UTM
        @test UTM(enu, z, h, lla_ref, wgs84) ≈ utm
        @test ENU(utm, z, h, lla_ref, wgs84) ≈ enu

        # LLA <-> UTMZ
        utmz = UTMZ(327412.48528248386, 4.692686244318043e6, 0.0, z, h)
        @test UTMZ(lla, wgs84) ≈ utmz
        @test LLA(utmz, wgs84) ≈ lla

        # ECEF <-> UTMZ
        @test UTMZ(ecef, wgs84) ≈ utmz
        @test ECEF(utmz, wgs84) ≈ ecef

        # ENU <-> UTMZ
        @test UTMZ(enu, lla_ref, wgs84) ≈ utmz
        @test ENU(utmz, lla_ref, wgs84) ≈ enu

        # UTM <-> UTMZ
        @test UTMZ(utm, z, h, wgs84) ≈ utmz
        @test UTM(utmz, z, h, wgs84) ≈ utm

        # Identity
        @test ECEF(ecef, wgs84) == ecef
        @test LLA(lla, wgs84) == lla
        @test ENU(enu, wgs84) == enu
        @test UTM(utm, wgs84) == utm
        @test UTMZ(utmz, wgs84) == utmz
    end

    @testset "Random distance tests" begin
        #############################
        ### Testing random errors ###
        #############################
        """
            randLLA(h_min = -5e6, h_max = 5e6)

        Random LLA coordinate with evenly-distributed angles (probabilty density
        uniform over the surface-area of a sphere)
        """
        function randLLA(h_min = -5e6, h_max = 5e6)
            lon = 180*(2*rand()-1)
            lat = asind(2*rand()-1)
            h = h_min + (h_max-h_min)*rand()
            return LLA(lat,lon,h)
        end

        number_of_utm_distance_tests = 0
        for _ = 1:5000
            lla1 = randLLA()
            lla2 = randLLA()

            d = distance(lla1, lla2)

            ecef1 = ECEF(lla1, wgs84)
            ecef2 = ECEF(lla2, wgs84)

            @test distance(ecef1, ecef2) ≈ sqrt((ecef1.x - ecef2.x)^2 + (ecef1.y - ecef2.y)^2 + (ecef1.z - ecef2.z)^2)
            @test distance(ecef1, ecef2) ≈ d

            utmz1 = UTMZ(lla1, wgs84)
            utmz2 = UTMZ(lla2, wgs84)

            @test distance(utmz1, utmz2) ≈ d

            if utmz1.zone == utmz2.zone && utmz1.isnorth == utmz2.isnorth
                number_of_utm_distance_tests += 1

                utm1 = UTM(utmz1)
                utm2 = UTM(utmz2)

                zone = utmz1.zone
                isnorth = utmz1.isnorth

                @test distance(utm1,  utm2,  zone, isnorth) ≈ d
                @test distance(utm1,  ecef2, zone, isnorth) ≈ d
                @test distance(ecef1, utm2,  zone, isnorth) ≈ d
            end

            enu000 = ENU(0.0, 0.0, 0.0)
            enu = ENU(ecef1, ecef2, wgs84)

            @test distance(enu, enu000) ≈ d
        end
        @test number_of_utm_distance_tests > 0
    end
end # @testset
