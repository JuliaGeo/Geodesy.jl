using CoordinateTransformations

@testset "Datum transforms" begin

@testset "GDA94 from ITRF" begin
    # Tests for ITRF -> GDA94 conversion.  Unfortunately there doesn't seem to be a
    # public testset for this.

    @testset "Basic example" begin
        # Basic worked example from Dawson and Woods
        x_ITRF  = ECEF(-4052052.3678, 4212836.0411, -2545105.1089)
        x_GDA94 = ECEF(-4052051.7615, 4212836.1945, -2545106.0145)
        @test norm(x_GDA94 - Geodesy.GDA94_from_ITRF_Dawson2010(2005, Date(2010,6,16))(x_ITRF)) < 2e-4

        # Test inverse
        epoch = Date(2010,6,16)
        T = datum_shift_ECEF(ITRF{2008}(epoch), GDA94) ∘ datum_shift_ECEF(GDA94, ITRF{2008}(epoch))
        @test transformation_matrix(T) ≈ eye(3)
        @test isapprox(translation_vector(T), zeros(3), atol=10*eps())
    end

    @testset "ITRS realizations" begin
        x_GDA94 = ECEF(-4052051.7615, 4212836.1945, -2545106.0145)
        epoch = Date(2010,6,16)

        # Test dispatch
        @test datum_shift_ECEF(GDA94(), ITRF{2008}(epoch))(x_GDA94) == Geodesy.GDA94_from_ITRF_Dawson2010(2008, epoch)(x_GDA94)
        @test datum_shift_ECEF(GDA94(), ITRF{1996}(epoch))(x_GDA94) == Geodesy.GDA94_from_ITRF_Dawson2010(1996, epoch)(x_GDA94)

        # Nonexistent ITRF realization
        @test_throws Exception datum_shift_ECEF(ITRF{10000}(epoch), GDA94())
    end

    @testset "GDA94 from ITRF" begin
        # The following test data was obtained by submitting a GPS processing
        # request to the Geosciences Australia AUSPOS service, extracting the
        # cartesian coordinates of the base stations used in the solution from
        # the pdf report.

        # Cartesian, GDA94
        stations_GDA94 = [
            # Station X (m)         Y (m)        Z (m)
            #=BDST=# ECEF(-5021920.612, 2559339.872, -2975290.670),
            #=BNDY=# ECEF(-5125976.801, 2688801.593, -2669891.529),
            #=CBLT=# ECEF(-5061144.438, 2584178.830, -2886586.871),
            #=CNDO=# ECEF(-4494145.811, 2901721.598, -3462001.509),
            #=CWRA=# ECEF(-4532143.250, 2755335.963, -3531098.620),
            #=DALB=# ECEF(-4979266.855, 2730160.233, -2895220.922),
            #=GATT=# ECEF(-5012218.430, 2628002.744, -2931853.186),
            #=GONG=# ECEF(-4601821.645, 2561462.987, -3585678.678),
            #=IPS2=# ECEF(-5028440.649, 2588779.894, -2938802.549),
            #=MSVL=# ECEF(-4571853.513, 2599983.302, -3597310.167),
            #=PTKL=# ECEF(-4599805.569, 2558779.075, -3590076.333),
            #=ROBI=# ECEF(-5034843.827, 2523322.876, -2984064.629),
            #=RSBY=# ECEF(-5121088.692, 2863243.621, -2493144.670),
            #=TOOW=# ECEF(-4994481.724, 2663618.280, -2931172.448),
            #=WARW=# ECEF(-4967998.840, 2638152.833, -2997613.670),
        ]

        # Cartesian, ITRF2008
        ITRF_instance = ITRF{2008}(Date(2012,05,08))
        stations_ITRF2008 = [
            # Station X (m)         Y (m)       Z (m)
            #=BDST=# ECEF(-5021921.133, 2559339.665, -2975289.754),
            #=BNDY=# ECEF(-5125977.300, 2688801.332, -2669890.587),
            #=CBLT=# ECEF(-5061144.951, 2584178.607, -2886585.948),
            #=CNDO=# ECEF(-4494146.434, 2901721.514, -3462000.602),
            #=CWRA=# ECEF(-4532143.863, 2755335.886, -3531097.728),
            #=DALB=# ECEF(-4979267.387, 2730160.019, -2895219.988),
            #=GATT=# ECEF(-5012218.954, 2628002.532, -2931852.262),
            #=GONG=# ECEF(-4601822.242, 2561462.912, -3585677.803),
            #=IPS2=# ECEF(-5028441.169, 2588779.681, -2938801.629),
            #=MSVL=# ECEF(-4571854.116, 2599983.231, -3597309.290),
            #=PTKL=# ECEF(-4599806.166, 2558779.000, -3590075.459),
            #=ROBI=# ECEF(-5034844.345, 2523322.669, -2984063.717),
            #=RSBY=# ECEF(-5121089.191, 2863243.334, -2493143.707),
            #=TOOW=# ECEF(-4994482.252, 2663618.069, -2931171.521),
            #=WARW=# ECEF(-4967999.373, 2638152.635, -2997612.749),
        ]

        to_GDA = datum_shift_ECEF(GDA94(), ITRF_instance)
        # Need a nonzero tolerance which is expected, but it should be half
        # this size if (1) GA correctly rounded the coordinates (2) we used
        # exactly the same time computation and (3) list of transformation
        # coefficients.
        tolerance = 1.1e-3

        for i = 1:size(stations_ITRF2008,1)
            x_ITRF = stations_ITRF2008[i]
            x_GDA  = stations_GDA94[i]
            @test norm(x_GDA - to_GDA(x_ITRF)) < tolerance
        end
    end

end


end # @testset
