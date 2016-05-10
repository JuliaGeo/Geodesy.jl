@testset "Co-ordinate system conversion" begin
    ###################################
    ### Testing fixed relationships ###
    ###################################

    lla = LLA(42.3673, -71.0960, 0.)
    lla_ref = LLA(42.36299, -71.09183, 0.)

    # LLA -> ECEF
    ecef_ref = ECEF(1529073.1560519305, -4465040.019013103, 4275835.339260309)
    @xyz_approx_eq ECEF(lla, wgs84) ecef_ref

    #LLA -> ENU
    enu_ref = ENU(-343.493749083977, 478.764855466788, -0.027242885224325164)
    @xyz_approx_eq_eps ENU(lla, lla_ref, wgs84) enu_ref 1e-8

    #############################
    ### Testing random errors ###
    #############################
    """
        randLLA(h_min = -5e6, h_max = 5e6)

    Random LLA coordinate with evenly-distributed angles (according to the
    Haar measure over the sphere)
    """
    function randLLA(h_min = -5e6, h_max = 5e6)
        lon = 180*(2*rand()-1)
        lat = asind(2*rand()-1)
        h = h_min + (h_max-h_min)*rand()
        return LLA(lat,lon,h)
    end

    for _ = 1:50_000
        lla = randLLA()
        lla2 = randLLA()

        ecef = ECEF(lla, wgs84)
        @xyz_approx_eq_eps LLA(ecef, wgs84) lla 1e-6

        enu000 = ENU(0.0, 0.0, 0.0)
        @xyz_approx_eq ENU(ecef, lla, wgs84) enu000

        enu2 = ENU(ecef, lla2, wgs84)
        @xyz_approx_eq ENU(lla, lla2, wgs84) enu2
    end
end # @testset
