using Geodesy
using Base.Test

################################################
### Helpers for testing approximate equality ###
################################################

macro xy_approx_eq(a, b)
    quote
        @test_approx_eq ($(esc(a))).(1) ($(esc(b))).(1)
        @test_approx_eq ($(esc(a))).(2) ($(esc(b))).(2)
    end
end
macro xy_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps ($(esc(a))).(1) ($(esc(b))).(1) $(esc(eps))
        @test_approx_eq_eps ($(esc(a))).(2) ($(esc(b))).(2) $(esc(eps))
    end
end

macro xyz_approx_eq(a, b)
    quote
        @test_approx_eq ($(esc(a))).(1) ($(esc(b))).(1)
        @test_approx_eq ($(esc(a))).(2) ($(esc(b))).(2)
        @test_approx_eq ($(esc(a))).(3) ($(esc(b))).(3)
    end
end
macro xyz_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps ($(esc(a))).(1) ($(esc(b))).(1) $(esc(eps))
        @test_approx_eq_eps ($(esc(a))).(2) ($(esc(b))).(2) $(esc(eps))
        @test_approx_eq_eps ($(esc(a))).(3) ($(esc(b))).(3) $(esc(eps))
    end
end

###################################
### Testing fixed relationships ###
###################################

lla = LLA(42.3673, -71.0960, 0)
lla_ref = LLA(42.36299, -71.09183, 0)

# LLA -> ECEF
ecef_ref = ECEF(1529073.1560519305, -4465040.019013103, 4275835.339260309)
@xyz_approx_eq ECEF(lla) ecef_ref

#LLA -> ENU
enu_ref = ENU(-343.493749083977, 478.764855466788, -0.027242885224325164)
@xyz_approx_eq_eps ENU(lla, lla_ref) enu_ref 1e-8

#############################
### Testing random errors ###
#############################

randLLA() = (rand() - .5) * 178, (rand() - .5) * 360, (rand() - .5) * 18000

for _ = 1:50_000
    y, x, z = randLLA()
    min_x = x < -179 ? x + 359 : x - 1
    max_x = x >  179 ? x - 359 : x + 1
    lla = LLA(y, x, z)

    y, x, z = randLLA()
    min_x = x < -179 ? x + 359 : x - 1
    max_x = x >  179 ? x - 359 : x + 1
    lla2 = LLA(y, x, z)

    ecef = ECEF(lla)

    @xyz_approx_eq_eps LLA(ecef) lla 1e-6

    enu000 = ENU(0.0, 0.0, 0.0)

    @xyz_approx_eq ENU(ecef, lla) enu000

    enu2 = ENU(ecef, lla2)

    @xyz_approx_eq ENU(lla, lla2) enu2

end
