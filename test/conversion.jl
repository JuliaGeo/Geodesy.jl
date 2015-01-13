using Geodesy
using Base.Test

################################################
### Helpers for testing approximate equality ###
################################################

macro type_approx_eq(a, b)
    quote
        @test names($(esc(a))) == names($(esc(b)))
        for n in names($(esc(a)))
            @test_approx_eq $(esc(a)).(n) $(esc(b)).(n)
        end
    end
end

macro xy_approx_eq(a, b)
    quote
        @test_approx_eq getX($(esc(a))) getX($(esc(b)))
        @test_approx_eq getY($(esc(a))) getY($(esc(b)))
    end
end
macro xy_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getX($(esc(a))) getX($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getY($(esc(a))) getY($(esc(b))) $(esc(eps))
    end
end

macro xyz_approx_eq(a, b)
    quote
        @test_approx_eq getX($(esc(a))) getX($(esc(b)))
        @test_approx_eq getY($(esc(a))) getY($(esc(b)))
        @test_approx_eq getZ($(esc(a))) getZ($(esc(b)))
    end
end
macro xyz_approx_eq_eps(a, b, eps, zeps)
    quote
        @test_approx_eq_eps getX($(esc(a))) getX($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getY($(esc(a))) getY($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getZ($(esc(a))) getZ($(esc(b))) $(esc(zeps))
    end
end

###################################
### Testing fixed relationships ###
###################################

lla = LLA(42.3673, -71.0960, 0)
lla_ref = LLA(42.36299, -71.09183, 0)

# LLA -> ECEF
ecef_ref = ECEF(1529073.1560519305, -4465040.019013103, 4275835.339260309)
@type_approx_eq ECEF(lla) ecef_ref

#LLA -> ENU
enu_ref = ENU(-343.49374908345493, 478.7648554687071, -0.027242884564486758)
@xyz_approx_eq_eps ENU(lla, lla_ref) enu_ref  1e-8 1e-10

bounds = Bounds(42.365, 42.3695, -71.1, -71.094)

# Bounds{LLA} -> Bounds{ENU}
bounds_enu_ref = Bounds{ENU}(-249.92653559014204, 249.9353534127322, -247.1091823451915, 247.1268196141778)
@type_approx_eq ENU(bounds) bounds_enu_ref

#############################
### Testing random errors ###
#############################

randLLA() = (rand() - .5) * 178, (rand() - .5) * 358, (rand() - .5) * 18000

for _ = 1:10_000
    y, x, z = randLLA()
    lla = LLA(y, x, z)
    lla_bounds = Bounds(y - 1, y + 1, x - 1, x + 1)

    y, x, z = randLLA()
    lla2 = LLA(y, x, z)
    lla2_bounds = Bounds(y - 1, y + 1, x - 1, x + 1)

    ecef = ECEF(lla)

    @xyz_approx_eq_eps LLA(ecef) lla 1e-7 1e-6

    @xy_approx_eq center(lla_bounds) lla

    enu000 = ENU(0.0, 0.0, 0.0)

    @xyz_approx_eq ENU(ecef, lla) enu000

    @xy_approx_eq_eps ENU(ecef, lla_bounds) enu000 1e-8
    # TODO: what's up with alt?

    @xy_approx_eq_eps ENU(lla, lla_bounds) enu000 1e-8

    enu2 = ENU(ecef, lla2)

    @xy_approx_eq_eps ENU(ecef, lla2_bounds) enu2 1e-8

    @xyz_approx_eq ENU(lla, lla2) enu2
    @xy_approx_eq_eps ENU(lla, lla2_bounds) enu2 1e-8

    @type_approx_eq ENU(lla_bounds) ENU(lla_bounds, center(lla_bounds))
end
