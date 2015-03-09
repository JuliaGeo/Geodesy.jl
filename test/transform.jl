using Geodesy
using Base.Test
using Compat

#############################################
### Decimal <=> Degrees, Minutes, Seconds ###
#############################################

for (decimal, d, m, s) in [(0.013, 0.0, 0.0, 46.8),
                           (-0.013, -0.0, 0.0, 46.8),
                           (-0.263, -0.0, 15.0, 46.8),
                           (-179.51, -179.0, 30.0, 36.0)]
    @test Geodesy.dms2decimal(d, m, s) === decimal
    d2, m2, s2 = Geodesy.decimal2dms(decimal)
    @test d2 === d
    @test m2 === m
    @test_approx_eq s2 s
end

################################################
### Helpers for testing approximate equality ###
################################################

# TODO: Move this to Compat.jl
if VERSION < v"0.4.0-dev+3616"
    fieldnames = names
end

macro type_approx_eq(a, b)
    quote
        @test fieldnames($(esc(a))) == fieldnames($(esc(b)))
        for n in fieldnames($(esc(a)))
            @test_approx_eq $(esc(a)).(n) $(esc(b)).(n)
        end
    end
end

Geodesy.getX(ecef::ECEF) = ecef.x
Geodesy.getY(ecef::ECEF) = ecef.y
Geodesy.getZ(ecef::ECEF) = ecef.z

macro xyz_approx_eq(a, b)
    quote
        @test_approx_eq getX($(esc(a))) getX($(esc(b)))
        @test_approx_eq getY($(esc(a))) getY($(esc(b)))
        @test_approx_eq getZ($(esc(a))) getZ($(esc(b)))
    end
end
macro xy_approx_eq(a, b)
    quote
        @test_approx_eq getX($(esc(a))) getX($(esc(b)))
        @test_approx_eq getY($(esc(a))) getY($(esc(b)))
    end
end

macro xyz_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getX($(esc(a))) getX($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getY($(esc(a))) getY($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getZ($(esc(a))) getZ($(esc(b))) $(esc(eps))
    end
end
macro xy_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getX($(esc(a))) getX($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getY($(esc(a))) getY($(esc(b))) $(esc(eps))
    end
end
macro z_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getZ($(esc(a))) getZ($(esc(b))) $(esc(eps))
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

# Bounds{LLA} -> Bounds{ENU}
bounds = Bounds(42.365, 42.3695, -71.1, -71.094)
bounds_enu_ref = Bounds{ENU}(-249.9308954374605, 249.9353534128848, -247.1268196136449, 247.1268196138187)
@type_approx_eq ENU(bounds) bounds_enu_ref

###################################
### Testing datum relationships ###
###################################

ecef = ECEF(5.953150599314804e6, 1.5951418955072558e6, 1.6403589592409942e6)
ecef2wgs = LLA{WGS84}(ecef)
ecef2nad = LLA{NAD27}(ecef)
ecef2osgb = LLA{OSGB36}(ecef)

@test LLA(ecef) == ecef2wgs
@test LL(ecef) == LL{WGS84}(ecef)

@test getX(ecef2wgs) == getX(ecef2nad)
@test abs(getY(ecef2wgs) - getY(ecef2nad)) > 1e-4
@test abs(getZ(ecef2wgs) - getZ(ecef2nad)) > 10

@test getX(ecef2wgs) == getX(ecef2osgb)
@test abs(getY(ecef2wgs) - getY(ecef2osgb)) > 1e-4
@test abs(getZ(ecef2wgs) - getZ(ecef2osgb)) > 10

#############################
### Testing random errors ###
#############################

randLLA() = (rand() - .5) * 178, (rand() - .5) * 360, (rand() - .5) * 18000

for _ = 1:50_000
    y, x, z = randLLA()
    min_x = x < -179 ? x + 359 : x - 1
    max_x = x >  179 ? x - 359 : x + 1
    lla = LLA(y, x, z)
    lla_bounds = Bounds{LLA}(y - 1, y + 1, min_x, max_x)
    ll = LL(y, x)
    ll_bounds = Bounds(y - 1, y + 1, min_x, max_x)

    y, x, z = randLLA()
    min_x = x < -179 ? x + 359 : x - 1
    max_x = x >  179 ? x - 359 : x + 1
    lla2 = LLA(y, x, z)
    lla2_bounds = Bounds{LLA}(y - 1, y + 1, min_x, max_x)
    ll2 = LL(y, x)
    ll2_bounds = Bounds(y - 1, y + 1, min_x, max_x)

    ecefa = ECEF(lla)
    ecef = ECEF(ll)
    @test_approx_eq_eps distance(ecef, ecefa) abs(getZ(lla)) 1e-8
    # TODO: could test proportionality

    @xyz_approx_eq_eps LLA(ecefa) lla 1e-6
    @xy_approx_eq_eps LL(ecefa) ll 1e-6

    @xy_approx_eq center(lla_bounds) lla
    @xy_approx_eq center(ll_bounds) ll

    enu000 = ENU(0.0, 0.0, 0.0)

    @xyz_approx_eq ENU(ecefa, lla) enu000

    @xy_approx_eq_eps ENU(ecefa, ll) enu000 1e-8
    @z_approx_eq_eps ENU(ecefa, ll) lla 1e-8

    ecefa2 = ECEF(lla2)
    ecef2 = ECEF(ll2)

    enu2 = ENU(ecefa, lla2)
    @xyz_approx_eq enu2 ENU(lla, lla2)

    @xy_approx_eq_eps enu2 ENU(ecefa, ll2) 1e-8
    zdiff = getZ(ENU(ecefa, ll2)) - getZ(enu2)
    @test_approx_eq_eps getZ(lla2) zdiff 1e-8

    # ECEF => ENU => ECEF w/ little change
    enu2v1 = ENU(ecef2, lla)
    @xyz_approx_eq_eps ECEF(enu2v1, lla) ecef2 1e-8

    # ENU => LL same as ENU => ECEF => LLA
    ecef2v1 = ECEF(enu2v1, lla)
    @xyz_approx_eq LLA(enu2v1, lla) LLA(ecef2v1)
    @xy_approx_eq LL(enu2v1, lla) LL(ecef2v1)

    dist_ecefa = distance(ecefa, ecefa2)
    dist_enua = distance(ENU(lla, lla), ENU(lla2, lla))
    @test_approx_eq dist_ecefa dist_enua

    dist_ecef = distance(ecef, ecef2)
    dist_enu = distance(ENU(ll, ll), ENU(ll2, ll))
    @test_approx_eq dist_ecef dist_enu
end
