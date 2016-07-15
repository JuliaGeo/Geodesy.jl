using BaseTestNext


# Worked example from the Dawson paper
x_ITRF  = Vec(-4052052.3678, 4212836.0411, -2545105.1089)
x_GDA94 = Vec(-4052051.7615, 4212836.1945, -2545106.0145)
@test maximum(abs(x_GDA94 - GDA94_from_ITRF(2005, Date(2010,6,16))(x_ITRF))) < 2e-4


# The following test data was obtained by submitting a GPS processing request
# to the Geosciences Australia AUSPOS service, and extracting the geodetic
# coordinates of the base stations used in the solution from the pdf report.

"""
    dms2deg(degrees, minutes, seconds)

Convert (degrees, minutes, seconds) to degrees.  Degrees holds the sign for all
three components.
"""
function dms2deg(degrees, minutes, seconds)
    degrees + (copysign(minutes, degrees) + copysign(seconds, degrees)/60)/60
end

# Geodetic, GRS80 ellipsoid, GDA94
stations_GDA94 = [
    #=        lat deg min sec   lon deg min sec   height =#
    #=BDST=# -27  59  13.56957  152  59  42.27814  101.0990
    #=BNDY=# -24  54  29.62405  152  19  15.60714  80.1208
    #=CBLT=# -27   5   3.97236  152  57   5.45879  83.9391
    #=CNDO=# -33   5   6.58765  147   9   3.61478  229.7387
    #=CWRA=# -33  49  52.48768  148  42   8.31799  333.3676
    #=DALB=# -27  10  13.97535  151  15  49.65026  394.6879
    #=GATT=# -27  32  38.17782  152  19  51.99962  140.5849
    #=GONG=# -34  25  38.01233  150  53  55.82790  75.5994
    #=IPS2=# -27  36  53.76284  152  45  33.62941  88.6452
    #=MSVL=# -34  33   1.95879  150  22  24.21243  703.1676
    #=PTKL=# -34  28  31.99648  150  54  49.30154  34.5429
    #=ROBI=# -28   4  37.08923  153  22  52.50845  65.3016
    #=RSBY=# -23   9  39.58837  150  47  24.28228  58.2422
    #=TOOW=# -27  32   4.00312  151  55  42.43296  685.7681
    #=WARW=# -28  12  48.54186  152   1  49.40192  507.4264
]

# Geodetic, GRS80 ellipsoid, ITRF2008, 2012-05-08
stations_ITRF2008 = [
    #=        lat deg min sec   lon deg min sec   height =#
    #=BDST=# -27  59  13.53765  152  59  42.29353  100.996
    #=BNDY=# -24  54  29.59189  152  19  15.62365  80.015
    #=CBLT=# -27   5   3.94041  152  57   5.47447  83.835
    #=CNDO=# -33   5   6.55452  147   9   3.63053  229.644
    #=CWRA=# -33  49  52.45488  148  42   8.33296  333.273
    #=DALB=# -27  10  13.94298  151  15  49.66638  394.585
    #=GATT=# -27  32  38.14572  152  19  52.01533  140.482
    #=GONG=# -34  25  37.98002  150  53  55.84187  75.504
    #=IPS2=# -27  36  53.73085  152  45  33.64499  88.542
    #=MSVL=# -34  33   1.92636  150  22  24.22655  703.073
    #=PTKL=# -34  28  31.96418  150  54  49.31548  34.448
    #=ROBI=# -28   4  37.05740  153  22  52.52371  65.199
    #=RSBY=# -23   9  39.55581  150  47  24.29963  58.135
    #=TOOW=# -27  32   3.97093  151  55  42.44879  685.665
    #=WARW=# -28  12  48.50971  152   1  49.41751  507.324
]

@testset "GDA94 from ITRF" begin
    to_GDA_LLA = LLAfromECEF(grs80) ∘
                GDA94_from_ITRF(2008, Date(2012,05,08)) ∘
                ECEFfromLLA(grs80)

    for i = 1:size(stations_ITRF2008,1)
        x_ITRF = LLA(lat=dms2deg(stations_ITRF2008[i,1:3]...),
                    lon=dms2deg(stations_ITRF2008[i,4:6]...),
                    alt=stations_ITRF2008[7])
        x_GDA  = LLA(lat=dms2deg(stations_GDA94[i,1:3]...),
                    lon=dms2deg(stations_GDA94[i,4:6]...),
                    alt=stations_GDA94[7])
        x_GDA2 = to_GDA_LLA(x_ITRF)
        @test abs(x_GDA.lat - x_GDA2.lat) < 8.5e-8
        @test abs(x_GDA.lon - x_GDA2.lon) < 8.5e-8
        @test abs(x_GDA.alt - x_GDA2.alt) < 2e-2 # FIXME ??
    end
end
