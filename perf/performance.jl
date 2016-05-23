# Measure the performance of the transverse Mercator projections

import Geodesy
import Proj4

function perf_geodesy(n)
    lla = Geodesy.LLA(1.0,1.0,0.0)
    trans1 = Geodesy.UTMfromLLA(Geodesy.utm_zone(lla)..., Geodesy.wgs84)
    trans2 = Geodesy.LLAfromUTM(Geodesy.utm_zone(lla)..., Geodesy.wgs84)
    for i = 1:n
        utm = Geodesy.transform(trans1, lla)
        lla2 = Geodesy.transform(trans2, utm)
        #if !(lla ≈ lla2)
    #        error("Not invertable")
#        end
    end
end


function perf_proj4(n)
    wgs84 = Proj4.Projection("+proj=longlat +datum=WGS84 +no_defs")
    utm56 = Proj4.Projection("+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs")
    lla = [150.0, -27.0, 0.0]
    for i = 1:n
        utm = Proj4.transform(wgs84, utm56, lla)
        lla2 = Proj4.transform(utm56, wgs84, utm)
#        if !(lla ≈ lla2)
#            error("Not invertable")
#        end
    end
end

println("Comparing LLA->UTM->LLA transformation speed for Geodesy and Proj4")
println("\nwarming up...")
perf_geodesy(1)
@time perf_geodesy(1)
perf_proj4(1)
@time perf_proj4(1)

n = 1_000_000
println("\nPerforming (2×) $n transformations with Geodesy")
perf_geodesy(1)
@time perf_geodesy(n)
println("\nPerforming (2×) $n transformations with Proj4")
perf_proj4(1)
@time perf_proj4(n)

# 5/5/2016 1:40pm (~14 × slower) (probably more like 25 × slower excluding ≈ test)
#
# Comparing LLA->UTM->LLA transformation speed for Geodesy and Proj4
#
# warming up...
#   0.000042 seconds (469 allocations: 19.448 KB)
#   0.000026 seconds (18 allocations: 736 bytes)
#
# Performing (2×) 1000000 transformations with Geodesy
#  16.281498 seconds (315.00 M allocations: 8.717 GB, 4.76% gc time)
#
# Performing (2×) 1000000 transformations with Proj4
#   1.126085 seconds (12.00 M allocations: 488.281 MB, 2.59% gc time)
#
# ------------------------------------------------------------------------------
# 12/5/2016 1:35pm (~2.4 × slower, ran excluding ≈ test)
#
# Comparing LLA->UTM->LLA transformation speed for Geodesy and Proj4
#
# warming up...
#   0.000067 seconds (856 allocations: 32.526 KB)
#   0.000022 seconds (9 allocations: 512 bytes)
#
# Performing (2×) 1000000 transformations with Geodesy
#   1.445283 seconds (712 allocations: 22.531 KB)
#
# Performing (2×) 1000000 transformations with Proj4
#   0.608355 seconds (3.00 M allocations: 274.658 MB, 1.71% gc time)
