# Quantify the error of the 6th order approximation to transverse Mercator
# projections by the invertability of points

using Geodesy
using PyPlot

z = 31
h = true
t1 = UTMfromLLA(z,h,wgs84)
t2 = LLAfromUTM(z,h,wgs84)

lats = Float64[0:1:80;]
lons = Float64[0:1:80;]
errors = Float64[distance(LLA(lats[i],lons[j]+3.0,0.0), transform(t2, transform(t1, LLA(lats[i],lons[j]+3.0,0.0)))) for i=1:length(lats), j in 1:length(lons)]
errors = convert(Matrix{Float64}, errors) # Compehension is broken (code_warntype is correct for distance but not for the array...)

figure()
imshow(log10(errors))

error_lon = [maximum(errors[:,i]) for i in 1:length(lons)]

figure()
plot(lons, error_lon)
