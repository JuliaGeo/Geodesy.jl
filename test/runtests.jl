using Geodesy
using Test

################################################
### Helpers for testing approximate equality ###
################################################

# Interesting that this isn't in Base...
Base.isapprox(a::T, b::T; kwargs...) where {T<:Union{Tuple,NamedTuple}} = all(ntuple(i->isapprox(a[i],b[i]), length(a)); kwargs...)

@testset "Geodesy" begin
    include("points.jl")
    include("datums.jl")
    include("transverse_mercator.jl")
    include("polar_stereographic.jl")
    include("utm.jl")
    include("transformations.jl")
    include("conversion.jl")
    include("datum_transformations.jl")
    include("geodesics.jl")
end # @testset "Geodesy"
