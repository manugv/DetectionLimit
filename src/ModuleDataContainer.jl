module Datacontainer
using ..Constants

mutable struct NO2toCO2
    co2_to_nox::Float64
    absolute::Float64
    relative::Float64
    flag::Bool
end

struct GeneralSimulationparameters
    simulationname::String
    gas::Constants.gasconst
    outputfile::String
    resolution::Vector{Float64}
    level2precision::Vector{Float64}
    no2_to_co2::NO2toCO2
end

struct LSQParams
    level4precision::Float64
    cdfflag::Bool
    threshold::Float64
end

struct PlumeParams
    cdfflag::Bool
    threshold::Float64
end

struct Fileandkeys
    filename::String
    timekeys::Vector{String}
end

mutable struct DetectionLimitVars
    const samples::Int64
    const cdfflag::Bool
    threshold::Float64
    const numberofpixels::Int64
    const plumelength::Float64
end

struct Grid3d
    x::AbstractRange
    y::AbstractRange
    z::AbstractRange
    xn::AbstractRange
    yn::AbstractRange
    zn::AbstractRange
    dx::Float64
    dy::Float64
    dz::Float64
    dxi::Float64
    dyi::Float64
    dzi::Float64
    nx::Int32
    ny::Int32
    nz::Int32
    function Grid3d(nos, ds)
        st = @. ds / 2
        return new(LinRange(st[1], ds[1] * nos[1] - st[1], nos[1]),
            LinRange(st[2], ds[2] * nos[2] - st[2], nos[2]),
            LinRange(st[3], ds[3] * nos[3] - st[3], nos[3]),
            0:ds[1]:nos[1]*ds[1],
            0:ds[2]:nos[2]*ds[2],
            0:ds[3]:nos[3]*ds[3],
            ds[1], ds[2], ds[3],
            1.0 / ds[1], 1.0 / ds[2], 1.0 / ds[3],
            nos[1], nos[2], nos[3])
    end
end

struct Grid2d
    x::Vector{Float64}
    y::Vector{Float64}
    xn::Vector{Float64}
    yn::Vector{Float64}
    dx::Float64
    dy::Float64
    nx::Int32
    ny::Int32
end


struct RectGrid2d
    x::Vector{Float64}
    y::Vector{Float64}
    xn::Vector{Float64}
    yn::Vector{Float64}
    ds::Float64
    dsi::Float64
    nx::Int32
    ny::Int32
end


mutable struct Level2Retrieval
    conc::Matrix{Float64}
    plumemask::BitMatrix
    segmentedimage::Matrix{Int64}
    function Level2Retrieval(sh)
        return new(zeros(Float64, sh), falses(sh), zeros(Int64, sh))
    end
end


mutable struct ResolutionData
    grid::RectGrid2d
    originindex::Vector{Int64}
    distancefromsource::Matrix{Float64}
    conc::Matrix{Float64}
end


mutable struct EmissionData
    emission::Float64
    conc::Matrix{Float64}
    plume::BitMatrix
    perpixelsigma::Matrix{Float64}
    function EmissionData(emis, sh)
        return new(emis, zeros(Float64, sh), falses(sh), zeros(Float64, sh))
    end
end


mutable struct DetectionLimitData
    actual::ResolutionData
    init_emission::Float64
    emisdata::EmissionData
    groupname::String
    allemissions::Vector{Float64}
    probabilities::Vector{Float64}
    level2::Array{Level2Retrieval}
    function DetectionLimitData(resdata)
        sh = size(resdata.conc)
        lvl2r = [Level2Retrieval(sh) for _ in 1:Threads.nthreads()]
        emis = EmissionData(0.0, sh)
        return new(resdata, 0.0, emis, "", [], [], lvl2r)
    end
end


mutable struct ConcData
    const source::Vector{Float64}
    const emissionstrength::Float64
    const grid::Grid2d
    time::Int64
    conc::Matrix{Float64}
    convolved_conc::Matrix{Float64}    
end


# mutable struct EmissionData
#     emission::Float64
#     actualconc::Matrix{Float64}
#     actualplume::BitMatrix
#     perpixelsigma::Matrix{Float64}
#     function EmissionData(emis, sh)
#         return new(emis, zeros(Float64, sh), falses(sh), zeros(Float64, sh))
#     end
# end


# mutable struct DetectionlimitData

#     level2retrieval::Matrix{Float64}
#     groupname::String
#     allemissions::Vector{Float64}
#     probabilities::Vector{Float64}
# end

# mutable struct SatelliteData
#     init_emission::Float64
#     emission::Float64
#     groupname::String
#     allemissions::Vector{Float64}
#     probabilities::Vector{Float64}
#     actualconc::Matrix{Float64}
#     actualplume::BitMatrix
#     perpixelsigma::Matrix{Float64}
#     level2retrieval::Matrix{Float64}
#     plumemask::BitMatrix
#     segmentedimage::Matrix{Int64}
#     function SatelliteData(emis, sh)
#         return new(0.0, emis, "", [], [], zeros(Float64, sh), falses(sh), zeros(Float64, sh), zeros(Float64, sh), falses(sh), zeros(Int64, sh))
#     end
# end

mutable struct CrossSectionalFluxData
    noise::Matrix{Float64}
    plumemask::BitMatrix
    function CrossSectionalFluxData(sh, plume)
        return new(zeros(Float64, sh), plume)
    end
end

end # module end
