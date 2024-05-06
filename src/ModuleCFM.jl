module cfm

using ..Plume: extractsingleplume
using ..Datacontainer: CrossSectionalFluxData, ResolutionData, PlumeParams, RectGrid2d
using ..SatelliteRetrieval: retrivalbasedonpixelsigma!
using StatsBase: mean
using ..Constants


struct line
    dist::Float64
    vel::Vector{Float64}
    ds::Float64
    coordindex::Vector{CartesianIndex{2}}
end

mutable struct cfmdatavars
    tlines::Vector{line}
    no_tlines::Int64
    lineflux::Vector{Float64}
end

mutable struct simulationvars
    plume::PlumeParams
    background::Float64
    ppm_to_gms::Float64
end


function transactionlinesinxdirection(plume::BitMatrix, data::ResolutionData, interpolate)
    # get all cartesian coordinate indices of the plume 
    ix = findall(plume)
    # get only x-direction indices and sort it  
    xids = sort(unique(getindex.(ix, 1))) # tlines ids
    # number of x-indices is also number of transaction lines
    no_tlines = length(xids)
    # create an empty array of tranaction lines
    tlines = Vector{line}(undef, no_tlines)
    for i = 1:no_tlines
        # get all xy-indices for a point in x
        c_idx = ix[(getindex.(ix, 1).==xids[i])]
        # distance from the source
        dist = data.grid.x[xids[i]] - data.source[1]
        # compute the velocity at the x,y location
        # loop over all carteitian indices
        nt = length(c_idx)  # number of pts on a tline
        vel = zeros(nt)
        for k = 1:nt
            x1 = data.grid.x[c_idx[k][1]]
            y1 = data.grid.y[c_idx[k][2]]
            vel[k] = interpolate(x1, y1)
        end
        # create a transaction line
        tlines[i] = line(dist, vel, data.grid.ds, c_idx)
    end
    return cfmdatavars(tlines, no_tlines, zeros(no_tlines))
end


function compute_emission(noisy::Matrix{Float64}, cfmdata::cfmdatavars, ppm_to_gms::Float64)
    # compute line emissions and then a mean
    for i = 1:cfmdata.no_tlines
        cfmdata.lineflux[i] = (sum(noisy[cfmdata.tlines[i].coordindex] .* cfmdata.tlines[i].vel) .*
                               (cfmdata.tlines[i].ds * ppm_to_gms))
    end
    return mean(cfmdata.lineflux)
end


function get_gasconstants(gas::String)
    if gas == "co2"
        return Constants.co2.background, Constants.co2.ppm_to_gms
    elseif gas == "ch4"
        return Constants.ch4.background, Constants.ch4.ppm_to_gms
    end
end


function emission_cfm(samples::Int64, level2precision::Float64, data::ResolutionData, params::PlumeParams, gas::String, interpolate)
    # get background and conversion factor
    background, ppm_to_gms = get_gasconstants(gas)
    # level2precision::Float64, params, plumemask::BitMatrix)
    cfmemission = zeros(samples)
    # Compute CFM emissions
    sh = size(data.conc)
    # define pixel noise
    perpixel_precision = (data.conc .+ background) * level2precision

    # create transaction lines based on plumemask
    # CFM data vars
    plume = extractsingleplume(data.conc, params, data.originindex, data.plume)
    cfmvars = transactionlinesinxdirection(plume, data, interpolate)

    # compute noise and the emission
    noisy = zeros(Float64, sh)
    for i = 1:samples
        # get noisy data
        retrivalbasedonpixelsigma!(i, perpixel_precision, data.conc, noisy)
        # get cfm emission
        cfmemission[i] = compute_emission(noisy, cfmvars, ppm_to_gms)
    end
    @. cfmemission *= Constants.gmspersec_2_MTperyear
    return cfmemission, plume
end


function find_emission_7pixels(conc, emis, thresh)
    val = partialsort(vec(conc), 7; rev=true)
    return emis*(thresh * 1.001)/val
end

function get_data_at_emission(data, emis)
    conc = data.conc.*(emis/data.emission)
    return  ResolutionData(data.grid, data.source, emis, conc, data.originindex, 
    data.distancefromsource, data.plume)
end

function cfm_differentemissions(samples::Int64, level2precision::Float64, data::ResolutionData, params::PlumeParams, gas::String, interpolate)
    st_emis = find_emission_7pixels(data.conc, data.emission, params.threshold)
    allemissions = st_emis:0.05:20.0
    cfmemis = zeros(allemissions.len,samples)
    for i = 1:allemissions.len
        dataemis = get_data_at_emission(data, allemissions[i])
        cfmemis[i,:], plume = emission_cfm(samples, level2precision, dataemis, params, gas, interpolate)
    end
    return collect(allemissions), cfmemis
end

function emission_wholedomain(data)
    sh = size(data.conc)
    i1 = findmin(abs.(data.grid.x .- data.source[1]))[2]
    em = sum(data.conc[i1:end,:] .* data.vel[i1:end,:])/(sh[1]-i1)  # emissions in ppm
    return em*Constants.co2.ppm_to_gms*Constants.gmspersec_2_MTperyear*data.grid.dx
end


end # end module
