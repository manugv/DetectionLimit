module LSQ
using ..Datacontainer: ResolutionData, NO2toCO2, LSQParams, GeneralSimulationparameters
using ..NO2Error
using ..Plume: cdfthresholdbasedplumemask, thresholdbasedplumemask, segmentimage!, plumepixelid
using LinearAlgebra: Diagonal


"""
    get_precision(back::Float64, sigmanoise::Float64, ymod::Vector{Float64})

Compute precision for Least square estimates.
"""
function compute_precision(variance_noise::Vector{Float64}, ymod::Vector{Float64})
    sy = Diagonal(variance_noise)
    syinv = inv(sy)
    K = reshape(ymod, length(ymod), 1)
    return sqrt(inv(K' * (syinv' * K))[1])
end


"""
    lsq_precision_emission(data::ResolutionData, simparams::GeneralSimulationparameters, lsq::LSQParams, sigmanoise::Float64)

TBW
"""
function lsq_level4precision(data::ResolutionData, simparams::GeneralSimulationparameters,
     lsq::LSQParams, sigmanoise::Float64)
    
     # compute plumemask for a given threshold
    if lsq.cdfflag
        plumemask = cdfthresholdbasedplumemask(data.conc, lsq.threshold)
    else
        plumemask = thresholdbasedplumemask(data.conc, lsq.threshold)
    end

    # compute error vector for given conc
    ydata = data.conc[plumemask]

    # Compute error
    sigma_gas = (ydata .+ simparams.gas.background) * sigmanoise
    # If no2 error is included then compute error
    if simparams.no2_to_co2.flag
        sigmaco2_no2 = NO2Error.get_co2_noise_from_no2(data.conc, data.distancefromsource, data.plume, simparams.no2_to_co2)
        variance_gas = NO2Error.inversevarianceweighting(sigma_gas, sigmaco2_no2[plumemask])
    else
        variance_gas = sigma_gas .^ 2
    end
    # compute the precision
    return compute_precision(variance_gas, ydata)
end


function getplume_connect(resdata, no_pixels)
    # get top 7 connected  pixels to the source.
    notfound = true
    # initialize
    i = no_pixels
    while notfound
        th = partialsort(vec(resdata.conc), i, rev=true)
        plumemask = resdata.conc .>= th
        segmentimg = zeros(Int64, size(resdata.conc))
        actualplume = resdata.conc .>= 1e-6
        segmentimage!(segmentimg, plumemask)
        goodpixels, _, plumelabel = plumepixelid(segmentimg, resdata.originindex, actualplume)
        if goodpixels >= no_pixels
            plume = zeros(Bool, size(segmentimg))
            plume = segmentimg .== plumelabel
            return plume
            notfound = false
        end
        i = i + 1
    end
end


function getplume(resdata, no_pixels)
    # get top 7 pixels.
    th = partialsort(vec(resdata.conc), no_pixels, rev=true)
    plumemask = resdata.conc .>= th
    return plumemask
end


function flux_at_level4precision(data::ResolutionData, simparams::GeneralSimulationparameters, param_lsq::LSQParams, sigmanoise::Float64, method="Top")
    # compute plumemask for a given threshold. This remains constant for changing emission.
    if method == "Connected"
        # get top 7 connected pixels with the source
        plumemask = getplume_connect(data, param_lsq.no_pixels)
    else
        # get top 7 pixels
        plumemask = getplume(data, param_lsq.no_pixels)
    end

    # Compute certain emission and its LSQ
    lsq = zeros(param_lsq.emission.len)
    conc_at_emission = similar(data.conc)
    actplume = falses(size(data.conc))
    for i=1:param_lsq.emission.len
        # compute emission
        @. conc_at_emission = data.conc * param_lsq.emission[i]
        # compute error vector for given conc
        ydata = conc_at_emission[plumemask]

        # Compute error
        sigma_gas = (ydata .+ simparams.gas.background) * sigmanoise
        # If no2 error is included then compute error
        if simparams.no2_to_co2.flag
            # actual plume for no2
            @. actplume = false
            @. actplume = conc_at_emission > 1e-6
            sigmaco2_no2 = NO2Error.get_co2_noise_from_no2(conc_at_emission, data.distancefromsource, actplume, simparams.no2_to_co2)
            variance_gas = NO2Error.inversevarianceweighting(sigma_gas, sigmaco2_no2[plumemask])
        else
            variance_gas = sigma_gas .^ 2
        end
        # compute the precision
        lsq[i] = compute_precision(variance_gas, ydata)
    end
    return lsq
end

end # end module
