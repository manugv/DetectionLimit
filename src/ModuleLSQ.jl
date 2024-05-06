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


function getplume(resdata)
    notfound = true
    i = 7
    while notfound
        th = partialsort(vec(resdata.conc), i, rev=true)
        plumemask = resdata.conc .>= th
        segmentimg = zeros(Int64, size(resdata.conc))
        actualplume = resdata.conc .>= 1e-6
        segmentimage!(segmentimg, plumemask)
        goodpixels, _, plumelabel = plumepixelid(segmentimg, resdata.originindex, actualplume)
        if goodpixels >= 7
            plume = zeros(Bool, size(segmentimg))
            plume = segmentimg .== plumelabel
            return plume
            notfound = false
        end
        i = i + 1
    end
end


function flux_at_level4precision(data::ResolutionData, simparams::GeneralSimulationparameters, param_lsq::LSQParams, sigmanoise::Float64, emission)
    # compute plumemask for a given threshold. This remains constant for changing emission.
    # get top 7 pixels
    plumemask = getplume(data)

    # Compute certain emission and its LSQ
    lsq = zeros(emission.len)
    conc_at_emission = similar(data.conc)
    actplume = falses(size(data.conc))
    for i=1:emission.len
        # compute emission
        @. conc_at_emission = data.conc * emission[i]
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
    return emission, lsq
end


"""
    get_data_for_lsq(data::Datacont, emission::Float64, conc=false)

TBW
"""
function get_data_for_lsq(data::ResolutionData, emission::Float64, conc=false)
    conc_at_emission = zeros(size(data.gas.conc))
    @. conc_at_emission = data.gas.conc * emission
    plumemask = get_plumemask(conc_at_emission, data.threshold, data.thresholdflag)
    if data.no2.flag
        no_error = get_co2_noise_from_no2(conc_at_emission, data.gas.dist, plumemask, data.no2)
        return conc_at_emission[plumemask], no_error[plumemask]
    else
        return conc_at_emission[plumemask], 0
    end
end


function fluxestimate_at_precision(sigmanoise::Float64, init_emission::Float64, data::ResolutionData, lvl4precision::Float64)
    emission = init_emission
    # Variables to compute detection limit
    opt_cont = optimizecontainer(zeros(2), zeros(2), falses(2), 0.0)
    found_emission = true
    kk = 0
    new_std = 0
    # Loop till the detection limit is found
    while found_emission
        kk += 1
        if kk == 1000
            println("large number of iterations, optimize not working")
            println(opt_cont.emis, "   ", opt_cont.sd)
            return emission, new_std
        end
        # some times things are not successful as the threshold for plumemask 
        # is much higher than the avaliable enhancements
        # If threshold is too high then success_info is false
        ydata, co2_error_no2 = get_data_for_lsq(data, emission)
        gas_noise_vector = get_noise(data.background_ppm, sigmanoise, ydata)
        if data.no2.flag
            noise_vector = inverse_variance(gas_noise_vector, co2_error_no2)
        else
            noise_vector = gas_noise_vector .^ 2
        end
        new_std = compute_precision(noise_vector, ydata)
        # make the l4 precision relative
        new_std /= emission
        # println(new_std, " ", emission, ";")
        println("        STD of LSQ = ", new_std, "  for emission= ", emission, "  ", kk)
        info = get_optimization_flag(emission, new_std, opt_cont, lvl4precision)
        if info
            found_emission = false
        else
            emission, opt_cont = optmize_emission_mid(emission, new_std, opt_cont, lvl4precision)
        end
    end
    return emission, new_std
end


end # end module
