module SatelliteRetrieval

using Random: randn!, MersenneTwister
using ..Datacontainer: EmissionData, NO2toCO2, ResolutionData
using ..NO2Error


"""
    get_noise(noi_p::Float64, noisy::Matrix{Float64})

Get normally distributed noise based on standard deviation 'sigma'
"""
function computenoise!(sigma::Float64, noisy::Matrix{Float64})
    randn!(noisy)
    @. noisy = sigma * noisy
    nothing
end


"""
    enhancementdependentretrival!(sigma::Float64, data_ppm::Matrix{Float64}, background::Float64, noisy::Matrix{Float64})

Get synthetic satellite retrival where noise is dependent on enhancement. 
"""
function retrivalbasedonpixelsigma!(i::Int64, perpixelsigma::Matrix{Float64}, data_ppm::Matrix{Float64}, noisy::Matrix{Float64})
    randn!(MersenneTwister(i), noisy)
    @. noisy = data_ppm +  perpixelsigma * noisy
end

"""
    enhancementdependentretrival!(sigma::Float64, data_ppm::Matrix{Float64}, background::Float64, noisy::Matrix{Float64})

Get synthetic satellite retrival where noise is dependent on enhancement. 
"""
function oldenhancementdependentretrival!(i::Int64, sigma::Float64, data_ppm::Matrix{Float64}, background::Float64, noisy::Matrix{Float64})
    randn!(MersenneTwister(i), noisy)
    @. noisy = sigma * noisy
    @. noisy = (data_ppm + background) * (noisy + 1) - background
    return noisy
end


"""
    enhancementdependentretrival!(sigma::Float64, data_ppm::Matrix{Float64}, background::Float64, noisy::Matrix{Float64})

Get synthetic satellite retrival where noise is dependent on enhancement. 
"""
function enhancementdependentretrival!(i::Int64, sigma::Matrix{Float64}, data_ppm::Matrix{Float64}, noisy::Matrix{Float64})
    randn!(MersenneTwister(i), noisy)
    @. noisy = data_ppm + (sigma * noisy)
end


"""
    synthetic_satellite_data_1(noi_p, data_ppm, back)

Get synthetic satellite retrival where noise is independent of enhancement. 
"""
function enhancementindependentretrival!(sigma::Float64, data_ppm::Matrix{Float64}, noisy::Matrix{Float64})
    randn!(noisy)
    @. noisy = sigma * noisy
    @. noisy = data_ppm + noisy
    return noisy
end


"""
    sigmanoise_perpixel(background_ppm::Float64, sigmanoise::Float64, ymod::Vector{Float64})

TBW
"""
function sigmanoise_perpixel(background_ppm::Float64, sigmanoise::Float64, ymod::Vector{Float64})
    noise = ((ymod .+ background_ppm) .* sigmanoise)
    return noise
end




"""
        compute_perpixel_sigma()
        Computed the sigma error per pixel
"""
function compute_perpixelsigma!(level2precision::Float64, background, no2_to_co2::NO2toCO2,
                                distancefromsource::Matrix{Float64}, data::EmissionData)
    @. data.perpixelsigma = (data.conc + background) * level2precision
    if no2_to_co2.flag
        sigmaco2_no2 = NO2Error.get_co2_noise_from_no2(data.conc, distancefromsource, data.plume, no2_to_co2)
        data.perpixelsigma = NO2Error.inversevarianceweight_sigma(data.perpixelsigma, sigmaco2_no2)
    end
end



function dataupdateatemission!(resdata::ResolutionData, emission::Float64, data::EmissionData, plumethreshold=1e-6)
    data.emission = emission
    @. data.conc = resdata.conc * emission
    @. data.plume = data.conc > plumethreshold
end


end # end module
