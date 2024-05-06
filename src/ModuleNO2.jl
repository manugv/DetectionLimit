module NO2Error

using ..Datacontainer: ResolutionData, NO2toCO2
using ..Constants


f(x) = (3.0528241104694387 * exp(-5.921587247929098 * x)
        + (0.010106636251216891 / x) + 1.315778787404395)

f1(x) = (-0.5269555789039001 * exp(-3.7262977658590692 * x)
         -
         (0.0001555336896657291 / x) + 0.012092214716475733 * x
         + 0.7306700581683978)


"""
    inversevarianceweighting(sigma1::Vector{Float64}, sigma2::Vector{Float64})

TBW
"""
function inversevarianceweighting(sigma1::Vector{Float64}, sigma2::Vector{Float64})
    var1 = 1.0 ./ (sigma1 .^ 2)
    var2 = 1.0 ./ (sigma2 .^ 2)
    return @. (1.0 / (var1 + var2))
end


"""
    inversevarianceweighting_sigma(sigma1::Matrix{Float64}, sigma2::Matrix{Float64})

TBW
"""
function inversevarianceweighting_sigma(sigma1::Matrix{Float64}, sigma2::Matrix{Float64})
    var1 = 1.0 ./ (sigma1 .^ 2)
    var2 = 1.0 ./ (sigma2 .^ 2)
    return @. sqrt(1.0 / (var1 + var2))
end


inversevarianceweight_sigma(x, y) = @. sqrt(1.0/((1.0/x^2) + (1.0/y^2)))



"""
    get_transform_matrix(dist)

dist is given in mts
"""
function get_transform_matrix(dist::Matrix{Float64})
    nox_to_no2 = f1.(dist / 10000)
    nox_to_no2[isinf.(nox_to_no2)] .= 0
    return nox_to_no2
end


"""
    get_co2_noise_from_no2(conc::Matrix{Float64}, distancefromsource::Matrix{Float64}, plume::BitMatrix, no2_co2_factors::NO2toCO2)

TBW
"""
function get_co2_noise_from_no2(conc::Matrix{Float64}, distancefromsource::Matrix{Float64}, 
    plume::BitMatrix, no2_co2_factors::NO2toCO2)
    # molecules of NO2 per gram 
    no2_molecules_per_gm = Constants.avogadro / Constants.no2.molarmass
    # get nox gms/m^2 to no2 molecules/m^2 conversion
    nox_to_no2 = get_transform_matrix(distancefromsource)
    @. nox_to_no2 *= no2_molecules_per_gm
    # create a matrix for all computations
    _data = similar(conc)
    # convert ppm to gms to NOx
    @. _data = conc * (Constants.co2.ppm_to_gms * no2_co2_factors.co2_to_nox)
    # convert NOx gms/m^2 to NO2 molecules/m^2
    @. _data *= nox_to_no2
    # Compute relative error in NO2
    @. _data *= no2_co2_factors.relative
    # Convert NO2 molecules/m^2 error to error in NOx in gms in the 
    # convert to NOX and add absolute error in plume and then add non-plume absolute error as constant
    @. _data = ((_data + no2_co2_factors.absolute * plume) / nox_to_no2) +
               (( no2_co2_factors.absolute * ~plume * 1.35) / no2_molecules_per_gm)
    @. _data[isinf(_data)] = 0.0
    # Convert NOx gms/m^2 to CO2 and then to ppm
    @. _data *= (1.0 / (Constants.co2.ppm_to_gms * no2_co2_factors.co2_to_nox))
    return _data
end


"""
    getdata_no2(conc::Matrix{Float64}, distancefromsource::Matrix{Float64}, plume::BitMatrix, no2_co2_factors::NO2toCO2)

TBW
"""
function getdata_no2(conc::Matrix{Float64}, distancefromsource::Matrix{Float64}, no2_co2_factors::NO2toCO2)
    # Create a conversion factor from NOX g/m^2 to NO2 molecules/m^2
    no2_molecules_per_gm = Constants.avogadro / Constants.no2.molarmass
    # get nox gms/m^2 to no2 molecules/m^2 conversion
    nox_to_no2 = get_transform_matrix(distancefromsource)
    @. nox_to_no2 *= no2_molecules_per_gm

    # create a matrix for all computations
    _data = similar(conc)
    # convert CO2 ppm to gms and then to NOx
    @. _data = conc * (Constants.co2.ppm_to_gms * no2_co2_factors.co2_to_nox)
    # convert NOx gms/m^2 to NO2 molecules/m^2
    @. _data *= nox_to_no2
    @. _data *= (1.0/no2_molecules_per_gm)
    @. _data[isinf(_data)] = 0.0
    return _data
end


end  # end module
