module MacroDetectionLimit

using ..SatelliteRetrieval: dataupdateatemission!, compute_perpixelsigma!
using ..Detectionlimit: computeprobability
using ..Datacontainer: GeneralSimulationparameters, DetectionLimitVars, ResolutionData, SatelliteData


function getdelatemission(ds::Float64, level2precision::Float64, gas::String)
    if (ds <= 120.0)
        if (level2precision <= 0.004) 
            delta_emis = 0.01
        else
            delta_emis = 0.025
        end
    else
        if (level2precision <= 0.004)
            delta_emis = 0.025
        else
            delta_emis = 0.05
        end
    end
    if gas =="co2"
        return delta_emis
    elseif gas == "ch4"
        return delta_emis/1000
    end
end


function computeforward(emis::Float64, delta_emis::Float64, level2precision::Float64, simparams::GeneralSimulationparameters,
                        detectionlimitparams::DetectionLimitVars, resdata::ResolutionData, satellitedata::SatelliteData, bound::Vector{Float64})
    emis_array = Float64[] 
    probability = Float64[]
    prob = 0.0
    while prob < bound
        # emission scale to 2.1
        dataupdateatemission!(resdata, emis, satellitedata)
        compute_perpixelsigma!(level2precision, simparams.gas.background, simparams.no2_to_co2,
                               resdata.distancefromsource, satellitedata)
        prob = computeprobability(resdata,detectionlimitparams, satellitedata)
        push!(probability, prob)
        push!(emis_array, emis)
        emis += delta_emis
    end
    return emis_array, probability
end


function computebackward(emis::Float64, delta_emis::Float64, level2precision::Float64, simparams::GeneralSimulationparameters,
                         detectionlimitparams::DetectionLimitVars, resdata::ResolutionData, satellitedata::SatelliteData, bound::Vector{Float64})
    emis_array = Float64[] 
    probability = Float64[]
    prob = 1.0
    while prob > bound
        # emission scale to 2.1
        dataupdateatemission!(resdata, emis, satellitedata)
        compute_perpixelsigma!(level2precision, simparams.gas.background, simparams.no2_to_co2,
                                                        resdata.distancefromsource, satellitedata)
        prob = computeprobability(resdata,detectionlimitparams, satellitedata)
        push!(probability, prob)
        push!(emis_array, emis)
        emis -= delta_emis
    end
    return reverse(emis_array), reverse(probability)
end

    
function getdatafornoise(simparams::GeneralSimulationparameters, detectionlimitparams::DetectionLimitVars, level2precision::Float64,
                         resdata::ResolutionData, satellitedata::SatelliteData, init_emis::Float64, bounds=[0.6, 0.995], delta_emis=0)
    # get delta value to hop
    if delta_emis == 0
        delta_emis = getdelatemission(resdata.grid.ds, level2precision, simparams.gas.name)
    end
    emis = init_emis
    # forward simulations
    femis, fprob = computeforward(emis, delta_emis, level2precision, simparams, detectionlimitparams, resdata, satellitedata, bounds[2])
    # backward simulations if probability is > 0.6 
    if fprob[1] > 0.6
        emis = femis[1] - delta_emis 
        bemis, bprob = computebackward(emis, delta_emis, level2precision, simparams, detectionlimitparams, resdata, satellitedata, bounds[1])
        new_emission = vcat(bemis, femis)
        new_prob = vcat(bprob, fprob)
    else
        new_emission = femis
        new_prob = fprob
    end
    append!(satellitedata.allemissions, new_emission)
    append!(satellitedata.probabilities, new_prob)
end

end # end module
