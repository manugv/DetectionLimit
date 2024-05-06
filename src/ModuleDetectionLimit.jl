module Detectionlimit

using ..Plume
using ..Datacontainer: DetectionLimitData, ResolutionData, DetectionLimitVars, GeneralSimulationparameters, NO2toCO2
using ..SatelliteRetrieval: retrivalbasedonpixelsigma!, dataupdateatemission!, compute_perpixelsigma!
using ..NO2Error
using Statistics: mean


function getdelatemission(ds::Float64, level2precision::Float64, gas::String)
    if (ds <= 120.0)
        if (level2precision == 0.0005) || (level2precision == 0.001)
            delta_emis = 0.001
        elseif (level2precision <= 0.004) & (level2precision > 0.001)
            delta_emis = 0.01
        else
            delta_emis = 0.025
        end
    elseif (ds > 120.0) & (ds <= 600.0)
        if (level2precision == 0.0005) || (level2precision == 0.001)
            delta_emis = 0.005
        elseif (level2precision <= 0.004) & (level2precision > 0.001)
            delta_emis = 0.025
        elseif (level2precision > 0.004) & (level2precision <= 0.01)
            delta_emis = 0.05
        else
            delta_emis = 0.1
        end
    else
        if (level2precision == 0.0005) || (level2precision == 0.001)
            delta_emis = 0.01
        elseif (level2precision <= 0.004) & (level2precision > 0.001)
            delta_emis = 0.05
        elseif (level2precision > 0.004) & (level2precision < 0.01)
            delta_emis = 0.1
        elseif (level2precision >= 0.01) & (level2precision < 0.02)
            delta_emis = 0.2
        else
            delta_emis = 0.4
        end
    end

    if gas =="co2"
        return delta_emis
    elseif gas == "ch4"
        return delta_emis/1000
    end
end


function computeforward(emis::Float64, delta_emis::Float64, level2precision::Float64, simparams::GeneralSimulationparameters,
                        detectionlimitparams::DetectionLimitVars, data::DetectionLimitData, bound::Float64)
    emis_array = Float64[] 
    probability = Float64[]
    prob = 0.0
    while prob < bound
        # emission scale to 2.1
        dataupdateatemission!(data.actual, emis, data.emisdata)
        compute_perpixelsigma!(level2precision, simparams.gas.background, simparams.no2_to_co2,
                               data.actual.distancefromsource, data.emisdata)
        prob = computeprobability(data, detectionlimitparams)
        push!(probability, prob)
        push!(emis_array, emis)
        println("         Emission and probability: ", emis, "  ", prob)
        emis += delta_emis
    end
    return emis_array, probability
end


function computebackward(emis::Float64, delta_emis::Float64, level2precision::Float64,
                         simparams::GeneralSimulationparameters, detectionlimitparams::DetectionLimitVars,
                         data::DetectionLimitData, bound::Float64)
    emis_array = Float64[] 
    probability = Float64[]
    prob = 1.0
    while prob > bound
        # emission scale to 2.1
        dataupdateatemission!(data.actual, emis, data.emisdata)
        compute_perpixelsigma!(level2precision, simparams.gas.background, simparams.no2_to_co2,
                               data.actual.distancefromsource, data.emisdata)
        prob = computeprobability(data, detectionlimitparams)
        push!(probability, prob)
        push!(emis_array, emis)
        println("         Emission and probability: ", emis, "  ", prob)        
        emis -= delta_emis
    end
    return reverse(emis_array), reverse(probability)
end


"""
   computedetectionprobability(data::ResolutionData, params::DetectionLimitVars, satellite::SatelliteData)
            
Computed the probability to detect a plume for a given data and level 2 precision 
    """
function computeprobability(data::DetectionLimitData, params::DetectionLimitVars)
    # calculate plume detection flag
    flag = zeros(Bool, params.samples)
    # Parameters to detect plume. A mutable struct to detect plumes for different retrievals
    # run loops
    Threads.@threads for i = 1:params.samples
        id = Threads.threadid()
        retrivalbasedonpixelsigma!(i, data.emisdata.perpixelsigma, data.emisdata.conc, data.level2[id].conc)
        flag[i] = Plume.detection(data.level2[id], data.emisdata.plume, data.actual, params)
    end
    # return the detection probability 
    return mean(flag)
end


function getdatafornoise(simparams::GeneralSimulationparameters, detectionlimitparams::DetectionLimitVars,
                         level2precision::Float64, data::DetectionLimitData, bounds=[0.6, 0.995], delta_emis=0)
    # get delta value to hop
    if delta_emis == 0
        delta_emis = getdelatemission(data.actual.grid.ds, level2precision, simparams.gas.name)
    end
    emis = data.init_emission
    # forward simulations
    femis, fprob = computeforward(emis, delta_emis, level2precision, simparams, detectionlimitparams, data, bounds[2])
    # backward simulations if probability is > 0.6 
    if fprob[1] > bounds[1]
        emis = femis[1] - delta_emis 
        bemis, bprob = computebackward(emis, delta_emis, level2precision, simparams, detectionlimitparams, data, bounds[1])
        new_emission = vcat(bemis, femis)
        new_prob = vcat(bprob, fprob)
    else
        new_emission = femis
        new_prob = fprob
    end
    append!(data.allemissions, new_emission)
    append!(data.probabilities, new_prob)
end


# function optimize_emis_prob(emis, prob, det_prob, bound, deltaemis)
#     if length(prob) == 1
#         if prob[1] <= det_prob
#             em1 = emis[1] + deltaemis
#         else
#             em1 = emis[1] - deltaemis
#         end
#     else
#         l = sortperm(prob)
#         pb_sort = prob[l]
#         em_sort = emis[l]

#         if (prob[end] >= (det_prob - bound)) & (prob <= (det_prob + bound))
#         end
#     end
#     return em1
# end

# function optimize_datafornoise(simparams::GeneralSimulationparameters, detectionlimitparams::DetectionLimitVars,
#                                level2precision::Float64, data::DetectionLimitData, det_prob = 0.68,
#                                bound=0.005, delta_emis=0.0)
#     # get delta value to hop
#     if delta_emis == 0
#         delta_emis = getdelatemission(data.actual.grid.ds, level2precision, simparams.gas.name)
#     end
#     emis = data.init_emission

#     # simulations
#     emis_array = Float64[] 
#     probability = Float64[]
#     simcontinue = true
#     while simcontinue
#         dataupdateatemission!(data.actual, emis, data.emisdata)
#         compute_perpixelsigma!(level2precision, simparams.gas.background, simparams.no2_to_co2,
#                                data.actual.distancefromsource, data.emisdata)
#         prob = computeprobability(data, detectionlimitparams)
        
#         push!(probability, prob)
#         push!(emis_array, emis)

#         emis = optimize_emis_prob(emis_array, probability, det_prob, bound, delta_emis)
#         println("         Emission and probability: ", emis, "  ", prob)
#         emis += delta_emis
#     end
    
    
#     # backward simulations if probability is > 0.6 
#     if fprob[1] > bounds[1]
#         emis = femis[1] - delta_emis 
#         bemis, bprob = computebackward(emis, delta_emis, level2precision, simparams, detectionlimitparams, data, bounds[1])
#         new_emission = vcat(bemis, femis)
#         new_prob = vcat(bprob, fprob)
#     else
#         new_emission = femis
#         new_prob = fprob
#     end
#     append!(data.allemissions, new_emission)
#     append!(data.probabilities, new_prob)
# end



end # end module
