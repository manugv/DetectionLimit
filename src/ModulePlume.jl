module Plume

using StatsBase: fit, Histogram
using ImageMorphology: label_components!
using LinearAlgebra: normalize
using ..Datacontainer: ResolutionData, DetectionLimitVars, Level2Retrieval, EmissionData

"""
    cdfbasedthreshold(plume::Matrix{Float64}, thresholdweight::Float64)
    
Get threshold based on CDF
    """
function cdfbasedthreshold(plume::Matrix{Float64}, thresholdweight::Float64)
    # compute histogram
    h = fit(Histogram, vec(plume), nbins=100)
    h = normalize(h, mode=:probability)
    wt1 = cumsum(h.weights)
    ix = findfirst(wt1 .> thresholdweight)
    thresh = (0.5 * (h.edges[1][ix-1] + h.edges[1][ix]) +
              (thresholdweight - wt1[ix-1]) *
              (h.edges[1][ix] - h.edges[1][ix-1]) / (wt1[ix] - wt1[ix-1]))
    return thresh
end


"""
    thresholdbasedplumemask(res_data_ppm, threshold)
  Compute plumemask based on a threshold value.
"""
function thresholdbasedplumemask(res_data_ppm::Matrix{Float64}, threshold::Float64)
    if threshold == 0
        plumemask = trues(size(res_data_ppm))
    else
        plumemask = falses(size(res_data_ppm))
        @. plumemask = res_data_ppm > threshold
    end
    return plumemask
end


function thresholdbasedplumemask!(plumemask::BitMatrix, res_data_ppm::Matrix{Float64},
                                  threshold::Float64)
    if threshold == 0
        @. plumemask = true
    else
        @. plumemask = res_data_ppm > threshold
    end
end


"""
    cdfthresholdbasedplumemask(res_data_ppm, threshold)

Compute threshold based on cumulative distribution function. 
    And use the threshold to compute plumemask
"""
function cdfthresholdbasedplumemask(res_data_ppm::Matrix{Float64}, thresholdweight::Float64)
    if thresholdweight == 0
        plumemask = trues(size(res_data_ppm))
    else
        threshold = cdfbasedthreshold(res_data_ppm, thresholdweight)
        plumemask = thresholdbasedplumemask(res_data_ppm, threshold)
    end
    return plumemask
end


"""
    cdfthresholdbasedplumemask(res_data_ppm, threshold)

Compute threshold based on cumulative distribution function. 
    And use the threshold to compute plumemask
"""
function cdfthresholdbasedplumemask!(plumemask::BitMatrix, res_data_ppm::Matrix{Float64},
                                     thresholdweight::Float64)
    if thresholdweight == 0
        @. plumemask = true
    else
        threshold = cdfbasedthreshold(res_data_ppm, thresholdweight)
        thresholdbasedplumemask!(plumemask, res_data_ppm, threshold)
    end
end


"""
    segmentimage!(segmentedimage::Matrix{Int64}, plumemask::BitMatrix)

TBW
"""
function segmentimage!(segmentedimage::Matrix{Int64}, plumemask::BitMatrix)
    # label the connected regions
    _conn = trues(3, 3)
    label_components!(segmentedimage, plumemask, _conn)
end


"""
    checkplumelength!(segmentedimage::Matrix{Int64}, plumemask::BitMatrix, plumelabel::Int64,
    distancefromsource::Matrix{Float64}, plumelength::Float64)

TBW
"""
function checkplumelength!(segmentedimage::Matrix{Int64}, plumemask::BitMatrix, plumelabel::Int64,
                           distancefromsource::Matrix{Float64}, plumelength::Float64)
    @. plumemask = segmentedimage == plumelabel
    # Check for plume
    val = findmax(distancefromsource[plumemask])
    return val[1] >= plumelength
end


"""
    singleplumebasedoncdfthreshold(res_data_ppm::Matrix{Float64}, params, id::Vector{Int64}, actualplume::BitMatrix)

Compute threshold based on cumulative distribution function. 
    And use the threshold to compute plumemask
"""
function extractsingleplume(res_data_ppm::Matrix{Float64}, params,
                            id::Vector{Int64}, actualplume::BitMatrix)
    if params.cdfflag
        plumemask = cdfthresholdbasedplumemask(res_data_ppm, params.threshold)
    else
        plumemask = thresholdbasedplumemask(res_data_ppm, params.threshold)
    end
    segmentedimage = zeros(Int64, size(res_data_ppm))
    segmentimage!(segmentedimage, plumemask)
    goodpixels, badpixels, plumelabel = plumepixelid(segmentedimage, id, actualplume)
    @. plumemask = segmentedimage == plumelabel
    return plumemask
end


"""
    get_plumepixelid(plm::Matrix{Int64}, id::Vector{Int64})

Find the label coresponding to plume  
"""
function plumepixelid(segmentedimg::Matrix{Int64}, id::Vector{Int64}, actualplume::BitMatrix)
    i1 = max(1, id[1] - 2)
    i2 = id[1] + 2
    # get plume ids near the source
    lbs = unique(segmentedimg[i1:i2, id[2]-2:id[2]+2])
    # find the largest plume near source
    pixel_nos = 0
    plume_id = 0
    for lb in lbs
        if lb != 0
            tmp = count(==(lb), segmentedimg)
            if tmp > pixel_nos
                pixel_nos = tmp
                plume_id = lb
            end
        end
    end
    # find true and false plume pixels
    if pixel_nos > 0
        true_pixels = sum(actualplume[segmentedimg.==plume_id])
        false_pixels = pixel_nos - true_pixels
    else
        true_pixels = 0
        false_pixels = 0
    end
    return true_pixels, false_pixels, plume_id
end



"""
    check_ifplumeexists(conc_ppm::Matrix{Float64}, segmentedimage::Matrix{Int64},
    plumemask::Matrix{Bool}, id, param, dist)

Check for a plume based on conditions
"""
function detection(satellite::Level2Retrieval, actualplume::BitMatrix, data::ResolutionData, params::DetectionLimitVars)
    goodpixels = 0
    badpixels = 0
    if params.cdfflag
        # Mask pixels based on cdf threshold
        cdfthresholdbasedplumemask!(satellite.plumemask, satellite.conc, params.threshold)
    else
        thresholdbasedplumemask!(satellite.plumemask, satellite.conc, params.threshold)
    end
    # Label the connected regions of the mask
    segmentimage!(satellite.segmentedimage, satellite.plumemask)

    # Identify largest segmented region around source and find good & bad pixels and the segment id
    goodpixels, badpixels, plumelabel = plumepixelid(satellite.segmentedimage, data.originindex, actualplume)

    # if number of pixels are more then the threshold pixels
    filter1 = goodpixels >= params.numberofpixels

    # check the maximum distance of plume from origin if filter1 is true
    if filter1 & (params.plumelength > 0.0)
        filter2 = checkplumelength!(satellite.plumemask, satellite.segmentedimage, plumelabel,
                                    data.distancefromsource, params.plumelength)
        return (filter1 & filter2)
    else
        return filter1
    end
end

end # module end
