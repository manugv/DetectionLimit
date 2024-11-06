if isdefined(@__MODULE__, :LanguageServer)
    include("../src/Tango.jl")
    using .Tango
else
    using Tango
end

using Interpolations
using Statistics: mean
using CSV, DataFrames
using HDF5, JLD2

###############################################################################
#                       some functions                                        #
###############################################################################
function get_co2interp()
    resolution = [60, 150, 300, 400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
    init_emis = [0.35, 0.9, 1.6, 2.4, 3.3, 4.0, 5.2, 6.4, 7.6, 9.7, 11., 12.2, 13.4]
    delta_emis = [0.0025, 0.005, 0.0075, 0.01, 0.01, 0.02, 0.02, 0.02, 0.04, 0.04, 0.06, 0.08, 0.1]
    fun_initialemission = linear_interpolation(resolution, init_emis)
    fun_deltaemission = linear_interpolation(resolution, delta_emis)
    return fun_deltaemission, fun_initialemission
end
    
function get_groupname(simname, res)
    return simname * "/Res_" * string(Int(res))
end


function compute_shifts(dx, dy, res)
    shx = Int(ceil(res/dx))
    shy = Int(ceil(res/dy))
    x1 = (1:shx)' .* ones(Int64, shy)
    y1 = ones(Int64, shx)' .* (1:shy)
    return hcat(vec(x1), vec(y1)), length(vec(x1))
end


function interpolate_for_detectionprob(prob, emis, detect_prob)
    k = searchsortedfirst(prob, detect_prob)
    if k == 1
        k1 = k
    else
        k1 = k - 1
    end
    k2 = k + 1
    if k2 > length(prob)
        k2 = length(prob)
    end
    fun_lin = LinearInterpolation(prob[k1:k2], emis[k1:k2])
    return fun_lin(detect_prob)
end

###############################################################################
#                                 main program                                #
###############################################################################
tomlfilename = "input.toml"

# Create interpolation functions
fun_deltaemission, fun_initialemission = get_co2interp()

# simulation parameters
simparams, detectionlimitparams = Tango.TOMLData.getparameters(tomlfilename)

# read data directory and time
dataparams = Tango.TOMLData.directory_simtimes(tomlfilename)

# level 2 precision change if needed
lvl2precision = 0.005
if ~detectionlimitparams.cdfflag
    detectionlimitparams.threshold = lvl2precision*simparams.gas.background
end

# read microhh data for a given time
concdata = Tango.SimData.get_concdata(dataparams.filename, dataparams.timekeys[1])

# write attributes to simulation group
attr_names = ["samples", "level2precision", "threshold", "numberofpixels", "plumelength", "cfdflag"]
attr_values = [detectionlimitparams.samples, lvl2precision, detectionlimitparams.threshold,
               detectionlimitparams.numberofpixels, detectionlimitparams.plumelength,
               detectionlimitparams.cdfflag]
Tango.ModReadWrite.write_grp_attributedata(simparams.outputfile, simparams.simulationname, attr_names, attr_values)


# pick a resolution
i = 5 # 300m
res = simparams.resolution[i]
println("   ", res,  " ", lvl2precision, "  ", detectionlimitparams.threshold)

# convolve data at a resolution
fwhm = res / concdata.grid.dx     # full mean half width
concdata.convolved_conc, _ = Tango.SimData.convolveimage(concdata.conc, fwhm)

# Define initial variables to compute the detection limit
init_emission = fun_initialemission(res)    # Initial emission
deltaemis = fun_deltaemission(res)          # Delta emission

# Data Sampling at a resolution
# compute variables defining all possible samples
shifts, nos_ensembles = compute_shifts(concdata.grid.dx, concdata.grid.dy, res)

# loop over all ensambles : for simplicity do one ensamble
ensemb = 1
println("                ensemble= ", ensemb)  

# Sample data and extract data at a resolution. This is actual data.
resdata = Tango.SimData.gridsandextract(res, concdata, shifts[ensemb,1], shifts[ensemb,2])
# Create a data structure to compute the detection limit.
detectdata = Tango.Datacontainer.DetectionLimitData(resdata)

# Compute few detection limits for a given range of detection probability
# range used here is 0.675 - 0.685 and minimal of 4 points are needed.
# Writing data: define group name 
detectdata.groupname = get_groupname(simparams.simulationname, res)

# initialize initial emission
detectdata.init_emission = init_emission

# Compute detection limits for the detection probability range.
Tango.Detectionlimit.getdatafornoise(simparams, detectionlimitparams, lvl2precision,
                                     detectdata, [0.675, 0.685], deltaemis)

# check if there are 4 points in the detection probability range.
# if not increase the sampling of emissions.
nos1 = sum((detectdata.probabilities .>= 0.675) .& (detectdata.probabilities .<= 0.685))
kl = 1
while nos1 < 4
    global deltaemis, nos1, kl
    detectdata.init_emission = interpolate_for_detectionprob(detectdata.probabilities, detectdata.allemissions, 0.68)
    if nos1 < 2
        deltaemis *= 0.2
    else
        deltaemis *= minimum([0.75, nos1/5.0])
    end
    detectdata.allemissions = []
    detectdata.probabilities = []
    Tango.Detectionlimit.getdatafornoise(simparams, detectionlimitparams, lvl2precision,
                                         detectdata, [0.675, 0.685], deltaemis)
    nos1 = sum((detectdata.probabilities .>= 0.675) .& (detectdata.probabilities .<= 0.685))
    kl += 1
    if kl > 4
        deltaemis *= 2
        break
    end
end
if nos1 > 8
    deltaemis *= maximum([1.2, nos1/8.0])
end

# Write data
Tango.ModReadWrite.write_output_data(simparams.outputfile, detectdata.groupname,
                                     "emission_"*string(ensemb), detectdata.allemissions)
Tango.ModReadWrite.write_output_data(simparams.outputfile, detectdata.groupname,
                                     "probabilities_"*string(ensemb), detectdata.probabilities)

# update the init_emission for faster computation in the sampling ensembles.
init_emission = interpolate_for_detectionprob(detectdata.probabilities, detectdata.allemissions, 0.68)
