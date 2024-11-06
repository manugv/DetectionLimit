if isdefined(@__MODULE__, :LanguageServer)
    include("../src/Tango.jl")
    using .Tango
else
    using Tango
end
using Plots

tomlfilename = "input_lsq.toml"

# simulation parameters
simparams, lsqparams = Tango.InputData.getparameters(tomlfilename)

# read microhh data directory and time
microhh_datadir = Tango.MicroHH.directory_simtimes(tomlfilename)

# read data
scale = 1.167759512561433
microhhdata = Tango.MicroHH.getsimulationparameters(microhh_params.datadir, simparams.gas.molarmass, scale)

p = plot()
for _time in microhh_params.filetimes
    println(_time)
    Tango.MicroHH.get_data!(microhhdata, _time, simparams.gas)
    resdata = Tango.MicroHH.convolveandsample(simparams.resolution[2], microhhdata);
    emission, lsq = Tango.LSQ.flux_at_level4precision(resdata, simparams, lsqparams, 0.005, 1)
    plot!(p, emission, (100.0.*lsq./emission), label=false)
end

display(p)

# for _time in microhh_params.filetimes
#     Tango.MicroHH.get_data!(microhhdata, _time, simparams.gas)
#     resdata = Tango.MicroHH.convolveandsample(simparams.resolution[2], microhhdata);
#     notfound = true
#     i = 7
#     while notfound
#         th = partialsort(vec(resdata.conc), i, rev=true)
#         plumemask = resdata.conc .>= th
#         segmentimg = zeros(Int64, size(resdata.conc))
#         actualplume = resdata.conc .>= 1e-6
#         Tango.Plume.segmentimage!(segmentimg, plumemask)
#         goodpixels, badpixels, plumelabel = Tango.Plume.plumepixelid(segmentimg, resdata.originindex, actualplume)
#         if goodpixels >= 7
#             plume = zeros(Bool, size(segmentimg))
#             plume = segmentimg .== plumelabel
#             #p = heatmap(resdata.grid.xn, resdata.grid.yn, plume')
#             #display(p)
#             notfound = false
#             println(_time,"   ", i)
#         end
#         i = i + 1
#     end
# end
