if isdefined(@__MODULE__, :LanguageServer)
    include("../src/Tango.jl")
    using .Tango
else
    using Tango
end
using Plots

# Input parameters file
tomlfilename = "input_lsq.toml"

# simulation parameters
simparams, lsqparams = Tango.TOMLData.getparameters(tomlfilename)
# read data directory and time
dataparams = Tango.TOMLData.directory_simtimes(tomlfilename)

plot_p = plot(xlabel = "Emission", ylabel = "Level 4 precision (%)")

# pick a resolution
res = simparams.resolution[2]
println("Resolution is ", res, "m and level-2 precision is ", simparams.level2precision[1]*100, "%")  # fix this
lsq_data = zeros(length(lsqparams.emission))

########## First loop over all plumes
for data_key in dataparams.timekeys
    global res, plot_p, lsq_data
    local lsq
    # read data at given key
    concdata = Tango.SimData.get_concdata(dataparams.filename, data_key)

    # convolve data at a resolution
    fwhm = res / concdata.grid.dx     # full mean half width
    concdata.convolved_conc, _ = Tango.SimData.convolveimage(concdata.conc, fwhm)

    # Data Sampling at a resolution
    # compute variables defining all possible samples
    shifts, nos_ensembles = Tango.SimData.compute_shifts(concdata.grid.dx, concdata.grid.dy, res)

    ######### Loop over different samples
    # for one ensemble
    lsq = zeros(length(lsqparams.emission), nos_ensembles)
    for ensemb = 1:nos_ensembles
        # Sample data and extract data at a resolution. This is actual data and is defined to 1 MT/y
        resdata = Tango.SimData.gridsandextract(res, concdata, shifts[ensemb,1], shifts[ensemb,2])
        # last value can either be "Connected" or "Top"
        # "Connected" means all connected pixels with source (defined in toml file)
        # "Top" means top number pixels (defined in toml file)
        lsq[:, ensemb] = Tango.LSQ.flux_at_level4precision(resdata, simparams, lsqparams, simparams.level2precision[1], "Top")
    end
    lsq_data .+= Tango.SimData.mean(lsq, dims=2)
    plot!(plot_p, lsqparams.emission, 100*Tango.SimData.mean(lsq, dims=2)./lsqparams.emission)
end
display(plot_p)
savefig(plot_p, "lsq_"*string(res)*"_th.png")

pp = plot(lsqparams.emission, 100*(lsq_data./length(dataparams.timekeys))./lsqparams.emission, xlabel = "Emission", ylabel = "Level 4 precision (%)")
savefig(pp, "lsq_final_"*string(res)*"_th.png")
