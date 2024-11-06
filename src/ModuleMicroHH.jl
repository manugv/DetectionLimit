# Recode to convert microhh code to data that can be read in detection limit algorithm.

module MicroHH

using CSV: read as csvread
using DataFrames: DataFrame
using IniFile
using NCDatasets
using ..Datacontainer: Grid3d, ResolutionData, RectGrid2d, Grid2d
using ..ExtractDataatResolution
using ..Constants
using ..TOMLData: TOML, get_var
using HDF5


mutable struct MicroHHData
    const datadirectory::String
    const source::Vector{Float64}
    const emissionstrength::Float64
    const grid::Grid3d
    time::Int64
    conc::Matrix{Float64}
    convolved_conc::Matrix{Float64}    
end

mutable struct MicroHHDataVel
    const datadirectory::String
    const source::Vector{Float64}
    const emissionstrength::Float64
    const grid::Grid3d
    time::Int64
    conc::Matrix{Float64}
    vel::Matrix{Float64}
    vel_height::Float64
end


function directory_simtimes(filename)
    df = TOML.parsefile(filename)
    dd = get_var(df, "MicroHHData")
    dr = get_var(dd, "inputdatafilename")
    _file_time = get_var(dd, "time")
    # Get array of times corresponding to plumes
    simtimes = csvread(_file_time, DataFrame)
    return Fileandkeys(dr, simtimes.time)
end


"""
    _get_simulation_parameters(datadir, gas_mol_mass, scale=1)

TBW
"""
function _get_simulation_parameters(datadir::String, gas_mol_mass::Float64, scale::Float64)
    filename = filter(x -> endswith(x, ".ini"), readdir(datadir))
    ini = read(Inifile(), datadir * filename[1])
    sh = [parse(Int, get(ini, "grid", "itot")),
        parse(Int, get(ini, "grid", "jtot")),
        parse(Int, get(ini, "grid", "ktot"))]
    ds = [parse(Float64, get(ini, "grid", "xsize")) / sh[1],
        parse(Float64, get(ini, "grid", "ysize")) / sh[2],
        parse(Float64, get(ini, "grid", "zsize")) / sh[3]]
    # emissions in kilo moles
    emis = parse(Float64, get(ini, "source", "strength"))
    emis *= 1000  # kilomoles to moles
    # convert moles/sec to gms/sec to MT/year
    emis_mtpy = scale * emis * gas_mol_mass * Constants.gmspersec_2_MTperyear
    src = [parse(Float64, get(ini, "source", "source_x0")),
        parse(Float64, get(ini, "source", "source_y0")),
        parse(Float64, get(ini, "source", "source_z0"))]
    return src, emis_mtpy, Grid3d(sh, ds)
end


"""
    getsimulationparameters(datadir::String, gas_mol_mass::Float64, scale::Float64=1.0)

Get simulation parameters for MicroHH. Remain fixed for a simulation.
    gas_mol_mass corresponds to the molar mass of the gas. This is for any gas.
"""
function getsimulationparameters(datadir::String, gas_mol_mass::Float64, scale::Float64=1.0)
    src, emis_mtpy, _grid = _get_simulation_parameters(datadir, gas_mol_mass, scale)
    return MicroHHData(datadir, src, emis_mtpy, _grid, 0, zeros(_grid.nx, _grid.ny), zeros(_grid.nx, _grid.ny))
end



"""
    getsimulationparameters_cfm(datadir::String, gas_mol_mass::Float64, scale::Float64=1.0)

Get simulation parameters for MicroHH for CFM simulations. Remain fixed for a simulation.
        gas_mol_mass corresponds to the molar mass of the gas. This is for any gas.
"""
function getsimulationparameters_cfm(datadir::String, gas_mol_mass::Float64, scale::Float64=1.0)
    src, emis_mtpy, _grid = _get_simulation_parameters(datadir, gas_mol_mass, scale)
    return MicroHHDataVel(datadir, src, emis_mtpy, _grid, 0, zeros(_grid.nx, _grid.ny), zeros(_grid.nx, _grid.ny), src[3])
end


"""
    get_variables_microhh(data_dir, file_prefix, time, sh)

TBW
"""
function get_variables_microhh(data_dir::String, file_prefix::String, time::Int64, sh)
    data = Array{Float64}(undef, sh[1] * sh[2] * sh[3])
    filename = data_dir * file_prefix * "." * lpad(time, 7, '0')
    read!(filename, data)
    return reshape(data, (sh[1], sh[2], sh[3]))
end


"""
    get_files_for_rho(dir_path)

TBW
"""
function get_files_for_rho(dir_path)
    return filter(x -> contains(x, "default"), readdir(dir_path))
end


"""
    get_rho(path, flzs, time)

TBW
"""
function get_rho(path, flzs, time)
    for fl in flzs
        ff = NCDataset(path * fl)
        # this is to supress the warning that one gets while defining the time
        tmp = cfvariable(ff, "time", units="second")[:]
        idx = findfirst(isapprox.(tmp, time))
        if !isnothing(idx)
            rho = ff.group["thermo"]["rho"][:, idx]
            close(ff)
            return rho
            break
        else
            close(ff)
        end
    end
    return -99
end


"""
    convert_conc_to_2d(rho, dz, conc3d_ppm)

TBW
"""
function convert_conc_to_2d(rho, dz, conc3d_ppm)
    # rho in kg/m^3
    # molar mass of air in kg/mol
    molar_mass_air = Constants.air.molarmass * 0.001
    sh = size(conc3d_ppm)
    conc3d_moles = zeros(sh)
    air_moles = zeros(sh[3])
    for i = 1:sh[3]
        # Air conc in moles/m^3
        air_moles[i] = rho[i] / molar_mass_air
        conc3d_moles[:, :, i] = conc3d_ppm[:, :, i] .* air_moles[i]
    end
    conc2d_moles = sum(conc3d_moles, dims=3)[:, :] .* dz
    return conc2d_moles
end


"""
    microhh_get_conc(data_dir, gas_prefix, time, simgrid, scale_emission)

Get concentration in moles from MicroHH
"""
function microhh_get_conc(data_dir, gas_prefix, time, simgrid)
    conc3d_ppm = get_variables_microhh(data_dir, gas_prefix, time, [simgrid.nx, simgrid.ny, simgrid.nz])
    rho_files = get_files_for_rho(data_dir)
    rho = get_rho(data_dir, rho_files, time)
    conc_moles = convert_conc_to_2d(rho, simgrid.dz, conc3d_ppm)
    conc3d_ppm = nothing
    # scale emission
    return conc_moles
end


"""
    get_loc_factor(x::Vector{Float64}, x1::Float64)

Factors for interpolation
"""
function get_loc_factor(x::Vector{Float64}, x1::Float64)
    ii = argmin(abs.(x .- x1))
    dx = x[ii+1] - x[ii]
    ds = (x1 - x[ii]) / dx
    return ii, ds
end


"""
    get_microhh_velocity(data_dir, simgrid, time, vel_prefix)

Get velocity microHH
"""
function microhh_get_velocity(data_dir, grid, time, vel_prefix, vel_height)
    u = get_variables_microhh(data_dir, vel_prefix, time, [grid.nx, grid.ny, grid.nz])
    # Interpolate the velocity at a given height in z
    ii, ds = get_loc_factor(collect(grid.z), vel_height)
    vel = (u[:, :, ii] .* (1.0 - ds)) + (u[:, :, ii+1] .* ds)
    return vel
end


"""
    get_data!(data::MicroHHData)

Get microHH concentration
"""
function get_data!(data::MicroHHData, time, gas)
    data.time = time
    conc_moles = microhh_get_conc(data.datadirectory, "co2", data.time, data.grid)
    # convert moles to gms to ppm
    data.conc = (conc_moles * gas.molarmass) / gas.ppm_to_gms
end


"""
    get_concdata(filename, time)

Get microHH concentration
"""
function get_concdata(filename, time)
    fid = h5open(filename, "r")
    conc = read(fid[string(time)])  # in ppm
    emissionstrength = read_attribute(fid[string(time)], "strength")
    x = read(fid["grid/x"])
    y = read(fid["grid/y"])
    xn = read(fid["grid/xn"])
    yn = read(fid["grid/yn"])
    grid = Grid2d(x, y, xn, yn, x[2]-x[1], y[2]-y[1], length(x), length(y))
    source = read(fid["source"])
    return ConcData(source, emissionstrength, grid, time, conc, zeros(size(conc)))
end


# """
#     get_data!(data::MicroHHData)

# Get microHH concentration
# """
# function get_data!(data::MicroHHData, time, gas, no2_to_co2)
#     data.time = time
#     conc_moles = microhh_get_conc(data.datadirectory, "co2", data.time, data.grid)
#     # convert moles to gms to ppm
#     dist = ExtractDataatResolution.distancefromsource(data.grid.x, data.grid.y, data.source)
#     conc = (conc_moles * gas.molarmass) / gas.ppm_to_gms
    
# end


"""
    get_data!(data::MicroHHDataVel)

Get microHH concentration and velocity
"""
function get_data!(data::MicroHHDataVel, time, gas)
    data.time = time
    conc_moles = microhh_get_conc(data.datadirectory, "co2", data.time, data.grid)
    # convert moles to gms to ppm
    data.conc = (conc_moles * gas.molarmass) / gas.ppm_to_gms
    data.vel = microhh_get_velocity(data.datadirectory, data.grid, data.time, "u", data.vel_height)
end


function samplegridatresolution(microhhdata, resolution::Float64, xshift=1, yshift=1)
    # grid data
    # sample the data
    ix, xc, xn = ExtractDataatResolution.indicesatresolution(microhhdata.grid.x, microhhdata.grid.dx, resolution, xshift)
    iy, yc, yn = ExtractDataatResolution.indicesatresolution(microhhdata.grid.y, microhhdata.grid.dy, resolution, yshift)
    grid2d = RectGrid2d(xc, yc, xn, yn, resolution, 1.0 / resolution, length(xc), length(yc))
    # find the origin id
    origin_id = [findmin(abs.(xc .- microhhdata.source[1]))[2], findmin(abs.(yc .- microhhdata.source[2]))[2]]
    # the distance array 
    dist = ExtractDataatResolution.distancefromsource(xc, yc, microhhdata.source)
    return grid2d, origin_id, dist, ix, iy
end


# function samplegridatresolution(microhhdata, resolution::Float64, xshift=1, yshift=1)
#     # grid data
#     # sample the data
#     _, xc, xn = ExtractDataatResolution.indicesatresolution1(microhhdata.grid.x, microhhdata.grid.dx, resolution,
#                                                              microhhdata.source[1], xshift)
#     _, yc, yn = ExtractDataatResolution.indicesatresolution1(microhhdata.grid.y, microhhdata.grid.dy, resolution,
#                                                              microhhdata.source[2], yshift)
#     grid2d = RectGrid2d(xc, yc, xn, yn, resolution, 1.0 / resolution, length(xc), length(yc))
#     # find the origin id
#     origin_id = [findmin(abs.(xc .- microhhdata.source[1]))[2], findmin(abs.(yc .- microhhdata.source[2]))[2]]
#     # the distance array 
#     dist = ExtractDataatResolution.distancefromsource(xc, yc, microhhdata.source)
#     return grid2d, origin_id, dist
# end


"""
    convolveandsample(resolution::Float64, data::MicroHHData)

TBW
"""
# function convolvedata(resolution::Float64, dx::Float64, conc::Matrix{Float64})
#     # convolve image
#     fwhm = resolution / dx
#     convolved, _ = ExtractDataatResolution.convolveimage(conc, fwhm)
#     # compute convolved conc
#     return convolved[xshift:ix:end, yshift:iy:end]
# end


"""
    convolveandsample(resolution::Float64, data::MicroHHData)

TBW
"""
function convolveandextract(resolution::Float64, dx::Float64, conc::Matrix{Float64}, xshift=1, yshift=1)
    # convolve image
    fwhm = resolution / dx
    ix = Int(resolution/dx)
    iy = Int(resolution/dx)
    convolved, _ = ExtractDataatResolution.convolveimage(conc, fwhm)
    # compute convolved conc
    return convolved[xshift:ix:end, yshift:iy:end]
end


"""
    sampleandconvolve(resolution::Float64, data::MicroHHData)

    conc corresponds to emission 1 
"""
function defineresolutiondata(resolution::Float64, microhhdata, xshift=1, yshift=1)
    grid2d, origin_id, dist = samplegridatresolution(microhhdata, resolution, xshift, yshift)
    conc = zeros(size(dist))
    return ResolutionData(grid2d, origin_id, dist, conc)
end


"""
    sampleandconvolve(resolution::Float64, data::MicroHHData)

    conc corresponds to emission 1 
"""
function convolveandscale!(resolution::Float64, microhhdata, resdata::ResolutionData, xshift=1, yshift=1)
    resdata.conc = convolveandextract(resolution, microhhdata.grid.dx, microhhdata.conc, xshift, yshift)
    resdata.conc *= 1.0/microhhdata.emissionstrength
end


"""
    gridsandconvolve(resolution::Float64, data::MicroHHData)

    conc corresponds to emission 1 
"""
function gridsandextract(resolution::Float64, microhhdata, xshift=1, yshift=1)
    grid2d, origin_id, dist, ix, iy = samplegridatresolution(microhhdata, resolution, xshift, yshift)
    conc = zeros(size(dist))
    resdata = ResolutionData(grid2d, origin_id, dist, conc)
    resdata.conc = microhhdata.convolved_conc[xshift:ix:end, yshift:iy:end]
    resdata.conc *= 1.0/microhhdata.emissionstrength
    return resdata
end


"""
    gridsandconvolve(resolution::Float64, data::MicroHHData)

    conc corresponds to emission 1 
"""
function gridsandextract_mean(resolution::Float64, microhhdata, xshift=1, yshift=1)
    grid2d, origin_id, dist, ix, iy = samplegridatresolution(microhhdata, resolution, xshift, yshift)
    conc = zeros(size(dist))
    resdata = ResolutionData(grid2d, origin_id, dist, conc)
    if resolution == microhhdata.grid.dx
        resdata.conc = microhhdata.conc
    else
        resdata.conc = ExtractDataatResolution.convolveimagemean(microhhdata.conc, ix, iy, xshift, yshift)
    end
    resdata.conc *= 1.0/microhhdata.emissionstrength
    return resdata
end

  
end # end module
