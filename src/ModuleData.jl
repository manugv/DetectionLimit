"""
     Module that consists of functions to read the data from a hdf5 file.
"""

module SimData
using ImageFiltering: Kernel, imfilter, BorderArray, Pad
using Statistics: mean
using ..Datacontainer: Grid2d, ResolutionData, RectGrid2d, ConcData, Fileandkeys
using HDF5


# CONVOLUTION
"""
    fwhm2sigma(fwhm::Float64)

Compute sigma from full width and half mean
"""
function fwhm2sigma(fwhm::Float64)
    return fwhm / (2 * sqrt(2 * log(2)))
end


"""
    compute_convolution(data, fwhm, fill_value)

Compute 2d Gaussian convolution for given fwhm
"""
function convolveimage(data, fwhm)
    # get the Gaussian kernel
    ker = Kernel.gaussian(fwhm2sigma(fwhm))
    # get the convolved image
    # here border is replicated (default behaviour) as aaa|abcd|ddd
    convolved = imfilter(data, ker)
    return convolved, ker
end


"""
    get_concdata(filename, time)

Get concentration from given file and its key
"""
function get_concdata(filename::String, _key::String)
    if ~isfile(filename)
        println("Data file not present")
    end        
    fid = h5open(filename, "r")
    conc = read(fid[_key])  # in ppm
    emissionstrength = read_attribute(fid[_key], "strength")
    x = read(fid["grid/x"])
    y = read(fid["grid/y"])
    xn = read(fid["grid/xn"])
    yn = read(fid["grid/yn"])
    grid = Grid2d(x, y, xn, yn, x[2]-x[1], y[2]-y[1], length(x), length(y))
    source = read(fid["source"])
    return ConcData(source, emissionstrength, grid, time, conc, zeros(size(conc)))
end



# Sample and Extract data at a given resolution
"""
    get_indices_at_resolution(x, y, res)

    Get indices for a resolution from high resolution data.
"""
function indicesatresolution(x, dx, res, shift::Int64)
    ix = Int(res / dx)
    xc = x[shift:ix:end]
    xn = zeros(length(xc) + 1)
    xn[1:end-1] = xc .- (res * 0.5)
    xn[end] = xn[end-1] + res
    return ix, xc, xn
end


"""
    distancefromsource(x, y, oid)

Generate a distance map for all pixels from the origin.
"""
function distancefromsource(x, y, src)
    x1 = x .- src[1]
    y1 = y .- src[2]
    xv = getindex.(Iterators.product(x1, y1), 1)
    yv = getindex.(Iterators.product(x1, y1), 2)
    dist = sqrt.(xv .^ 2 + yv .^ 2)
    return dist
end


function samplegridatresolution(conc_data, resolution::Float64, xshift=1, yshift=1)
    # grid data
    # sample the data
    ix, xc, xn = indicesatresolution(conc_data.grid.x, conc_data.grid.dx, resolution, xshift)
    iy, yc, yn = indicesatresolution(conc_data.grid.y, conc_data.grid.dy, resolution, yshift)
    grid2d = RectGrid2d(xc, yc, xn, yn, resolution, 1.0 / resolution, length(xc), length(yc))
    # find the origin id
    origin_id = [findmin(abs.(xc .- conc_data.source[1]))[2], findmin(abs.(yc .- conc_data.source[2]))[2]]
    # the distance array 
    dist = distancefromsource(xc, yc, conc_data.source)
    return grid2d, origin_id, dist, ix, iy
end

"""
    gridsandconvolve(resolution::Float64, data::MicroHHData)

    conc corresponds to emission 1 
"""
function gridsandextract(resolution::Float64, conc_data, xshift=1, yshift=1)
    grid2d, origin_id, dist, ix, iy = samplegridatresolution(conc_data, resolution, xshift, yshift)
    conc = zeros(size(dist))
    resdata = ResolutionData(grid2d, origin_id, dist, conc)
    resdata.conc = conc_data.convolved_conc[xshift:ix:end, yshift:iy:end]
    resdata.conc *= 1.0/conc_data.emissionstrength
    return resdata
end


end
