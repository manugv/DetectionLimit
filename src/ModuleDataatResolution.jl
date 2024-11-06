# Flag TO BE REMOVED.

module ExtractDataatResolution

using ImageFiltering: Kernel, imfilter, BorderArray, Pad
using Statistics: mean

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


# convolve mean
ff(i, j, x, nx, ny) = mean(x[i-nx:i+nx, j-ny:j+ny])
"""
This is only valid for uniform distance in x and y
"""
function convolveimagemean(conc, nx, ny, sx, sy)
    sh = size(conc)
    halfnx = Int((nx - 1) / 2)
    halfny = Int((ny - 1) / 2)
    ar1 = BorderArray(conc, Pad((halfnx, halfny)))
    g(i, j) = ff(i, j, ar1, halfnx, halfny)
    x1 = collect(sx:nx:sh[1])
    y1 = collect(sy:ny:sh[2])
    kx = ones(Int, length(y1))' .* x1
    ky = y1' .* ones(Int, length(x1))
    return g.(kx, ky)
end



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
    get_indices_at_resolution_test(x, y, res)

    Get indices for a resolution from high resolution data.
"""
function indicesatresolution1(x, dx, res, src, shift)
    sx = Int(res / dx)
    ix1 = argmin(abs.(x .- src))
    ix = [reverse(collect(ix1:-sx:1))[1:end-1]; collect(ix1:sx:length(x))]
    xc = x[ix]
    xn = collect(xc[1]-res*0.5:res:xc[end]+res*0.5)
    # xn = zeros(length(xc) + 1)
    # xn[1:end-1] = xc .- (res*0.5)
    # xn[end] = xn[end-1] + res
    return ix, xc, xn
end


"""
    get_origin_distance(x, y, oid)

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

end # module end
