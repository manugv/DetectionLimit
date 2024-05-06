module ModReadWrite

using HDF5

# function get_groupname(root_grp, grp)
#     if grp in keys(root_grp)
#         return root_grp[grp]
#     else
#         return create_group(root_grp, grp)
#     end
# end

function get_groupname(root_grp, grps)
    _ky = string(grps[1])
    if haskey(root_grp, _ky)
        tmp = root_grp[_ky]
        if length(grps) == 1
            return tmp
        else
            return get_groupname(tmp, grps[2:end])
        end
    else
        return create_group(root_grp, join(grps, "/"))
    end
end


function write_output_data(fname, group_name, data_name, data_value, compress=false)
    h5open(fname, "cw") do fid
        grpnames = collect(eachsplit(group_name, "/"))
        g = get_groupname(fid, grpnames)
        if compress
            g[data_name, chunk=(1, 10000), shuffle=(), deflate=3] = data_value
        else
            g[data_name] = data_value
        end
    end
end


function write_output_data(fname, data_name, data_value)
    h5open(fname, "cw") do fid
        fid[data_name] = data_value
    end
end


function write_grp_attributedata(fname, group_name, attrsname, attrsvalue)
    h5open(fname, "cw") do fid
        grpnames = collect(eachsplit(group_name, "/"))
        g = get_groupname(fid, grpnames)
        for i in eachindex(attrsname, attrsvalue)
            attributes(g)[attrsname[i]] = attrsvalue[i]
        end
    end
end


function write_data_attributedata(fname, group_name, dataname, attrsname, attrsvalue)
    grpnames = collect(eachsplit(group_name, "/"))
    h5open(fname, "cw") do fid
        g = get_groupname(fid, grpnames)
        dset = g[dataname]
        for i in eachindex(attrsname, attrsvalue)
            attributes(dset)[attrsname[i]] = attrsvalue[i]
        end
    end
end


function write_data_attributedata(fname, dataname, attrsname, attrsvalue)
    h5open(fname, "cw") do fid
        dset = fid[dataname]
        for i in eachindex(attrsname, attrsvalue)
            attributes(dset)[attrsname[i]] = attrsvalue[i]
        end
    end
end


function get_existing_group(root_grp, grp)
    if grp in keys(root_grp)
        return root_grp[grp]
    else
        print("ERROR")
    end
end


function get_data(fname, group_name, data_name)
    fid = h5open(fname, "r")
    g = get_existing_group(fid, group_name)
    data_value = read(g, data_name)
    close(fid)
    return data_value
end

function get_data(fname, data_name)
    fid = h5open(fname, "r")
    data_value = read(fid, data_name)
    close(fid)
    return data_value
end


end # end module
