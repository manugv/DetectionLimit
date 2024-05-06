module InputData

using TOML
using ..Datacontainer
using ..Constants

################################################################################

"""
    get_var(df, var)

TBW
"""
function get_var(df, var)
    if var in keys(df)
        return df[var]
    else
        println(var * " variable not defined")
    end
end

################################################################################

"""
    get_res_lvl2(df)

TBW
"""
function get_res_lvl2(df)
    res = get_var(df, "resolution")
    lvl2p = get_var(df, "level2precision")
    return res, lvl2p
end


function get_no2_co2(df)
    var1 = get_var(df, "nox_to_co2")
    var2 = get_var(df, "absoluteerror")
    var3 = get_var(df, "relativeerror")
    var4 = get_var(df, "flag")
    return Datacontainer.NO2toCO2(var1, var2, var3, var4)
end


"""
    read_simulation_parameters(df)

Read simulation parameters
"""
function general_simulation_parameters(df1)
    df = get_var(df1, "GeneralSimulationParameters")
    # get the gas to be studied
    gasname = get_var(df, "gas")
    if gasname == "CO2"
        gas = Constants.co2
    elseif gasname == "NO2"
        gas = Constants.no2
    elseif gasname == "CH4"
        gas = Constants.ch4
    end
    # get simulation name
    simname = get_var(df, "simulationname")
    # get output file location
    outputfilename = get_var(df, "outputfile")
    # get resolution and level 2 precision
    res, lvl2 = get_res_lvl2(df)
    # Get no2_to_co2 conversion data
    no2 = get_var(df, "NO2toCO2")
    no2_to_co2 = get_no2_co2(no2)
    Datacontainer.GeneralSimulationparameters(simname, gas, outputfilename, res, lvl2, no2_to_co2)
end


#############################################################################
# function get_air_const(df)
#     gas_var = get_var(df, "AIR")
#     mass = get_var(gas_var, "mol_mass")
#     tot_mass = get_var(gas_var, "integrated_column_mass")
#     tot_airmass_moles = tot_mass / mass
#     return Datacontainer.airconst("air", mass, tot_airmass_moles)
# end


"""
    get_gas_const(df, gas, tot_airmass_moles)

TBW
"""
# function get_gas_const(df, gas, tot_airmass_moles)
#     gas_var = get_var(df, gas)
#     mol_mass = get_var(gas_var, "mol_mass")
#     back = get_var(gas_var, "background")
#     ppm_to_gms = tot_airmass_moles * mol_mass / 1e6
#     return Datacontainer.gasconst(gas, mol_mass, back, ppm_to_gms)
# end


"""
    get_gasconstants(df)

TBW
"""
# function get_gasconstants(df)
#     av = df["avogadro"]
#     air = get_air_const(df)
#     co2 = get_gas_const(df, "CO2", air.column_mass_moles)
#     ch4 = get_gas_const(df, "CH4", air.column_mass_moles)
#     no2 = get_gas_const(df, "NO2", air.column_mass_moles)
#     return Datacontainer.GasConstants(av, air, co2, ch4, no2)
# end
################################################################################


function get_lsqparams(df)
    v1 = get_var(df, "level4precision")
    v2 = get_var(df, "cdfflag")
    v3 = get_var(df, "threshold")
    return Datacontainer.LSQParams(v1, v2, v3)
end


function get_plumeparams(df)
    v1 = get_var(df, "cdfflag")
    v2 = get_var(df, "threshold")
    return Datacontainer.PlumeParams(v1, v2)
end

function get_detectionlimitparams(df)
    v1 = get_var(df, "samples")
    v2 = get_var(df, "cdfflag")
    v3 = get_var(df, "threshold")
    v4 = get_var(df, "numberofpixels")
    v5 = get_var(df, "plumelength")    
    return Datacontainer.DetectionLimitVars(v1, v2, v3, v4, v5)
end


################################################################################


"""
    getparameters(filename)

TBW
"""
function getparameters(filename)
    df = TOML.parsefile(filename)
    # SimulationParameters
    simparams = general_simulation_parameters(df)
    # specific simparams for LSQ, CFM or Detection limit
    if "LSQ" in keys(df)
        simvar = get_lsqparams(df["LSQ"])
    elseif "CFM" in keys(df)
        simvar = get_plumeparams(df["CFM"])
    elseif "DetectionLimit" in keys(df)
        simvar = get_detectionlimitparams(df["DetectionLimit"])
    end
    return simparams, simvar
end


end # end module
