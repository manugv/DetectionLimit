module Constants

avogadro::Float64 = 6.0221408e+23
gravity::Float64 = 9.81
surfacepressure::Float64 = 101325.0

struct gasconst
    name::String
    molarmass::Float64
    background::Float64 # ppm
    ppm_to_gms::Float64
    function gasconst(name, molarmass, background, totairmass)
        ppm = totairmass * molarmass / 1e6
        return new(name, molarmass, background, ppm)
    end
end

struct airconst
    name::String
    molarmass::Float64 # g/mol
    column_mass_moles::Float64 # ((101325.0 /(9.81*28.9647)) * 1000) in moles
    function airconst(name, molarmass)
        tmp = (1000 * surfacepressure / gravity) / molarmass
        return new(name, molarmass, tmp)
    end
end


# Define different gases and constants
air = airconst("air", 28.9647)
co2 = gasconst("co2", 44.01, 412, air.column_mass_moles)
ch4 = gasconst("ch4", 16.04, 1.8957, air.column_mass_moles)
no2 = gasconst("no2", 46.01, 0, air.column_mass_moles)

gmspersec_2_MTperyear = 365 * 24 * 3600 / 1e12

end  # end module
