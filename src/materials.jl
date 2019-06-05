"""
Struct material
Used for solid materials in pipes. So far includes just name, thermal
conductivity, and surface roughness.
"""
struct material
    name::String
    thermal_conductivity::Function
    roughness::Real
end

function stainless_steel_thermal_conductivity(T)
    # Source: http://www.mace.manchester.ac.uk/project/research/structures/strucfire/materialInFire/Steel/StainlessSteel/thermalProperties.htm
    return 14.6 + 1.27e-2 *(T-273.15)
end

function copper_thermal_conductivity(T)
    # Source: http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/cu.html
    return 14.6 + 1.27e-2 *(T-273.15)
end

ss_roughness_by_finish = Dict("2D" => 1E-6,
                              "2B" => 0.5E-6,
                              "2R" => 0.2E-6,
                              "BA" => 0.2E-6,
                              "2BB" => 0.1E-6)
# These surface roughness values are in m and were obtained from the following
# corporate site: http://www.outokumpu.com/en/products-properties/
# more-stainless/stainless-steel-surface-finishes/cold-rolled-finishes/Pages/default.aspx
# Note that the given values are the maximum roughness values in the range on
# the website, except for finish 2BB, which didn't have a roughness value
# so I took the lowest value for finish 2B since the compnay says this:
# "Due to the fact that surface roughness of the finish 2BB is lower than that
# of 2B, some modifications on lubrication during forming might be needed."

ss_mat = material("Stainless Steel", stainless_steel_thermal_conductivity,
                  ss_roughness_by_finish["2B"])

cu_roughness = 0.03e-3
# Copper roughness.
# Source: http://www.pressure-drop.com/Online-Calculator/rauh.html

cu_mat = material("Copper", copper_thermal_conductivity, cu_roughness)
