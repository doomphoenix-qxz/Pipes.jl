#pipes.jl
include("substances.jl")
module pipes
import NLsolve
import subst
struct pipe{T<:Real}
    length::T
    diameter::T
    radius::T
    flowarea::T
    thickness::T
    outer_diameter::T
    outer_radius::T
    totalarea::T
    roughness::T
    elevation_change::T
    therm_cond::Function
end

valid_pipe_diams = [0.125, 0.25, 0.375, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5,
                    3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0]

function pipe(l, d, t, ϵ, z, kfunc)
    # Assumes a circular pipe of diameter d
    return pipe(l,d,d/2,pi*(d/2)^2,t,d+2t,(d+2t)/2,pi*((d+2t)/2)^2,ϵ,z,kfunc)
end

mutable struct pipeProps{T<:Real}
    dᵢ::T
    dₒ::T
    thickness::T
    rᵢ::T
    rₒ::T
    units::String
end

function pipeProps(d, t, units="in")
  di = d - 2t
  ro = d/2
  ri = di/2
  return pipeProps(di, d, t, ri, ro, units)
end

conversions = Dict(["cm","m"] => 1e-2,
                   ["m","cm"] => 100,
                   ["in","cm"] => 2.54,
                   ["cm","in"] => 1/2.54,
                   ["in","m"] => 2.54e-2,
                   ["m","in"] => 1/2.54e-2,
                   ["m","m"] => 1.0,
                   ["cm","cm"] => 1.0,
                   ["in","in"] => 1.0)

function convert(p::pipeProps, newUnits::String)
    key_ = [p.units, newUnits]
    conversion_factor = conversions[key_]
    p.dᵢ *= conversion_factor
    p.dₒ *= conversion_factor
    p.rᵢ *= conversion_factor
    p.rₒ *= conversion_factor
    p.units = newUnits
end

schedule40_props = [pipeProps(0.405, 0.068), pipeProps(0.540, 0.088),
                    pipeProps(0.675, 0.091), pipeProps(0.840, 0.109),
                    pipeProps(1.050, 0.113), pipeProps(1.315, 0.133),
                    pipeProps(1.660, 0.140), pipeProps(1.900, 0.145),
                    pipeProps(2.375, 0.154), pipeProps(2.875, 0.203),
                    pipeProps(3.500, 0.216), pipeProps(4.000, 0.226),
                    pipeProps(4.500, 0.237), pipeProps(5.563, 0.258),
                    pipeProps(6.625, 0.280), pipeProps(8.625, 0.322),
                    pipeProps(10.750, 0.365), pipeProps(12.750, 0.406),
                    pipeProps(16.000, 0.500)]
# Note: These outer diameters and wall thicknesses were obtained at this web
# site: http://www.engineeringtoolbox.com/steel-pipes-dimensions-d_43.html

for props in schedule40_props
    convert(props, "m")
end

props40_by_nominal_diameter = Dict(zip(valid_pipe_diams, schedule40_props))

schedule_props = Dict("40" => props40_by_nominal_diameter)

pipe_materials = Dict("Stainless Steel" => subst.ss_mat,
                      "Copper"          => subst.cu_mat)

function pipe(schedule_info::Tuple, mat::subst.material, l, Δz)
  # schedule_info is a Tuple: the schedule number comes first and the nominal
  # diameter comes second
  props = schedule_props[schedule_info[1]][schedule_info[2]]
  return pipe(l,props.dᵢ,props.dᵢ/2,pi*(props.rᵢ)^2,props.thickness,props.dₒ,
            props.rₒ,pi*(props.rₒ)^2,mat.roughness,Δz,mat.thermal_conductivity)
end

function pipe(schedule_info::Tuple, mat::String, l, Δz)
  return pipe(schedule_info, pipe_materials[mat], l, Δz)
end

function colebrook(Re, p::pipes.pipe, f_g)
    if f_g < 0
        return f_g * 1e5
    else
        return 10^(-1 / (2 * √(f_g))) - (p.roughness / (3.7 * p.diameter) + 2.51 / (Re * √(f_g)))
    end
end

function getf(Re, p::pipes.pipe)
    function dummysolve!(f_guess, fvec)
        fvec[1] = colebrook(Re, p, f_guess[1])
    end
    return NLsolve.nlsolve(dummysolve!, [0.001]).zero[1]
end

tube = pipe(("40", 1.0), subst.ss_mat, 0.8, 0.0)

end # ends pipes module
