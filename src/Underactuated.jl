module Underactuated

using MuJoCo
using MuJoCo.PythonCall
using LinearAlgebra

include("types.jl")
include("step.jl")

include("robots/simple_pendulum.jl")
include("robots/acrobot.jl")

include("controllers/energy_shaping.jl")

# Robot types
export SimplePendulumType,
        AcrobotType

# Robots 
export SimplePendulum,
        Acrobot

export create_robot,
       step,
       get_generalized_coordinates,
       compute_energy,
       apply_torque!,
       get_linearized_dynamics

# Controllers 
export energy_shaping_controller

end
