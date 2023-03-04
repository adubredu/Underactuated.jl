module Underactuated

using MuJoCo
using MuJoCo.PythonCall
using LinearAlgebra
using ControlSystems 
using JuMP 
using OSQP

include("types.jl")
include("step.jl")

include("robots/simple_pendulum.jl")
include("robots/acrobot.jl")
 

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

 
end
